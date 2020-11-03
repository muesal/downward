/***
 * cpddl
 * -------
 * Copyright (c)2019 Daniel Fiser <danfis@danfis.cz>,
 * Faculty of Electrical Engineering, Czech Technical University in Prague.
 * All rights reserved.
 *
 * This file is part of cpddl.
 *
 * Distributed under the OSI-approved BSD License (the "License");
 * see accompanying file BDS-LICENSE for details or see
 * <http://www.opensource.org/licenses/bsd-license.php>.
 *
 * This software is distributed WITHOUT ANY WARRANTY; without even the
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the License for more information.
 */

#include <boruvka/sort.h>
#include "pddl/fdr.h"
#include "assert.h"

static void stripsToFDRState(const pddl_fdr_vars_t *fdr_var,
                             const bor_iset_t *state,
                             int *fdr_state);
static int stripsToFDRPartState(const pddl_fdr_vars_t *fdr_var,
                                const bor_iset_t *part_state,
                                pddl_fdr_part_state_t *fdr_ps);
static void addOp(pddl_fdr_ops_t *fdr_ops,
                  const pddl_fdr_vars_t *fdr_var,
                  const pddl_strips_t *strips,
                  const pddl_mutex_pairs_t *mutex,
                  int op_id);
static void printOp(const pddl_fdr_op_t *op, FILE *fout);

int pddlFDRInitFromStrips(pddl_fdr_t *fdr,
                          const pddl_strips_t *strips,
                          const pddl_mgroups_t *mg,
                          const pddl_mutex_pairs_t *mutex,
                          unsigned fdr_var_flags,
                          bor_err_t *err)
{
    bzero(fdr, sizeof(*fdr));

    BOR_INFO2(err, "Translation to FDR...");

    // variables
    if (pddlFDRVarsInitFromStrips(&fdr->var, strips, mg, mutex,
                                  fdr_var_flags) != 0){
        return -1;
    }
    BOR_INFO(err, "  Created %d variables.", fdr->var.var_size);
    int num_none_of_those = 0;
    for (int vi = 0; vi < fdr->var.var_size; ++vi){
        if (fdr->var.var[vi].val_none_of_those != -1)
            ++num_none_of_those;
    }
    BOR_INFO(err, "  Created %d none-of-those values.", num_none_of_those);

    fdr->goal_is_unreachable = strips->goal_is_unreachable;

    // Initial state
    fdr->init = BOR_ALLOC_ARR(int, fdr->var.var_size);
    stripsToFDRState(&fdr->var, &strips->init, fdr->init);

    // Goal
    pddlFDRPartStateInit(&fdr->goal);
    stripsToFDRPartState(&fdr->var, &strips->goal, &fdr->goal);

    // Operators
    pddlFDROpsInit(&fdr->op);
    for (int op_id = 0; op_id < strips->op.op_size; ++op_id)
        addOp(&fdr->op, &fdr->var, strips, mutex, op_id);

    return 0;
}

void pddlFDRFree(pddl_fdr_t *fdr)
{
    if (fdr->init != NULL)
        BOR_FREE(fdr->init);
    pddlFDRPartStateFree(&fdr->goal);
    pddlFDROpsFree(&fdr->op);
    pddlFDRVarsFree(&fdr->var);
}

void pddlFDRReduce(pddl_fdr_t *fdr,
                   const bor_iset_t *del_vars,
                   const bor_iset_t *_del_facts,
                   const bor_iset_t *del_ops)
{
    if (del_ops != NULL && borISetSize(del_ops) > 0)
        pddlFDROpsDelSet(&fdr->op, del_ops);

    BOR_ISET(del_facts);
    if (_del_facts != NULL && borISetSize(_del_facts) > 0)
        borISetUnion(&del_facts, _del_facts);

    if (del_vars != NULL && borISetSize(del_vars) > 0){
        int var;
        BOR_ISET_FOR_EACH(del_vars, var){
            for (int val = 0; val < fdr->var.var[var].val_size; ++val)
                borISetAdd(&del_facts, fdr->var.var[var].val[val].global_id);
        }
    }

    if (borISetSize(&del_facts) > 0){
        int old_var_size = fdr->var.var_size;

        pddl_fdr_vars_remap_t remap;
        // Delete facts
        pddlFDRVarsDelFacts(&fdr->var, &del_facts, &remap);
        // Remap facts in operators
        pddlFDROpsRemapFacts(&fdr->op, &remap);

        // Remap the initial state
        for (int v = 0; v < old_var_size; ++v){
            if (remap.remap[v][fdr->init[v]] != NULL){
                const pddl_fdr_val_t *val = remap.remap[v][fdr->init[v]];
                fdr->init[val->var_id] = val->val_id;
            }
        }

        // Remap goal
        pddlFDRPartStateRemapFacts(&fdr->goal, &remap);

        // Remove operators with empty effects
        BOR_ISET(useless_ops);
        for (int op_id = 0; op_id < fdr->op.op_size; ++op_id){
            const pddl_fdr_op_t *op = fdr->op.op[op_id];
            if (op->eff.fact_size == 0 && op->cond_eff_size == 0)
                borISetAdd(&useless_ops, op_id);
        }
        if (borISetSize(&useless_ops) > 0)
            pddlFDROpsDelSet(&fdr->op, &useless_ops);
        borISetFree(&useless_ops);

        // Set cond-eff flag
        fdr->has_cond_eff = 0;
        for (int op_id = 0; op_id < fdr->op.op_size; ++op_id){
            if (fdr->op.op[op_id]->cond_eff_size > 0){
                fdr->has_cond_eff = 1;
                break;
            }
        }

        pddlFDRVarsRemapFree(&remap);
    }
    borISetFree(&del_facts);
}

void pddlFDRPrintFD(const pddl_fdr_t *fdr,
                    const pddl_mgroups_t *mg,
                    FILE *fout)
{
    fprintf(fout, "begin_version\n3\nend_version\n");
    fprintf(fout, "begin_metric\n1\nend_metric\n");

    // variables
    fprintf(fout, "%d\n", fdr->var.var_size);
    for (int vi = 0; vi < fdr->var.var_size; ++vi){
        const pddl_fdr_var_t *var = fdr->var.var + vi;
        fprintf(fout, "begin_variable\n");
        fprintf(fout, "var%d\n", vi);
        fprintf(fout, "-1\n");
        fprintf(fout, "%d\n", var->val_size);
        for (int vali = 0; vali < var->val_size; ++vali)
            fprintf(fout, "%s\n", var->val[vali].name);
        fprintf(fout, "end_variable\n");
    }

    // mutex groups
    if (mg == NULL){
        fprintf(fout, "0\n");
    }else{
        fprintf(fout, "%d\n", mg->mgroup_size);
        for (int mi = 0; mi < mg->mgroup_size; ++mi){
            const pddl_mgroup_t *m = mg->mgroup + mi;
            fprintf(fout, "begin_mutex_group\n");
            fprintf(fout, "%d\n", borISetSize(&m->mgroup));
            int fact_id;
            BOR_ISET_FOR_EACH(&m->mgroup, fact_id){
                // TODO
                int val_id = borISetGet(&fdr->var.strips_id_to_val[fact_id], 0);
                const pddl_fdr_val_t *v = fdr->var.global_id_to_val[val_id];
                fprintf(fout, "%d %d\n", v->var_id, v->val_id);
            }
            fprintf(fout, "end_mutex_group\n");
        }
    }

    // initial state
    fprintf(fout, "begin_state\n");
    for (int vi = 0; vi < fdr->var.var_size; ++vi)
        fprintf(fout, "%d\n", fdr->init[vi]);
    fprintf(fout, "end_state\n");

    // goal
    fprintf(fout, "begin_goal\n");
    fprintf(fout, "%d\n", fdr->goal.fact_size);
    for (int i = 0; i < fdr->goal.fact_size; ++i){
        const pddl_fdr_fact_t *f = fdr->goal.fact + i;
        fprintf(fout, "%d %d\n", f->var, f->val);
    }
    fprintf(fout, "end_goal\n");

    // operators
    fprintf(fout, "%d\n", fdr->op.op_size);
    for (int op_id = 0; op_id < fdr->op.op_size; ++op_id)
        printOp(fdr->op.op[op_id], fout);

    // axioms
    fprintf(fout, "0\n");
}

static void stripsToFDRState(const pddl_fdr_vars_t *fdr_var,
                             const bor_iset_t *state,
                             int *fdr_state)
{
    for (int vi = 0; vi < fdr_var->var_size; ++vi)
        fdr_state[vi] = -1;

    int fact_id;
    BOR_ISET_FOR_EACH(state, fact_id){
        int val_id;
        BOR_ISET_FOR_EACH(&fdr_var->strips_id_to_val[fact_id], val_id){
            const pddl_fdr_val_t *v = fdr_var->global_id_to_val[val_id];
            fdr_state[v->var_id] = v->val_id;
        }
    }

    for (int vi = 0; vi < fdr_var->var_size; ++vi){
        if (fdr_state[vi] == -1){
            ASSERT(fdr_var->var[vi].val_none_of_those >= 0);
            fdr_state[vi] = fdr_var->var[vi].val_none_of_those;
        }
    }
}

static void setDelEffFact(const pddl_mutex_pairs_t *mutex,
                          const pddl_fdr_vars_t *fdr_var,
                          const bor_iset_t *pre,
                          const bor_iset_t *ce_pre,
                          pddl_fdr_part_state_t *eff,
                          int fact_id,
                          const pddl_fdr_val_t *v)
{
    const pddl_fdr_var_t *var = fdr_var->var + v->var_id;
    if (!pddlMutexPairsIsMutexFactSet(mutex, fact_id, pre)
            && (ce_pre == NULL
                    || !pddlMutexPairsIsMutexFactSet(mutex, fact_id, ce_pre))
            && var->val_none_of_those >= 0){
        pddlFDRPartStateSet(eff, var->var_id, var->val_none_of_those);
    }
}

static void stripsToFDRDelEff(const pddl_fdr_vars_t *fdr_var,
                              const bor_iset_t *eff,
                              pddl_fdr_part_state_t *fdr_eff,
                              const pddl_mutex_pairs_t *mutex,
                              const bor_iset_t *pre,
                              const bor_iset_t *ce_pre)
{
    int fact_id;
    BOR_ISET_FOR_EACH(eff, fact_id){
        int val_id;
        BOR_ISET_FOR_EACH(&fdr_var->strips_id_to_val[fact_id], val_id){
            const pddl_fdr_val_t *v = fdr_var->global_id_to_val[val_id];
            setDelEffFact(mutex, fdr_var, pre, ce_pre, fdr_eff, fact_id, v);
        }
    }
}

static int stripsToFDRPartState(const pddl_fdr_vars_t *fdr_var,
                                const bor_iset_t *part_state,
                                pddl_fdr_part_state_t *fdr_ps)
{
    int ret = 0;
    int fact_id;
    BOR_ISET_FOR_EACH(part_state, fact_id){
        int val_id;
        BOR_ISET_FOR_EACH(&fdr_var->strips_id_to_val[fact_id], val_id){
            const pddl_fdr_val_t *v = fdr_var->global_id_to_val[val_id];
            if (pddlFDRPartStateIsSet(fdr_ps, v->var_id))
                ret = -1;
            pddlFDRPartStateSet(fdr_ps, v->var_id, v->val_id);
        }
    }

    return ret;
}

static int cmpCondEff(const void *a, const void *b, void *_)
{
    const pddl_fdr_op_cond_eff_t *ce1 = a;
    const pddl_fdr_op_cond_eff_t *ce2 = b;
    return pddlFDRPartStateCmp(&ce1->pre, &ce2->pre);
}

static void addOp(pddl_fdr_ops_t *fdr_ops,
                  const pddl_fdr_vars_t *fdr_var,
                  const pddl_strips_t *strips,
                  const pddl_mutex_pairs_t *mutex,
                  int op_id)
{
    const pddl_strips_op_t *op = strips->op.op[op_id];
    pddl_fdr_op_t *fdr_op = pddlFDROpNewEmpty();
    pddl_fdr_part_state_t pre;

    if (op->name != NULL)
        fdr_op->name = BOR_STRDUP(op->name);
    fdr_op->cost = op->cost;

    pddlFDRPartStateInit(&pre);
    if (stripsToFDRPartState(fdr_var, &op->pre, &pre) == 0){
        pddlFDRPartStateInitCopy(&fdr_op->pre, &pre);
    }else{
        pddlFDRPartStateFree(&pre);
        pddlFDROpDel(fdr_op);
        return;
    }
    pddlFDRPartStateFree(&pre);

    stripsToFDRPartState(fdr_var, &op->pre, &fdr_op->pre);
    stripsToFDRDelEff(fdr_var, &op->del_eff, &fdr_op->eff,
                      mutex, &op->pre, NULL);
    stripsToFDRPartState(fdr_var, &op->add_eff, &fdr_op->eff);

    for (int cei = 0; cei < op->cond_eff_size; ++cei){
        const pddl_strips_op_cond_eff_t *ce = op->cond_eff + cei;
        pddl_fdr_op_cond_eff_t *fdr_ce = pddlFDROpAddEmptyCondEff(fdr_op);
        pddlFDRPartStateInit(&pre);
        if (stripsToFDRPartState(fdr_var, &ce->pre, &pre) == 0){
            pddlFDRPartStateInitCopy(&fdr_ce->pre, &pre);
        }else{
            pddlFDRPartStateFree(&pre);
            pddlFDROpDel(fdr_op);
            return;
        }
        pddlFDRPartStateFree(&pre);
        stripsToFDRDelEff(fdr_var, &ce->del_eff, &fdr_ce->eff,
                          mutex, &op->pre, &ce->pre);
        stripsToFDRPartState(fdr_var, &ce->add_eff, &fdr_ce->eff);
    }

    if (fdr_op->cond_eff_size > 1){
        borSort(fdr_op->cond_eff, fdr_op->cond_eff_size,
                sizeof(pddl_fdr_op_cond_eff_t), cmpCondEff, NULL);
    }

    pddlFDROpsAddSteal(fdr_ops, fdr_op);
}

static void printOp(const pddl_fdr_op_t *op, FILE *fout)
{
    pddl_fdr_part_state_t prevail;

    fprintf(fout, "begin_operator\n");
    fprintf(fout, "%s\n", op->name);

    BOR_ISET(eff_var);
    for (int i = 0; i < op->eff.fact_size; ++i)
        borISetAdd(&eff_var, op->eff.fact[i].var);
    for (int cei = 0; cei < op->cond_eff_size; ++cei){
        const pddl_fdr_op_cond_eff_t *ce = op->cond_eff + cei;
        for (int i = 0; i < ce->eff.fact_size; ++i)
            borISetAdd(&eff_var, ce->eff.fact[i].var);
    }

    pddlFDRPartStateInit(&prevail);
    for (int pi = 0; pi < op->pre.fact_size; ++pi){
        const pddl_fdr_fact_t *f = op->pre.fact + pi;
        if (!borISetIn(f->var, &eff_var))
            pddlFDRPartStateSet(&prevail, f->var, f->val);
    }

    fprintf(fout, "%d\n", prevail.fact_size);
    for (int i = 0; i < prevail.fact_size; ++i)
        fprintf(fout, "%d %d\n", prevail.fact[i].var, prevail.fact[i].val);
    pddlFDRPartStateFree(&prevail);
    borISetFree(&eff_var);

    int num_effs = op->eff.fact_size;
    for (int cei = 0; cei < op->cond_eff_size; ++cei)
        num_effs += op->cond_eff[cei].eff.fact_size;

    fprintf(fout, "%d\n", num_effs);
    for (int i = 0; i < op->eff.fact_size; ++i){
        const pddl_fdr_fact_t *f = op->eff.fact + i;
        int pre = pddlFDRPartStateGet(&op->pre, f->var);
        fprintf(fout, "0 %d %d %d\n", f->var, pre, f->val);
    }

    for (int cei = 0; cei < op->cond_eff_size; ++cei){
        const pddl_fdr_op_cond_eff_t *ce = op->cond_eff + cei;
        for (int i = 0; i < ce->eff.fact_size; ++i){
            const pddl_fdr_fact_t *f = ce->eff.fact + i;
            int num_prevails = ce->pre.fact_size;
            fprintf(fout, "%d", num_prevails);
            for (int pi = 0; pi < ce->pre.fact_size; ++pi){
                const pddl_fdr_fact_t *p = ce->pre.fact + pi;
                fprintf(fout, " %d %d", p->var, p->val);
            }
            int pre = pddlFDRPartStateGet(&op->pre, f->var);
            fprintf(fout, " %d %d %d\n", f->var, pre, f->val);
        }
    }

    fprintf(fout, "%d\n", op->cost);
    fprintf(fout, "end_operator\n");
}
