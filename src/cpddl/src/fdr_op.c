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

#include <boruvka/alloc.h>
#include "pddl/fdr_op.h"
#include "assert.h"

static void condEffFree(pddl_fdr_op_cond_eff_t *ce)
{
    pddlFDRPartStateFree(&ce->pre);
    pddlFDRPartStateFree(&ce->eff);
}

pddl_fdr_op_t *pddlFDROpNewEmpty(void)
{
    pddl_fdr_op_t *op = BOR_ALLOC(pddl_fdr_op_t);
    bzero(op, sizeof(*op));
    return op;
}

void pddlFDROpDel(pddl_fdr_op_t *op)
{
    if (op->name != NULL)
        BOR_FREE(op->name);
    pddlFDRPartStateFree(&op->pre);
    pddlFDRPartStateFree(&op->eff);
    for (int cei = 0; cei < op->cond_eff_size; ++cei)
        condEffFree(op->cond_eff + cei);
    if (op->cond_eff != NULL)
        BOR_FREE(op->cond_eff);
    BOR_FREE(op);
}

pddl_fdr_op_cond_eff_t *pddlFDROpAddEmptyCondEff(pddl_fdr_op_t *op)
{
    if (op->cond_eff_size >= op->cond_eff_alloc){
        if (op->cond_eff_alloc == 0)
            op->cond_eff_alloc = 1;
        op->cond_eff_alloc *= 2;
        op->cond_eff = BOR_REALLOC_ARR(op->cond_eff, pddl_fdr_op_cond_eff_t,
                                       op->cond_eff_alloc);
    }

    pddl_fdr_op_cond_eff_t *ce = op->cond_eff + op->cond_eff_size++;
    bzero(ce, sizeof(*ce));
    return ce;
}

void pddlFDROpRemapFacts(pddl_fdr_op_t *op, const pddl_fdr_vars_remap_t *rmp)
{
    pddlFDRPartStateRemapFacts(&op->pre, rmp);
    pddlFDRPartStateRemapFacts(&op->eff, rmp);

    int ins = 0;
    for (int cei = 0; cei < op->cond_eff_size; ++cei){
        pddl_fdr_op_cond_eff_t *ce = op->cond_eff + cei;
        pddlFDRPartStateRemapFacts(&ce->pre, rmp);
        pddlFDRPartStateRemapFacts(&ce->eff, rmp);
        if (ce->pre.fact_size == 0){
            for (int fi = 0; fi < ce->eff.fact_size; ++fi){
                const pddl_fdr_fact_t *f = ce->eff.fact + fi;
                ASSERT_RUNTIME(!pddlFDRPartStateIsSet(&ce->eff, f->var));
                if (!pddlFDRPartStateIsSet(&ce->pre, f->var)
                        || pddlFDRPartStateGet(&ce->pre, f->var) != f->val){
                    pddlFDRPartStateSet(&ce->eff, f->var, f->val);
                }
            }
            condEffFree(ce);

        }else{
            op->cond_eff[ins++] = *ce;
        }
    }
    op->cond_eff_size = ins;
}

void pddlFDROpApplyOnState(const pddl_fdr_op_t *op, int *state)
{
    for (int fi = 0; fi < op->eff.fact_size; ++fi)
        state[op->eff.fact[fi].var] = op->eff.fact[fi].val;
}

void pddlFDROpApplyOnState2(const pddl_fdr_op_t *op,
                            int num_vars,
                            const int *in_state,
                            int *out_state)
{
    memcpy(out_state, in_state, sizeof(int) * num_vars);
    pddlFDROpApplyOnState(op, out_state);
}


void pddlFDROpsInit(pddl_fdr_ops_t *ops)
{
    bzero(ops, sizeof(*ops));
}

void pddlFDROpsFree(pddl_fdr_ops_t *ops)
{
    for (int op_id = 0; op_id < ops->op_size; ++op_id){
        if (ops->op[op_id] != NULL)
            pddlFDROpDel(ops->op[op_id]);
    }
    if (ops->op != NULL)
        BOR_FREE(ops->op);
}

void pddlFDROpsDelSet(pddl_fdr_ops_t *ops, const bor_iset_t *set)
{
    int size = borISetSize(set);
    int cur = 0;
    int ins = 0;
    for (int op_id = 0; op_id < ops->op_size; ++op_id){
        if (cur < size && borISetGet(set, cur) == op_id){
            pddlFDROpDel(ops->op[op_id]);
            ++cur;
        }else{
            ops->op[op_id]->id = ins;
            ops->op[ins++] = ops->op[op_id];
        }
    }

    ops->op_size = ins;
}

void pddlFDROpsRemapFacts(pddl_fdr_ops_t *ops, const pddl_fdr_vars_remap_t *r)
{
    for (int op_id = 0; op_id < ops->op_size; ++op_id)
        pddlFDROpRemapFacts(ops->op[op_id], r);
}

void pddlFDROpsAddSteal(pddl_fdr_ops_t *ops, pddl_fdr_op_t *op)
{
    if (ops->op_size >= ops->op_alloc){
        if (ops->op_alloc == 0)
            ops->op_alloc = 8;
        ops->op_alloc *= 2;
        ops->op = BOR_REALLOC_ARR(ops->op, pddl_fdr_op_t *, ops->op_alloc);
    }

    op->id = ops->op_size;
    ops->op[ops->op_size++] = op;
}
