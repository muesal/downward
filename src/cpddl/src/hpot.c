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

#include <boruvka/hashset.h>
#include <boruvka/rand.h>
#include "pddl/hpot.h"
#include "pddl/pot.h"
#include "pddl/critical_path.h"
#include "pddl/random_walk.h"
#include "assert.h"

#define ROUND_EPS 0.001

static const uint32_t rand_sampler_seed = 524287;
static const uint32_t rand_diverse_seed = 131071;


static void init(pddl_hpot_t *hpot, int var_size)
{
    bzero(hpot, sizeof(*hpot));
    hpot->pot_alloc = hpot->pot_size = 0;
    hpot->pot = NULL;
    hpot->var_size = var_size;
    hpot->func = BOR_ALLOC_ARR(double, var_size);
}

static void addFunc2(pddl_hpot_t *hpot, const double *p)
{
    if (hpot->pot_size == hpot->pot_alloc){
        int old_size = hpot->pot_alloc;
        if (hpot->pot_alloc == 0)
            hpot->pot_alloc = 2;
        hpot->pot_alloc *= 2;
        hpot->pot = BOR_REALLOC_ARR(hpot->pot, double *, hpot->pot_alloc);
        for (int i = old_size; i < hpot->pot_alloc; ++i)
            hpot->pot[i] = BOR_ALLOC_ARR(double, hpot->var_size);
    }
    double *dst = hpot->pot[hpot->pot_size++];
    memcpy(dst, p, sizeof(double) * hpot->var_size);
}

static void addFunc(pddl_hpot_t *hpot)
{
    addFunc2(hpot, hpot->func);
}

static int solve2(pddl_hpot_t *hpot, pddl_pot_t *pot, double *w)
{
    return pddlPotSolve(pot, w, hpot->var_size, 0);
}

static int solve(pddl_hpot_t *hpot, pddl_pot_t *pot)
{
    return solve2(hpot, pot, hpot->func);
}


static int roundOff(double z)
{
    return ceil(z - ROUND_EPS);
}

static double fdrStateEstimateDbl(const double *pot,
                                  const pddl_fdr_vars_t *vars,
                                  const int *state)
{
    double p = 0;
    for (int var = 0; var < vars->var_size; ++var)
        p += pot[vars->var[var].val[state[var]].global_id];
    if (p < 0.)
        return 0;
    if (p > 1E8)
        return PDDL_COST_DEAD_END;
    return p;
}

static int fdrStateEstimate(const double *pot,
                            const pddl_fdr_vars_t *vars,
                            const int *state)
{
    return roundOff(fdrStateEstimateDbl(pot, vars, state));
}

static void initPot(pddl_hpot_t *hpot,
                    pddl_pot_t *pot,
                    const pddl_fdr_t *fdr,
                    const pddl_mg_strips_t *mg_strips,
                    const pddl_mutex_pairs_t *mutex,
                    const pddl_hpot_config_t *cfg,
                    bor_err_t *err)
{

    if (cfg->weak_disambiguation){
        pddlPotInitMGStripsSingleFactDisamb(pot, mg_strips, mutex);
        BOR_INFO(err, "Pot: Initialized with weak-disambiguation."
                      " vars: %d, op-constr: %d,"
                      " goal-constr: %d, maxpots: %d",
                      pot->var_size,
                      pot->constr_op.size,
                      pot->constr_goal.size,
                      pot->maxpot_size);

    }else if (cfg->disambiguation){
        pddlPotInitMGStrips(pot, mg_strips, mutex);
        BOR_INFO(err, "Pot: Initialized with disambiguation."
                      " vars: %d, op-constr: %d,"
                      " goal-constr: %d, maxpots: %d",
                      pot->var_size,
                      pot->constr_op.size,
                      pot->constr_goal.size,
                      pot->maxpot_size);

    }else{
        pddlPotInitFDR(pot, fdr);
        BOR_INFO(err, "Pot: Initialized without disambiguation."
                      " vars: %d, op-constr: %d,"
                      " goal-constr: %d, maxpots: %d",
                      pot->var_size,
                      pot->constr_op.size,
                      pot->constr_goal.size,
                      pot->maxpot_size);
    }
}

static int addInitConstr(pddl_hpot_t *hpot,
                         pddl_pot_t *pot,
                         const pddl_fdr_t *fdr,
                         const pddl_hpot_config_t *cfg,
                         bor_err_t *err)
{
    pddlPotResetLowerBoundConstr(pot);
    pddlPotSetObjFDRState(pot, &fdr->var, fdr->init);
    int ret = solve(hpot, pot);
    if (ret != 0){
        BOR_INFO2(err, "Pot: No optimal solution for the initial state");
        return ret;
    }

    double rhs = fdrStateEstimateDbl(hpot->func, &fdr->var, fdr->init);
    BOR_INFO(err, "Pot: Solved for the initial state: %.4f", rhs);
    // make sure it is feasible
    rhs = floor((rhs - ROUND_EPS) * 100.) / 100.;
    rhs *= cfg->init_constr_coef;

    BOR_ISET(vars);
    for (int var = 0; var < fdr->var.var_size; ++var){
        int v = fdr->var.var[var].val[fdr->init[var]].global_id;
        borISetAdd(&vars, v);
    }
    pddlPotSetLowerBoundConstr(pot, &vars, rhs);
    BOR_INFO(err, "Pot: added lower bound constraint with rhs: %.2f", rhs);
    borISetFree(&vars);

    return 0;
}




#define STATE_SAMPLER_SYNTACTIC 0
#define STATE_SAMPLER_SYNTACTIC_MUTEX 1
#define STATE_SAMPLER_RANDOM_WALK 2
struct state_sampler {
    int type;
    const pddl_fdr_t *fdr;
    const pddl_mutex_pairs_t *mutex;
    pddl_random_walk_t random_walk;
    int random_walk_max_steps;
    bor_rand_mt_t *rnd;
    int *state;
};
typedef struct state_sampler state_sampler_t;

static void stateSamplerInit(state_sampler_t *s,
                             const pddl_hpot_config_t *cfg,
                             const pddl_fdr_t *fdr,
                             const pddl_mutex_pairs_t *mutex,
                             pddl_hpot_t *hpot,
                             pddl_pot_t *pot,
                             bor_err_t *err)
{
    bzero(s, sizeof(*s));
    s->fdr = fdr;
    s->state = BOR_ALLOC_ARR(int, fdr->var.var_size);
    if (cfg->samples_random_walk){
        s->type = STATE_SAMPLER_RANDOM_WALK;
        //pddlRandomWalkInit(&s->random_walk, fdr, NULL);
        pddlRandomWalkInitSeed(&s->random_walk, fdr, NULL, rand_sampler_seed);

        pddlPotSetObjFDRState(pot, &fdr->var, fdr->init);
        int ret = solve(hpot, pot);
        if (ret != 0){
            BOR_INFO2(err, "Pot: No optimal solution for the initial state");
            s->random_walk_max_steps = 0;
        }

        int hinit = fdrStateEstimate(hpot->func, &fdr->var, fdr->init);
        double avg_op_cost = 0.;
        for (int oi = 0; oi < fdr->op.op_size; ++oi)
            avg_op_cost += fdr->op.op[oi]->cost;
        avg_op_cost /= fdr->op.op_size;
        if (avg_op_cost < 1E-2){
            s->random_walk_max_steps = 10;
        }else{
            s->random_walk_max_steps = (ceil(hinit / avg_op_cost) + .5) * 4;
        }

    }else{
        //s->rnd = borRandMTNewAuto();
        s->rnd = borRandMTNew(rand_sampler_seed);
        if (mutex != NULL){
            s->mutex = mutex;
            s->type = STATE_SAMPLER_SYNTACTIC_MUTEX;
        }else{
            s->type = STATE_SAMPLER_SYNTACTIC;
        }
    }
}

static void stateSamplerFree(state_sampler_t *s)
{
    if (s->type == STATE_SAMPLER_RANDOM_WALK)
        pddlRandomWalkFree(&s->random_walk);
    if (s->state != NULL)
        BOR_FREE(s->state);
    if (s->rnd != NULL)
        borRandMTDel(s->rnd);
}

static void stateSamplerSample(state_sampler_t *s, bor_err_t *err)
{
    if (s->type == STATE_SAMPLER_SYNTACTIC){
        for (int var = 0; var < s->fdr->var.var_size; ++var){
            int val = borRandMT(s->rnd, 0, s->fdr->var.var[var].val_size);
            val = BOR_MIN(val, s->fdr->var.var[var].val_size - 1);
            s->state[var] = val;
        }

    }else if (s->type == STATE_SAMPLER_SYNTACTIC_MUTEX){
        ASSERT(s->mutex != NULL);
        BOR_ISET(state);
        unsigned long count = 0UL;
        do {
            borISetEmpty(&state);
            for (int var = 0; var < s->fdr->var.var_size; ++var){
                int val = borRandMT(s->rnd, 0, s->fdr->var.var[var].val_size);
                val = BOR_MIN(val, s->fdr->var.var[var].val_size - 1);
                s->state[var] = val;
                borISetAdd(&state, s->fdr->var.var[var].val[val].global_id);
            }
            if (++count % 100000UL == 0UL)
                BOR_INFO(err, "Pot: tried %lu random states", count);
        } while (pddlMutexPairsIsMutexSet(s->mutex, &state));
        borISetFree(&state);

    }else if (s->type == STATE_SAMPLER_RANDOM_WALK){
        pddlRandomWalkSampleState(&s->random_walk,
                                  s->fdr->init,
                                  s->random_walk_max_steps,
                                  s->state);
    }
}


static double countStatesMutex(const pddl_mgroups_t *mgs,
                               const pddl_mutex_pairs_t *mutex,
                               const bor_iset_t *fixed)
{
    if (pddlMutexPairsIsMutexSet(mutex, fixed))
        return 0.;

    if (fixed == NULL || borISetSize(fixed) == 0){
        double num = borISetSize(&mgs->mgroup[0].mgroup);
        for (int i = 1; i < mgs->mgroup_size; ++i)
            num *= borISetSize(&mgs->mgroup[i].mgroup);
        return num;
    }

    double num = 1.;
    for (int mgi = 0; mgi < mgs->mgroup_size; ++mgi){
        int mg_size = 0;
        int fact;
        BOR_ISET_FOR_EACH(&mgs->mgroup[mgi].mgroup, fact){
            if (!pddlMutexPairsIsMutexFactSet(mutex, fact, fixed))
                mg_size += 1;
        }
        num *= (double)mg_size;
    }
    return num;
}

static void setObjAllStatesMutex1(pddl_pot_t *pot,
                                  const pddl_mgroups_t *mgs,
                                  const pddl_mutex_pairs_t *mutex)
{
    double *coef = BOR_CALLOC_ARR(double, pot->var_size);
    BOR_ISET(fixed);

    for (int mgi = 0; mgi < mgs->mgroup_size; ++mgi){
        const pddl_mgroup_t *mg = mgs->mgroup + mgi;
        double sum = 0.;
        int fixed_fact;
        BOR_ISET_FOR_EACH(&mg->mgroup, fixed_fact){
            borISetEmpty(&fixed);
            borISetAdd(&fixed, fixed_fact);
            coef[fixed_fact] = countStatesMutex(mgs, mutex, &fixed);
            sum += coef[fixed_fact];
        }
        BOR_ISET_FOR_EACH(&mg->mgroup, fixed_fact){
            coef[fixed_fact] /= sum;
            if (coef[fixed_fact] < 1E-6)
                coef[fixed_fact] = 0.;
        }
    }

    pddlPotSetObj(pot, coef);

    borISetFree(&fixed);
    if (coef != NULL)
        BOR_FREE(coef);
}

static void setObjAllStatesMutex2(pddl_pot_t *pot,
                                  const pddl_mgroups_t *mgs,
                                  int fact_size,
                                  const pddl_mutex_pairs_t *mutex)
{
    double *coef = BOR_CALLOC_ARR(double, pot->var_size);
    BOR_ISET(fixed);

    for (int mgi = 0; mgi < mgs->mgroup_size; ++mgi){
        const pddl_mgroup_t *mg = mgs->mgroup + mgi;
        double sum = 0.;
        int fixed_fact;
        BOR_ISET_FOR_EACH(&mg->mgroup, fixed_fact){
            coef[fixed_fact] = 0.;
            for (int f = 0; f < fact_size; ++f){
                if (f == fixed_fact)
                    continue;

                borISetEmpty(&fixed);
                borISetAdd(&fixed, fixed_fact);
                borISetAdd(&fixed, f);
                ASSERT(borISetSize(&fixed) == 2);
                coef[fixed_fact] += countStatesMutex(mgs, mutex, &fixed);
            }
            sum += coef[fixed_fact];
        }
        BOR_ISET_FOR_EACH(&mg->mgroup, fixed_fact){
            coef[fixed_fact] /= sum;
            if (coef[fixed_fact] < 1E-6)
                coef[fixed_fact] = 0.;
        }
    }

    pddlPotSetObj(pot, coef);

    borISetFree(&fixed);
    if (coef != NULL)
        BOR_FREE(coef);
}


static void setObjAllStatesMutex(pddl_pot_t *pot,
                                 const pddl_mg_strips_t *s,
                                 const pddl_mutex_pairs_t *mutex,
                                 int mutex_size)
{
    if (mutex_size == 1){
        setObjAllStatesMutex1(pot, &s->mg, mutex);
    }else if (mutex_size == 2){
        setObjAllStatesMutex2(pot, &s->mg, s->strips.fact.fact_size, mutex);
    }else{
        ASSERT_RUNTIME_M(0, "mutex-size >= 3 is not supported!");
    }
}

static void setObjAllStatesMutexConditioned(pddl_pot_t *pot,
                                            const bor_iset_t *cond,
                                            const pddl_mg_strips_t *s,
                                            const pddl_mutex_pairs_t *mutex,
                                            int mutex_size)
{
    pddl_mgroups_t mgs;
    pddlMGroupsInitEmpty(&mgs);
    BOR_ISET(mg);
    for (int mgi = 0; mgi < s->mg.mgroup_size; ++mgi){
        int fact_id;
        borISetEmpty(&mg);
        BOR_ISET_FOR_EACH(&s->mg.mgroup[mgi].mgroup, fact_id){
            if (!pddlMutexPairsIsMutexFactSet(mutex, fact_id, cond))
                borISetAdd(&mg, fact_id);
        }
        pddlMGroupsAdd(&mgs, &mg);
        ASSERT_RUNTIME(borISetSize(&mg) > 0);
    }
    borISetFree(&mg);

    if (mutex_size == 1){
        setObjAllStatesMutex1(pot, &mgs, mutex);
    }else if (mutex_size == 2){
        setObjAllStatesMutex2(pot, &mgs, s->strips.fact.fact_size, mutex);
    }else{
        ASSERT_RUNTIME_M(0, "mutex-size >= 3 is not supported!");
    }
    pddlMGroupsFree(&mgs);
}

static int allStatesMutexCond(pddl_hpot_t *hpot,
                              pddl_pot_t *pot,
                              const pddl_mg_strips_t *mg_strips,
                              const pddl_mutex_pairs_t *mutex,
                              int mutex_size,
                              const bor_iset_t *facts,
                              bor_err_t *err)
{
    BOR_ISET(cond);
    int fact_id;
    int count = 0;
    BOR_ISET_FOR_EACH(facts, fact_id){
        borISetEmpty(&cond);
        borISetAdd(&cond, fact_id);
        setObjAllStatesMutexConditioned(pot, &cond, mg_strips, mutex,
                                        mutex_size);
        if (solve(hpot, pot) == 0)
            addFunc(hpot);
        if (++count % 10 == 0){
            BOR_INFO(err, "Computed conditioned func %d/%d and generated %d"
                          " potential functions",
                     count, borISetSize(facts), hpot->pot_size);
        }
    }
    BOR_INFO(err, "Computed conditioned func %d/%d and generated %d"
                  " potential functions",
             count, borISetSize(facts), hpot->pot_size);
    borISetFree(&cond);

    if (hpot->pot_size > 0)
        return 0;
    return -1;
}

static int allStatesMutexCond2(pddl_hpot_t *hpot,
                               pddl_pot_t *pot,
                               const pddl_mg_strips_t *mg_strips,
                               const pddl_mutex_pairs_t *mutex,
                               int mutex_size,
                               int num_samples,
                               bor_err_t *err)
{
    bor_rand_mt_t *rnd = borRandMTNew(rand_sampler_seed);
    BOR_ISET(cond);
    int count = 0;
    int fact_size = mg_strips->strips.fact.fact_size;
    for (int i = 0; i < num_samples; ++i){
        int f1 = borRandMT(rnd, 0, fact_size);
        int f2 = borRandMT(rnd, 0, fact_size);
        if (pddlMutexPairsIsMutex(mutex, f1, f2))
            continue;
        borISetEmpty(&cond);
        borISetAdd(&cond, f1);
        borISetAdd(&cond, f2);

        setObjAllStatesMutexConditioned(pot, &cond, mg_strips, mutex,
                                        mutex_size);
        if (solve(hpot, pot) == 0)
            addFunc(hpot);
        if (++count % 10 == 0){
            BOR_INFO(err, "Computed conditioned func^2 %d and generated %d"
                          " potential functions",
                     count, hpot->pot_size);
        }
    }
    BOR_INFO(err, "Computed conditioned func^2 %d and generated %d"
                  " potential functions",
             count, hpot->pot_size);
    borISetFree(&cond);
    borRandMTDel(rnd);

    if (hpot->pot_size > 0)
        return 0;
    return -1;
}


static int samples(pddl_hpot_t *hpot,
                   pddl_pot_t *pot,
                   const pddl_fdr_t *fdr,
                   const pddl_mutex_pairs_t *mutex,
                   const pddl_hpot_config_t *cfg,
                   bor_err_t *err)
{
    BOR_INFO(err, "Pot: generating %d samples (mutex: %d, random-walk: %d)...",
             cfg->num_samples,
             (mutex != NULL),
             cfg->samples_random_walk);

    state_sampler_t sampler;
    stateSamplerInit(&sampler, cfg, fdr, mutex, hpot, pot, err);

    int num_states = 0;
    double *coef = BOR_CALLOC_ARR(double, pot->var_size);
    for (int si = 0; si < cfg->num_samples; ++si){
        stateSamplerSample(&sampler, err);
        if (cfg->obj == PDDL_HPOT_OBJ_SAMPLES_MAX){
            bzero(coef, sizeof(double) * pot->var_size);
            for (int var = 0; var < fdr->var.var_size; ++var)
                coef[fdr->var.var[var].val[sampler.state[var]].global_id] = 1.;

            pddlPotSetObj(pot, coef);
            // Dead-ends are simply skipped
            if (solve(hpot, pot) == 0){
                int h = fdrStateEstimate(hpot->func, &fdr->var, sampler.state);
                if (h != PDDL_COST_DEAD_END){
                    addFunc(hpot);
                    ++num_states;
                    if ((si + 1) % 100 == 0){
                        BOR_INFO(err, "Pot: Solved for state: %d/%d",
                                 num_states, cfg->num_samples);
                    }
                }
            }

        }else{
            for (int var = 0; var < fdr->var.var_size; ++var)
                coef[fdr->var.var[var].val[sampler.state[var]].global_id] += 1.;
        }
    }

    int ret = 0;
    if (cfg->obj == PDDL_HPOT_OBJ_SAMPLES_SUM){
        pddlPotSetObj(pot, coef);
        if (solve(hpot, pot) == 0){
            addFunc(hpot);
            BOR_INFO(err, "Pot: Solved for a sum of %d/%d states",
                     num_states, cfg->num_samples);
        }else{
            BOR_INFO(err, "Pot: No solution for sum of %d/%d states",
                     num_states, cfg->num_samples);
            ret = -1;
        }
    }

    if (coef != NULL)
        BOR_FREE(coef);
    stateSamplerFree(&sampler);

    return ret;
}

static void setStateToFDRState(const bor_iset_t *state,
                               int *fdr_state,
                               const pddl_fdr_t *fdr)
{
    int fact_id;
    BOR_ISET_FOR_EACH(state, fact_id){
        const pddl_fdr_val_t *v = fdr->var.global_id_to_val[fact_id];
        fdr_state[v->var_id] = v->val_id;
    }
}

struct diverse_pot {
    double *coef;
    double **func;
    double *avg_func;
    int *state_est;
    bor_hashset_t states;
    int active_states;
    bor_rand_mt_t *rnd;
};
typedef struct diverse_pot diverse_pot_t;

static void diverseInit(diverse_pot_t *div,
                        const pddl_pot_t *pot,
                        const pddl_fdr_t *fdr,
                        int num_samples)
{
    div->coef = BOR_ALLOC_ARR(double, pot->var_size);
    div->func = BOR_ALLOC_ARR(double *, num_samples);
    for (int i = 0; i < num_samples; ++i)
        div->func[i] = BOR_ALLOC_ARR(double, fdr->var.global_id_size);
    div->avg_func = BOR_ALLOC_ARR(double, fdr->var.global_id_size);
    div->state_est = BOR_CALLOC_ARR(int, num_samples);
    borHashSetInitISet(&div->states);
    div->active_states = 0;
    //div->rnd = borRandMTNewAuto();
    div->rnd = borRandMTNew(rand_diverse_seed);
}

static void diverseFree(diverse_pot_t *div,
                        const pddl_pot_t *pot,
                        const pddl_fdr_t *fdr,
                        int num_samples)
{
    BOR_FREE(div->coef);
    for (int i = 0; i < num_samples; ++i)
        BOR_FREE(div->func[i]);
    BOR_FREE(div->func);
    BOR_FREE(div->avg_func);
    BOR_FREE(div->state_est);
    borHashSetFree(&div->states);
    borRandMTDel(div->rnd);
}

static void diverseGenStates(diverse_pot_t *div,
                             pddl_hpot_t *hpot,
                             pddl_pot_t *pot,
                             const pddl_fdr_t *fdr,
                             const pddl_hpot_config_t *_cfg,
                             bor_err_t *err)
{
    pddl_hpot_config_t cfg = *_cfg;
    // force random walk
    cfg.samples_random_walk = 1;
    ASSERT_RUNTIME(cfg.num_samples > 0);

    BOR_INFO(err, "Pot: generating %d samples with random walk and"
                  " computing potentials...", cfg.num_samples);
    state_sampler_t sampler;
    stateSamplerInit(&sampler, &cfg, fdr, NULL, hpot, pot, err);

    // Samples states, filter out dead-ends and compute estimate for each
    // state
    int num_states = 0;
    int num_dead_ends = 0;
    int num_duplicates = 0;
    BOR_ISET(state);
    for (int si = 0; si < cfg.num_samples; ++si){
        stateSamplerSample(&sampler, err);

        borISetEmpty(&state);
        bzero(div->coef, sizeof(double) * pot->var_size);
        for (int var = 0; var < fdr->var.var_size; ++var){
            int id = fdr->var.var[var].val[sampler.state[var]].global_id;
            div->coef[id] = 1.;
            borISetAdd(&state, id);
        }

        if (borHashSetFind(&div->states, &state) >= 0){
            // Ignore duplicates
            ++num_duplicates;
            continue;
        }

        // Compute heuristic estimate
        pddlPotSetObj(pot, div->coef);
        if (solve2(hpot, pot, div->func[num_states]) == 0){
            int h = fdrStateEstimate(div->func[num_states], &fdr->var,
                                     sampler.state);
            if (h != PDDL_COST_DEAD_END){
                // Add state to the set of states and store heuristic estimate
                int state_id = borHashSetAdd(&div->states, &state);
                ASSERT(state_id == num_states);
                div->state_est[state_id] = h;
                ASSERT_RUNTIME(div->state_est[state_id] >= 0);
                ++num_states;

                if ((si + 1) % 100 == 0){
                    BOR_INFO(err, "Pot: Diverse: %d/%d (dead-ends: %d)",
                             num_states, cfg.num_samples, num_dead_ends);
                }

            }else{
                ++num_dead_ends;
            }
        }else{
            // Dead-ends are simply skipped
            ++num_dead_ends;
        }
    }
    BOR_INFO(err, "Pot: Detected dead-ends: %d", num_dead_ends);
    BOR_INFO(err, "Pot: Detected duplicates: %d", num_duplicates);
    ASSERT(num_states == div->states.size);
    div->active_states = div->states.size;
    borISetFree(&state);
    stateSamplerFree(&sampler);
}



static int diverseAvg(diverse_pot_t *div,
                      pddl_hpot_t *hpot,
                      pddl_pot_t *pot,
                      bor_err_t *err)
{
    bzero(div->coef, sizeof(double) * pot->var_size);
    for (int i = 0; i < div->states.size; ++i){
        if (div->state_est[i] < 0)
            continue;
        const bor_iset_t *state = borHashSetGet(&div->states, i);
        int fact_id;
        BOR_ISET_FOR_EACH(state, fact_id)
            div->coef[fact_id] += 1.;
    }
    pddlPotSetObj(pot, div->coef);
    return solve2(hpot, pot, div->avg_func);
}

static const double *diverseSelectFunc(diverse_pot_t *div,
                                       pddl_hpot_t *hpot,
                                       pddl_pot_t *pot,
                                       const pddl_fdr_t *fdr,
                                       bor_err_t *err)
{
    if (diverseAvg(div, hpot, pot, err) != 0)
        return NULL;

    int *fdr_state = BOR_ALLOC_ARR(int, hpot->var_size);
    for (int si = 0; si < div->states.size; ++si){
        if (div->state_est[si] < 0)
            continue;
        const bor_iset_t *state = borHashSetGet(&div->states, si);
        setStateToFDRState(state, fdr_state, fdr);

        int hest = fdrStateEstimate(div->avg_func, &fdr->var, fdr_state);
        if (hest == div->state_est[si]){
            BOR_FREE(fdr_state);
            return div->avg_func;
        }
    }

    int sid = borRandMT(div->rnd, 0, div->active_states);
    for (int si = 0; si < div->states.size; ++si){
        if (div->state_est[si] < 0)
            continue;
        if (sid-- == 0){
            BOR_FREE(fdr_state);
            return div->func[si];
        }
    }
    ASSERT_RUNTIME_M(0, "The number of active states is invalid!");
    return NULL;
}

static void diverseFilterOutStates(diverse_pot_t *div,
                                   const pddl_fdr_t *fdr,
                                   const double *func,
                                   bor_err_t *err)
{
    int *fdr_state = BOR_ALLOC_ARR(int, fdr->var.var_size);
    for (int si = 0; si < div->states.size; ++si){
        if (div->state_est[si] < 0)
            continue;
        const bor_iset_t *state = borHashSetGet(&div->states, si);
        setStateToFDRState(state, fdr_state, fdr);
        int hest = fdrStateEstimate(func, &fdr->var, fdr_state);
        if (hest >= div->state_est[si]){
            div->state_est[si] = -1;
            --div->active_states;
        }
    }
    BOR_FREE(fdr_state);
}

static int diverse(pddl_hpot_t *hpot,
                   pddl_pot_t *pot,
                   const pddl_fdr_t *fdr,
                   const pddl_hpot_config_t *cfg,
                   bor_err_t *err)
{
    BOR_INFO(err, "Pot: Diverse potentials with %d samples", cfg->num_samples);
    ASSERT_RUNTIME(cfg->num_samples > 0);
    diverse_pot_t div;
    diverseInit(&div, pot, fdr, cfg->num_samples);
    diverseGenStates(&div, hpot, pot, fdr, cfg, err);
    while (div.active_states > 0){
        const double *func = diverseSelectFunc(&div, hpot, pot, fdr, err);
        if (func == NULL)
            return -1;
        addFunc2(hpot, func);
        diverseFilterOutStates(&div, fdr, func, err);
    }
    diverseFree(&div, pot, fdr, cfg->num_samples);
    BOR_INFO(err, "Pot: Computed diverse potentials with %d functions",
             hpot->pot_size);
    return 0;
}

int pddlHPotInit(pddl_hpot_t *hpot,
                 const pddl_fdr_t *fdr,
                 const pddl_hpot_config_t *cfg,
                 bor_err_t *err)
{
    int ret = 0;
    if (fdr->has_cond_eff){
        BOR_INFO2(err, "Pot: Conditional effects are not supported");
        return -1;
    }

    init(hpot, fdr->var.global_id_size);

    // Construct MG-Strips and compute h^2 mutexes if necessary
    pddl_mg_strips_t mg_strips;
    pddl_mutex_pairs_t mutex;
    int need_mutex = 0;
    if (cfg->disambiguation
            || cfg->weak_disambiguation
            || cfg->samples_use_mutex
            || cfg->obj == PDDL_HPOT_OBJ_ALL_STATES_MUTEX
            || cfg->obj == PDDL_HPOT_OBJ_ALL_STATES_MUTEX_CONDITIONED
            || cfg->obj == PDDL_HPOT_OBJ_ALL_STATES_MUTEX_CONDITIONED_RAND
            || cfg->obj == PDDL_HPOT_OBJ_ALL_STATES_MUTEX_CONDITIONED_RAND2){
        need_mutex = 1;
        pddlMGStripsInitFDR(&mg_strips, fdr);
        pddlMutexPairsInitStrips(&mutex, &mg_strips.strips);
        pddlMutexPairsAddMGroups(&mutex, &mg_strips.mg);
        pddlH2(&mg_strips.strips, &mutex, NULL, NULL, err);
    }

    // Initialize potential heuristic
    pddl_pot_t pot;
    initPot(hpot, &pot, fdr, &mg_strips, &mutex, cfg, err);

    if (cfg->add_init_constr){
        // Add constraint on the initial state
        if (addInitConstr(hpot, &pot, fdr, cfg, err) != 0){
            pddlPotFree(&pot);
            return -1;
        }
    }

    if (cfg->obj == PDDL_HPOT_OBJ_INIT){
        pddlPotSetObjFDRState(&pot, &fdr->var, fdr->init);
        if ((ret = solve(hpot, &pot)) == 0)
            addFunc(hpot);
        BOR_INFO(err, "Pot: Solved for the initial state: %d", ret);

    }else if (cfg->obj == PDDL_HPOT_OBJ_ALL_STATES){
        pddlPotSetObjFDRAllSyntacticStates(&pot, &fdr->var);
        if ((ret = solve(hpot, &pot)) == 0)
            addFunc(hpot);
        BOR_INFO(err, "Pot: Solved for all states: %d", ret);

    }else if (cfg->obj == PDDL_HPOT_OBJ_MAX_INIT_ALL_STATES){
        pddlPotSetObjFDRState(&pot, &fdr->var, fdr->init);
        if ((ret = solve(hpot, &pot)) == 0)
            addFunc(hpot);
        BOR_INFO(err, "Pot: Solved for the initial state: %d", ret);

        pddlPotSetObjFDRAllSyntacticStates(&pot, &fdr->var);
        if ((ret = solve(hpot, &pot)) == 0)
            addFunc(hpot);
        BOR_INFO(err, "Pot: Solved for all states: %d", ret);

    }else if (cfg->obj == PDDL_HPOT_OBJ_SAMPLES_MAX
                || cfg->obj == PDDL_HPOT_OBJ_SAMPLES_SUM){
        const pddl_mutex_pairs_t *m = NULL;
        if (cfg->samples_use_mutex)
            m = &mutex;
        ret = samples(hpot, &pot, fdr, m, cfg, err);

    }else if (cfg->obj == PDDL_HPOT_OBJ_ALL_STATES_MUTEX){
        setObjAllStatesMutex(&pot, &mg_strips, &mutex,
                             cfg->all_states_mutex_size);
        if (cfg->all_states_mutex_size < 1 || cfg->all_states_mutex_size > 2){
            BOR_FATAL("all-states-mutex with size %d unsupported!",
                      cfg->all_states_mutex_size);
        }
        if ((ret = solve(hpot, &pot)) == 0)
            addFunc(hpot);

    }else if (cfg->obj == PDDL_HPOT_OBJ_ALL_STATES_MUTEX_CONDITIONED){
        BOR_ISET(facts);
        for (int f = 0; f < mg_strips.strips.fact.fact_size; ++f)
            borISetAdd(&facts, f);
        ret = allStatesMutexCond(hpot, &pot, &mg_strips, &mutex,
                                 cfg->all_states_mutex_size, &facts, err);
        borISetFree(&facts);

    }else if (cfg->obj == PDDL_HPOT_OBJ_ALL_STATES_MUTEX_CONDITIONED_RAND){
        bor_rand_mt_t *rnd = borRandMTNew(rand_sampler_seed);
        BOR_ISET(facts);
        int fact_size = mg_strips.strips.fact.fact_size;
        for (int i = 0; i < cfg->num_samples; ++i)
            borISetAdd(&facts, borRandMT(rnd, 0, fact_size));

        ret = allStatesMutexCond(hpot, &pot, &mg_strips, &mutex,
                                 cfg->all_states_mutex_size, &facts, err);
        borISetFree(&facts);

    }else if (cfg->obj == PDDL_HPOT_OBJ_ALL_STATES_MUTEX_CONDITIONED_RAND2){
        ret = allStatesMutexCond2(hpot, &pot, &mg_strips, &mutex,
                                  cfg->all_states_mutex_size,
                                  cfg->num_samples, err);

    }else if (cfg->obj == PDDL_HPOT_OBJ_DIVERSE){
        ret = diverse(hpot, &pot, fdr, cfg, err);

    }else{
        if (need_mutex){
            pddlMutexPairsFree(&mutex);
            pddlMGStripsFree(&mg_strips);
        }
        pddlPotFree(&pot);
        BOR_ERR_RET(err, -1, "Unkown objective function for potential"
                             " heuristic: %d", cfg->obj);
    }

    if (need_mutex){
        pddlMutexPairsFree(&mutex);
        pddlMGStripsFree(&mg_strips);
    }
    pddlPotFree(&pot);

    if (ret != 0)
        BOR_INFO2(err, "Pot: No optimal solution found");

    return ret;
}

double pddlHPotFDRStateEstimateDbl(const pddl_hpot_t *hpot,
                                   const pddl_fdr_vars_t *vars,
                                   const int *state)
{
    if (hpot->pot_size <= 0)
        return -1;

    double est = fdrStateEstimateDbl(hpot->pot[0], vars, state);
    for (int p = 1; p < hpot->pot_size; ++p){
        double e = fdrStateEstimateDbl(hpot->pot[p], vars, state);
        if (e > est)
            est = e;
    }
    return est;
}

int pddlHPotFDRStateEstimate(const pddl_hpot_t *hpot,
                             const pddl_fdr_vars_t *vars,
                             const int *state)
{
    if (hpot->pot_size <= 0)
        return -1;

    int est = fdrStateEstimate(hpot->pot[0], vars, state);
    for (int p = 1; p < hpot->pot_size; ++p){
        int e = fdrStateEstimate(hpot->pot[p], vars, state);
        if (e > est)
            est = e;
    }
    return est;
}

void pddlHPotFree(pddl_hpot_t *hpot)
{
    for (int i = 0; i < hpot->pot_alloc; ++i)
        BOR_FREE(hpot->pot[i]);
    if (hpot->pot != NULL)
        BOR_FREE(hpot->pot);
    if (hpot->func != NULL)
        BOR_FREE(hpot->func);
}
