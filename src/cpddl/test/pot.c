#include <sys/time.h>
#include <sys/resource.h>
#include <boruvka/alloc.h>
#include <cu/cu.h>
#include "pddl/pddl.h"

#define ROUND_EPS 0.001

static int optimalCost(const char *problem_fn)
{
    int cost = -1;
    int fnlen = strlen(problem_fn);
    char *fn = BOR_ALLOC_ARR(char, fnlen + 2);
    strcpy(fn, problem_fn);
    strcpy(fn + fnlen - 4, "plan");

    FILE *fin = fopen(fn, "r");
    if (fin == NULL){
        BOR_FREE(fn);
        return -1;
    }

    size_t len = 0;
    char *line = NULL;
    ssize_t nread;

    while ((nread = getline(&line, &len, fin)) != -1){
        char *f = strstr(line, "Optimal cost:");
        if (f != NULL){
            char *p = f + strlen("Optimal cost:");
            for (; *p == ' '; ++p);
            char *end = p;
            for (; *end >= '0' && *end <= '9'; ++end);
            *end = 0;
            cost = atoi(p);
            break;
        }
    }
    fclose(fin);

    if (line != NULL)
        free(line);

    BOR_FREE(fn);
    return cost;
}

static int roundOff(double z)
{
    return ceil(z - ROUND_EPS);
}

static int potFDRState(const pddl_pot_t *pot,
                       const pddl_fdr_t *fdr,
                       const int *state,
                       const double *w)
{
    double p = 0;
    for (int var = 0; var < fdr->var.var_size; ++var)
        p += w[fdr->var.var[var].val[state[var]].global_id];
    if (p < 0.)
        return 0;
    return roundOff(p);
}

static int potStripsState(const pddl_pot_t *pot,
                          const bor_iset_t *state,
                          const double *w)
{
    double p = 0.;
    int fact_id;
    BOR_ISET_FOR_EACH(state, fact_id)
        p += w[fact_id];
    if (p <= 0.)
        return 0;
    return roundOff(p);
}

static void potFDR(const char *domain_fn, const char *problem_fn,
                   const char *outfn,
                   int compile_away_cond_eff_strips)
{
    pddl_config_t cfg = PDDL_CONFIG_INIT;
    pddl_ground_config_t ground_cfg = PDDL_GROUND_CONFIG_INIT;
    pddl_t pddl;
    bor_err_t err = BOR_ERR_INIT;
    pddl_strips_t strips;
    int ret;

    //borErrInfoEnable(&err, stderr);
    cfg.force_adl = 1;
    ret = pddlInit(&pddl, domain_fn, problem_fn, &cfg, &err);
    assertEquals(ret, 0);
    if (ret != 0){
        borErrPrint(&err, 1, stderr);
        return;
    }
    pddlCheckSizeTypes(&pddl);

    pddlNormalize(&pddl);

    pddl_lifted_mgroups_infer_limits_t lifted_mgroups_limits
            = PDDL_LIFTED_MGROUPS_INFER_LIMITS_INIT;
    pddl_lifted_mgroups_t lifted_mgroups;
    pddlLiftedMGroupsInit(&lifted_mgroups);
    pddlLiftedMGroupsInferFAMGroups(&pddl, &lifted_mgroups_limits,
                                    &lifted_mgroups, &err);
    pddlLiftedMGroupsSetExactlyOne(&pddl, &lifted_mgroups, &err);
    pddlLiftedMGroupsSetStatic(&pddl, &lifted_mgroups, &err);

    ground_cfg.lifted_mgroups = &lifted_mgroups;
    ground_cfg.prune_op_pre_mutex = 1;
    ground_cfg.prune_op_dead_end = 1;
    ret = pddlStripsGround(&strips, &pddl, &ground_cfg, &err);
    assertEquals(ret, 0);
    if (compile_away_cond_eff_strips)
        pddlStripsCompileAwayCondEff(&strips);

    BOR_ISET(irr_facts);
    BOR_ISET(irr_ops);
    pddlIrrelevanceAnalysis(&strips, &irr_facts, &irr_ops, NULL, &err);
    pddlStripsReduce(&strips, &irr_facts, &irr_ops);
    borISetFree(&irr_facts);
    borISetFree(&irr_ops);

    pddl_mgroups_t mgroups;
    pddl_mutex_pairs_t mutex;

    pddlMGroupsGround(&mgroups, &pddl, &lifted_mgroups, &strips);
    pddlMGroupsSetExactlyOne(&mgroups, &strips);
    pddlMGroupsSetGoal(&mgroups, &strips);

    pddlMutexPairsInitStrips(&mutex, &strips);
    pddlMutexPairsAddMGroups(&mutex, &mgroups);

    BOR_ISET(unreachable_fact);
    BOR_ISET(unreachable_op);
    pddlH2FwBw(&strips, &mgroups, &mutex,
               &unreachable_fact, &unreachable_op, &err);
    pddlStripsReduce(&strips, &unreachable_fact, &unreachable_op);
    pddlMGroupsReduce(&mgroups, &unreachable_fact);
    pddlMutexPairsReduce(&mutex, &unreachable_fact);
    borISetFree(&unreachable_fact);
    borISetFree(&unreachable_op);

    pddl_fdr_t fdr;
    unsigned fdr_var_flag = PDDL_FDR_VARS_LARGEST_FIRST;
    pddlFDRInitFromStrips(&fdr, &strips, &mgroups, &mutex, fdr_var_flag, &err);

    int optimal_cost = optimalCost(problem_fn);
    fprintf(stdout, "Optimal cost: %d\n", optimal_cost);

    pddl_pot_t pot;
    pddlPotInitFDR(&pot, &fdr);
    pddlPotSetObjFDRState(&pot, &fdr.var, fdr.init);

    double *w = BOR_ALLOC_ARR(double, pot.var_size);
    ret = pddlPotSolve(&pot, w, pot.var_size, 0);
    assertEquals(ret, 0);

    int val_init_state = potFDRState(&pot, &fdr, fdr.init, w);
    assertTrue(val_init_state <= optimal_cost);
    fprintf(stdout, "Init state: %d\n", val_init_state);


    pddlPotSetObjFDRAllSyntacticStates(&pot, &fdr.var);
    ret = pddlPotSolve(&pot, w, pot.var_size, 0);
    assertEquals(ret, 0);

    int val_all_synt_states = potFDRState(&pot, &fdr, fdr.init, w);
    assertTrue(val_all_synt_states <= optimal_cost);
    assertTrue(val_all_synt_states <= val_init_state);
    fprintf(stdout, "All syntactic states: %d\n", val_all_synt_states);

    pddlPotFree(&pot);
    BOR_FREE(w);

    pddl_mg_strips_t mg_strips2;
    pddlMGStripsInitFDR(&mg_strips2, &fdr);
    pddlMutexPairsFree(&mutex);
    pddlMutexPairsInitStrips(&mutex, &mg_strips2.strips);
    //pddlH2(&mg_strips2.strips, &mutex, NULL, NULL, &err);
    pddlMutexPairsAddMGroups(&mutex, &mg_strips2.mg);
    pddlPotInitMGStrips(&pot, &mg_strips2, &mutex);
    pddlPotSetObjStripsState(&pot, &mg_strips2.strips.init);
    w = BOR_ALLOC_ARR(double, pot.var_size);
    pddlPotSolve(&pot, w, pot.var_size, 0);
    int val_init_fdr_strips = potStripsState(&pot, &mg_strips2.strips.init, w);
    assertTrue(val_init_fdr_strips >= val_init_state);
    if (val_init_fdr_strips != val_init_state)
        fprintf(stdout, "Init state 2: %d\n", val_init_fdr_strips);
    pddlMGStripsFree(&mg_strips2);
    BOR_FREE(w);

    pddl_mg_strips_t mg_strips;
    pddlMGStripsInit(&mg_strips, &strips, &mgroups);

    pddlMutexPairsFree(&mutex);
    pddlMutexPairsInitStrips(&mutex, &mg_strips.strips);
    pddlMutexPairsAddMGroups(&mutex, &mg_strips.mg);

    //borErrInfoEnable(&err, stdout);
    ret = pddlH2(&mg_strips.strips, &mutex, NULL, NULL, &err);

    ret = pddlPotInitMGStrips(&pot, &mg_strips, &mutex);
    assertTrue(ret == 0);
    w = BOR_ALLOC_ARR(double, pot.var_size);

    pddlPotSetObjStripsState(&pot, &mg_strips.strips.init);
    ret = pddlPotSolve(&pot, w, pot.var_size, 0);
    assertEquals(ret, 0);

    int val_mgs_init_state = potStripsState(&pot, &mg_strips.strips.init, w);
    assertTrue(val_mgs_init_state >= val_init_state);
    fprintf(stdout, "Init state strips: %d\n", val_mgs_init_state);

    pddlMGStripsFree(&mg_strips);

    BOR_FREE(w);
    pddlPotFree(&pot);

    pddlFDRFree(&fdr);
    pddlMutexPairsFree(&mutex);
    pddlMGroupsFree(&mgroups);
    pddlLiftedMGroupsFree(&lifted_mgroups);
    pddlFree(&pddl);
    pddlStripsFree(&strips);
}

static void setMemLimit(void)
{
    struct rlimit mem_limit;

    mem_limit.rlim_cur = mem_limit.rlim_max = 512 * 1024 * 1024;
    mem_limit.rlim_cur = mem_limit.rlim_max = 2024 * 1024 * 1024;
    setrlimit(RLIMIT_AS, &mem_limit);
}

TEST(testPotSetUp)
{
    setMemLimit();
}

#define P(N, P) \
TEST(testPotFDR_##N) \
{ \
    pddl_files_t files; \
    bor_err_t err = BOR_ERR_INIT; \
    if (pddlFiles(&files, "pddl-data/", P, &err) != 0){ \
        borErrPrint(&err, 1, stderr); \
        return; \
    } \
    potFDR(files.domain_pddl, files.problem_pddl, #N, 0); \
}

#define NCE(N, P) \
TEST(testPotFDR_##N##_noce_strips) \
{ \
    pddl_files_t files; \
    bor_err_t err = BOR_ERR_INIT; \
    if (pddlFiles(&files, "pddl-data/", P, &err) != 0){ \
        borErrPrint(&err, 1, stderr); \
        return; \
    } \
    potFDR(files.domain_pddl, files.problem_pddl, #N, 1); \
}
#include "pot_prob.h"
