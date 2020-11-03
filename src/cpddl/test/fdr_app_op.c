#include <sys/time.h>
#include <sys/resource.h>
#include <cu/cu.h>
#include <boruvka/rand.h>
#include "pddl/pddl.h"

static void findOps(const pddl_fdr_app_op_t *app_op,
                    const pddl_fdr_ops_t *ops,
                    const int *state,
                    int depth)
{
    BOR_ISET(app);
    int ret = pddlFDRAppOpFind(app_op, state, &app);
    assertEquals(ret, borISetSize(&app));

    /*
    printf("Init: %d:", ret);
    int op_id;
    BOR_ISET_FOR_EACH(&app, op_id)
        printf(" %d", op_id);
    printf("\n");
    */

    BOR_ISET(app2);
    for (int op_id = 0; op_id < ops->op_size; ++op_id){
        const pddl_fdr_op_t *op = ops->op[op_id];
        int applicable = 1;
        for (int fi = 0; fi < op->pre.fact_size; ++fi){
            const pddl_fdr_fact_t *fact = op->pre.fact + fi;
            if (state[fact->var] != fact->val){
                applicable = 0;
                break;
            }
        }

        if (applicable)
            borISetAdd(&app2, op_id);
    }

    assertTrue(borISetEq(&app, &app2));

    if (depth > 0){
        int *next_state = BOR_ALLOC_ARR(int, app_op->var_size);
        int op_id;
        BOR_ISET_FOR_EACH(&app, op_id){
            memcpy(next_state, state, sizeof(int) * app_op->var_size);

            const pddl_fdr_op_t *op = ops->op[op_id];
            for (int fi = 0; fi < op->eff.fact_size; ++fi)
                next_state[op->eff.fact[fi].var] = op->eff.fact[fi].val;
            findOps(app_op, ops, next_state, depth - 1);
        }
        BOR_FREE(next_state);
    }

    borISetFree(&app2);
    borISetFree(&app);
}

static void findOpsRand(const pddl_fdr_app_op_t *app_op,
                        const pddl_fdr_vars_t *vars,
                        const pddl_fdr_ops_t *ops,
                        int num_samples)
{
    bor_rand_t rnd;
    borRandInit(&rnd);

    int *state = BOR_ALLOC_ARR(int, vars->var_size);
    for (int sample = 0; sample < num_samples; ++sample){
        for (int var = 0; var < vars->var_size; ++var){
            int val = borRand(&rnd, 0, vars->var[var].val_size);
            val = BOR_MIN(val, vars->var[var].val_size - 1);
            state[var] = val;
        }
        findOps(app_op, ops, state, 0);
    }
    BOR_FREE(state);
}

static void run(const char *domain_fn,
                const char *problem_fn,
                unsigned fdr_flag)
{
    pddl_config_t cfg = PDDL_CONFIG_INIT;
    pddl_ground_config_t ground_cfg = PDDL_GROUND_CONFIG_INIT;
    pddl_t pddl;
    bor_err_t err = BOR_ERR_INIT;
    pddl_strips_t strips;
    int ret;

    cfg.force_adl = 1;
    ret = pddlInit(&pddl, domain_fn, problem_fn, &cfg, &err);
    assertEquals(ret, 0);
    if (ret != 0){
        borErrPrint(&err, 1, stderr);
        return;
    }
    pddlCheckSizeTypes(&pddl);
    pddlNormalize(&pddl);

    pddl_lifted_mgroups_infer_limits_t infer_limit
                = PDDL_LIFTED_MGROUPS_INFER_LIMITS_INIT;
    pddl_lifted_mgroups_t lmgs;
    pddlLiftedMGroupsInit(&lmgs);
    pddlLiftedMGroupsInferFAMGroups(&pddl, &infer_limit, &lmgs, &err);
    pddlLiftedMGroupsSetExactlyOne(&pddl, &lmgs, &err);

    ret = pddlStripsGround(&strips, &pddl, &ground_cfg, &err);
    assertEquals(ret, 0);
    pddlStripsRemoveStaticFacts(&strips, &err);


    pddl_mgroups_t mgs;
    pddlMGroupsGround(&mgs, &pddl, &lmgs, &strips);
    pddlMGroupsSetExactlyOne(&mgs, &strips);

    pddl_mutex_pairs_t mutex;
    pddlMutexPairsInit(&mutex, strips.fact.fact_size);
    pddlMutexPairsAddMGroups(&mutex, &mgs);


    pddl_fdr_t fdr;
    pddlFDRInitFromStrips(&fdr, &strips, &mgs, &mutex, fdr_flag, &err);

    pddl_fdr_app_op_t app_op;
    pddlFDRAppOpInit(&app_op, &fdr.var, &fdr.op, &fdr.goal);
    findOps(&app_op, &fdr.op, fdr.init, 2);
    findOpsRand(&app_op, &fdr.var, &fdr.op, 5000);

    pddlFDRAppOpFree(&app_op);
    pddlFDRFree(&fdr);

    pddlMutexPairsFree(&mutex);
    pddlMGroupsFree(&mgs);
    pddlStripsFree(&strips);
    pddlLiftedMGroupsFree(&lmgs);
    pddlFree(&pddl);
}

static void setMemLimit(void)
{
    struct rlimit mem_limit;

    mem_limit.rlim_cur = mem_limit.rlim_max = 512 * 1024 * 1024;
    mem_limit.rlim_cur = mem_limit.rlim_max = 2024 * 1024 * 1024;
    setrlimit(RLIMIT_AS, &mem_limit);
}

TEST(testFDRAppOpSetUp)
{
    setMemLimit();
}

#define P(N, P) \
TEST(testFDRAppOp_##N) \
{ \
    unsigned fdr_flag = PDDL_FDR_VARS_LARGEST_FIRST; \
    pddl_files_t files; \
    bor_err_t err = BOR_ERR_INIT; \
    if (pddlFiles(&files, "pddl-data/", P, &err) != 0){ \
        borErrPrint(&err, 1, stderr); \
        return; \
    } \
    run(files.domain_pddl, files.problem_pddl, fdr_flag); \
}
#include "fdr_app_op_prob.h"
