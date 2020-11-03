#include <stdio.h>
#include <pddl/pddl.h>
#include <opts.h>

struct options {
    int help;
    int fd_fam_groups;
    int fam_groups;
    char *fdr_out;
    char *eval_plan;
    pddl_files_t files;
} opt;

static void usage(const char *name)
{
    fprintf(stderr, "pddl-pot for computing potential heuristic.\n");
    fprintf(stderr, "Usage: %s [OPTIONS] [domain.pddl] problem.pddl\n", name);
    fprintf(stderr, "  OPTIONS:\n");
    optsPrint(stderr, "    ");
    fprintf(stderr, "\n");
}

static int readOpts(int *argc,
                    char *argv[],
                    pddl_hpot_config_t *pot_cfg,
                    bor_err_t *err)
{
    int disamb = 0;
    int weak_disamb = 0;
    int no_disamb = 0;
    int obj_init = 0;
    int obj_all_states = 0;
    int obj_max_init_all_states = 0;
    double add_init_constr = -1.;
    int obj_samples_sum = -1;
    int obj_samples_max = -1;
    int obj_samples_mutex = 0;
    int obj_all_states_mutex = 0;
    int obj_all_states_mutex_cond = 0;
    int obj_all_states_mutex_cond_rand[2];
    int obj_all_states_mutex_cond_rand2[2];
    int obj_diverse = 0;

    obj_all_states_mutex_cond_rand[0] = 0;
    obj_all_states_mutex_cond_rand[1] = 0;
    obj_all_states_mutex_cond_rand2[0] = 0;
    obj_all_states_mutex_cond_rand2[1] = 0;

    optsAddDesc("help", 'h', OPTS_NONE, &opt.help, NULL,
                "Print this help.");
    optsAddDesc("output", 'o', OPTS_STR, &opt.fdr_out, NULL,
                "Output filename (default: stdout)");
    optsAddDesc("eval-plan", 0x0, OPTS_STR, &opt.eval_plan, NULL,
                "Evaluate given plan");
    optsAddDesc("fam-groups", 'f', OPTS_NONE, &opt.fam_groups, NULL,
                "Use LP to infer all maximal fam-groups.");
    optsAddDesc("fd", 0x0, OPTS_NONE, &opt.fd_fam_groups, NULL,
                "Use Fast-Downward fam-groups.");

    optsAddDesc("disamb", 'd', OPTS_NONE, &disamb, NULL,
                "Enable disambiguation.");
    optsAddDesc("weak-disamb", 'w', OPTS_NONE, &weak_disamb, NULL,
                "Enable weak disambiguation.");
    optsAddDesc("no-disamb", 'n', OPTS_NONE, &no_disamb, NULL,
                "Disable disambiguation.");

    optsAddDesc("init-state", 'I', OPTS_NONE, &obj_init, NULL,
                "Optimize for the initial state");
    optsAddDesc("all-states", 'A', OPTS_NONE, &obj_all_states, NULL,
                "Optimize for all syntactic states");
    optsAddDesc("max-init-all-states", 'X', OPTS_NONE,
                &obj_max_init_all_states, NULL,
                "Maximum of -I and -A");
    optsAddDesc("samples-sum", 'S', OPTS_INT, &obj_samples_sum, NULL,
                "Optimize for the sum of samples.");
    optsAddDesc("samples-max", 'T', OPTS_INT, &obj_samples_max, NULL,
                "Optimize for each samples.");
    optsAddDesc("samples-use-mutex", 0x0, OPTS_NONE, &obj_samples_mutex, NULL,
                "Use mutexes to filter out unreachable states.");
    optsAddDesc("all-states-mutex", 'M', OPTS_INT, &obj_all_states_mutex, NULL,
                "Optimize for all syntactic states using mutexes");
    optsAddDesc("all-states-mutex-cond", 'C', OPTS_INT,
                &obj_all_states_mutex_cond, NULL,
                "Optimize for all syntactic states conditioned on facts.");
    optsAddDesc("all-states-mutex-cond-rand", 'R', OPTS_INT_ARR(2),
                obj_all_states_mutex_cond_rand, NULL,
                "Optimize for all syntactic states conditioned on facts.");
    optsAddDesc("all-states-mutex-cond-rand2", 'P', OPTS_INT_ARR(2),
                obj_all_states_mutex_cond_rand2, NULL,
                "Optimize for all syntactic states conditioned on"
                " pair of facts.");
    optsAddDesc("diverse", 'D', OPTS_INT, &obj_diverse, NULL,
                "Diverse potentials");

    optsAddDesc("add-init-constr", 'L', OPTS_DOUBLE, &add_init_constr, NULL,
                "Add lower bound constraint on the initial state.");


    if (opts(argc, argv) != 0 || opt.help || (*argc != 2 && *argc != 3)){
        if (*argc <= 1)
            fprintf(stderr, "Error: Missing input file(s)\n\n");

        if (*argc > 3){
            for (int i = 0; i < *argc; ++i){
                if (argv[i][0] == '-'){
                    fprintf(stderr, "Error: Unrecognized option '%s'\n",
                            argv[i]);
                }
            }
        }

        usage(argv[0]);
        return -1;
    }

    if (disamb + weak_disamb + no_disamb != 1){
        fprintf(stderr, "Error: One of -d/-w/-n must be specified!\n\n");
        usage(argv[0]);
        return -1;
    }else if (disamb){
        pot_cfg->disambiguation = 1;
        pot_cfg->weak_disambiguation = 0;
    }else if (weak_disamb){
        pot_cfg->disambiguation = 0;
        pot_cfg->weak_disambiguation = 1;
    }else{
        pot_cfg->disambiguation = 0;
        pot_cfg->weak_disambiguation = 0;
    }

    pot_cfg->samples_use_mutex = 0;
    if (obj_init
            + obj_all_states
            + obj_max_init_all_states
            + (obj_samples_sum > 0 ? 1 : 0)
            + (obj_samples_max > 0 ? 1 : 0)
            + (obj_all_states_mutex > 0 ? 1 : 0)
            + (obj_all_states_mutex_cond > 0 ? 1 : 0)
            + (obj_all_states_mutex_cond_rand[1] > 0 ? 1 : 0)
            + (obj_all_states_mutex_cond_rand2[1] > 0 ? 1 : 0)
            + (obj_diverse > 0 ? 1 : 0) != 1){
        fprintf(stderr, "Error: One of -I/-A/-X/-S/-T/-M/-D/-C/-R/-P"
                        " must be specified!\n\n");
        usage(argv[0]);
        return -1;

    }else if (obj_init){
        pot_cfg->obj = PDDL_HPOT_OBJ_INIT;

    }else if (obj_all_states){
        pot_cfg->obj = PDDL_HPOT_OBJ_ALL_STATES;

    }else if (obj_max_init_all_states){
        pot_cfg->obj = PDDL_HPOT_OBJ_MAX_INIT_ALL_STATES;

    }else if (obj_samples_sum > 0){
        pot_cfg->obj = PDDL_HPOT_OBJ_SAMPLES_SUM;
        pot_cfg->num_samples = obj_samples_sum;
        if (obj_samples_mutex){
            pot_cfg->samples_use_mutex = 1;
        }else{
            pot_cfg->samples_random_walk = 1;
        }

    }else if (obj_samples_max > 0){
        pot_cfg->obj = PDDL_HPOT_OBJ_SAMPLES_MAX;
        pot_cfg->num_samples = obj_samples_max;
        if (obj_samples_mutex){
            pot_cfg->samples_use_mutex = 1;
        }else{
            pot_cfg->samples_random_walk = 1;
        }

    }else if (obj_all_states_mutex > 0){
        pot_cfg->obj = PDDL_HPOT_OBJ_ALL_STATES_MUTEX;
        pot_cfg->all_states_mutex_size = obj_all_states_mutex;

    }else if (obj_all_states_mutex_cond > 0){
        pot_cfg->obj = PDDL_HPOT_OBJ_ALL_STATES_MUTEX_CONDITIONED;
        pot_cfg->all_states_mutex_size = obj_all_states_mutex_cond;

    }else if (obj_all_states_mutex_cond_rand[0] > 0){
        pot_cfg->obj = PDDL_HPOT_OBJ_ALL_STATES_MUTEX_CONDITIONED_RAND;
        pot_cfg->all_states_mutex_size = obj_all_states_mutex_cond_rand[0];
        pot_cfg->num_samples = obj_all_states_mutex_cond_rand[1];

    }else if (obj_all_states_mutex_cond_rand2[0] > 0){
        pot_cfg->obj = PDDL_HPOT_OBJ_ALL_STATES_MUTEX_CONDITIONED_RAND2;
        pot_cfg->all_states_mutex_size = obj_all_states_mutex_cond_rand2[0];
        pot_cfg->num_samples = obj_all_states_mutex_cond_rand2[1];

    }else if (obj_diverse > 0){
        pot_cfg->obj = PDDL_HPOT_OBJ_DIVERSE;
        pot_cfg->num_samples = obj_diverse;
    }


    if (add_init_constr < 0.){
        pot_cfg->add_init_constr = 0;
    }else{
        pot_cfg->add_init_constr = 1;
        pot_cfg->init_constr_coef = add_init_constr;
    }


    if (*argc == 2){
        BOR_INFO(err, "Input file: '%s'", argv[1]);
        if (pddlFiles1(&opt.files, argv[1], err) != 0)
            BOR_TRACE_RET(err, -1);
    }else{ // *argc == 3
        BOR_INFO(err, "Input files: '%s' and '%s'", argv[1], argv[2]);
        if (pddlFiles(&opt.files, argv[1], argv[2], err) != 0)
            BOR_TRACE_RET(err, -1);
    }
    BOR_INFO(err, "PDDL files: '%s' '%s'\n",
             opt.files.domain_pddl, opt.files.problem_pddl);

    return 0;
}

static void printPotentials(const pddl_fdr_t *fdr,
                            const pddl_hpot_t *hpot,
                            FILE *fout)
{
    fprintf(fout, "%d\n", hpot->pot_size);
    for (int pi = 0; pi < hpot->pot_size; ++pi){
        const double *w = hpot->pot[pi];
        fprintf(fout, "begin_potentials\n");
        for (int fi = 0; fi < fdr->var.global_id_size; ++fi){
            const pddl_fdr_val_t *fval = fdr->var.global_id_to_val[fi];
            fprintf(fout, "%d %d %.20f\n",
                    fval->var_id, fval->val_id, w[fi]);
        }
        fprintf(fout, "end_potentials\n");
    }
}

static void evalPlan(const pddl_fdr_t *fdr,
                     const pddl_hpot_t *hpot)
{
    FILE *fin = fopen(opt.eval_plan, "r");
    if (fin == NULL){
        fprintf(stderr, "Error: Could not open '%s'\n", opt.eval_plan);
        exit(-1);
    }

    char *line = NULL;
    size_t len = 0;
    ssize_t nread;

    int *state = BOR_CALLOC_ARR(int, fdr->var.var_size);
    memcpy(state, fdr->init, sizeof(int) * fdr->var.var_size);
    int est = pddlHPotFDRStateEstimate(hpot, &fdr->var, state);
    double est2 = 0.;
    for (int v = 0; v < fdr->var.var_size; ++v){
        const pddl_fdr_val_t *val = fdr->var.var[v].val + state[v];
        est2 += hpot->pot[0][val->global_id];
    }
    printf("S: %d | %f", est, est2);
    /*
    for (int v = 0; v < fdr->var.var_size; ++v)
        printf(" %d", state[v]);
    */
    printf("\n");
    double last_est = est2;

    while ((nread = getline(&line, &len, fin)) != -1) {
        if (nread == 0)
            continue;
        if (line[0] != '(' || line[nread - 2] != ')')
            continue;
        line[nread - 2] = 0x0;
        const char *op_name = line + 1;
        int op_id = -1;
        for (int oi = 0; oi < fdr->op.op_size; ++oi){
            if (strcmp(fdr->op.op[oi]->name, op_name) == 0){
                if (op_id >= 0){
                    fprintf(stderr, "Two operators with the same name!\n");
                    exit(-1);
                }
                op_id = oi;
            }
        }

        const pddl_fdr_op_t *op = fdr->op.op[op_id];
        printf("O: (%s) %d, cost: %d | pre:", op->name, op_id, op->cost);
        for (int i = 0; i < op->pre.fact_size; ++i){
            const pddl_fdr_fact_t *f = op->pre.fact + i;
            if (state[f->var] != f->val){
                printf("Operator not applicable!\n");
                exit(-1);
            }
            printf(" %d:%d", f->var, f->val);
        }
        printf(" | eff:");

        for (int i = 0; i < op->eff.fact_size; ++i){
            const pddl_fdr_fact_t *f = op->eff.fact + i;
            state[f->var] = f->val;
            printf(" %d:%d", f->var, f->val);
        }
        printf("\n");

        int est = pddlHPotFDRStateEstimate(hpot, &fdr->var, state);
        double est2 = 0.;
        for (int v = 0; v < fdr->var.var_size; ++v){
            const pddl_fdr_val_t *val = fdr->var.var[v].val + state[v];
            est2 += hpot->pot[0][val->global_id];
        }
        if (last_est > est2 + op->cost)
            printf("********\n");
        printf("S: %d | %f", est, est2);
        for (int v = 0; v < fdr->var.var_size; ++v)
            printf(" %d:%d", v, state[v]);
        printf("\n");
        last_est = est2;
    }

    BOR_FREE(state);
}

int main(int argc, char *argv[])
{
    pddl_hpot_config_t hpot_cfg = PDDL_HPOT_CONFIG_INIT;

    bor_err_t err = BOR_ERR_INIT;
    borErrWarnEnable(&err, stderr);
    borErrInfoEnable(&err, stderr);

    if (readOpts(&argc, argv, &hpot_cfg, &err) != 0){
        borErrPrint(&err, 1, stderr);
        return -1;
    }

    // Parse PDDL
    pddl_config_t pddl_cfg = PDDL_CONFIG_INIT;
    pddl_cfg.force_adl = 1;
    pddl_t pddl;
    if (pddlInit(&pddl, opt.files.domain_pddl, opt.files.problem_pddl,
                 &pddl_cfg, &err) != 0){
        fprintf(stderr, "Error: ");
        borErrPrint(&err, 1, stderr);
        return -1;
    }
    pddlNormalize(&pddl);
    pddlCheckSizeTypes(&pddl);

    // Lifted mgroups
    pddl_lifted_mgroups_infer_limits_t lifted_mgroups_limits
                = PDDL_LIFTED_MGROUPS_INFER_LIMITS_INIT;
    pddl_lifted_mgroups_t lifted_mgroups;
    pddlLiftedMGroupsInit(&lifted_mgroups);
    if (opt.fd_fam_groups){
        pddlLiftedMGroupsInferMonotonicity(&pddl, &lifted_mgroups_limits,
                                           NULL, &lifted_mgroups, &err);
    }else{
        pddlLiftedMGroupsInferFAMGroups(&pddl, &lifted_mgroups_limits,
                                        &lifted_mgroups, &err);
    }
    pddlLiftedMGroupsSetExactlyOne(&pddl, &lifted_mgroups, &err);
    pddlLiftedMGroupsSetStatic(&pddl, &lifted_mgroups, &err);

    // Ground to STRIPS
    pddl_ground_config_t ground_cfg = PDDL_GROUND_CONFIG_INIT;
    ground_cfg.lifted_mgroups = &lifted_mgroups;
    ground_cfg.prune_op_pre_mutex = 1;
    ground_cfg.prune_op_dead_end = 1;
    pddl_strips_t strips;
    if (pddlStripsGround(&strips, &pddl, &ground_cfg, &err) != 0){
        BOR_INFO2(&err, "Grounding failed.");
        fprintf(stderr, "Error: ");
        borErrPrint(&err, 1, stderr);
        return -1;
    }
    //if (strips.has_cond_eff)
    //    pddlStripsCompileAwayCondEff(&strips);
    if (strips.has_cond_eff){
        BOR_INFO2(&err, "Has conditional effects -- terminating...");
        return -1;
    }

    // Ground mutex groups
    pddl_mgroups_t mgroups;
    pddlMGroupsGround(&mgroups, &pddl, &lifted_mgroups, &strips);
    pddlMGroupsSetExactlyOne(&mgroups, &strips);
    pddlMGroupsSetGoal(&mgroups, &strips);

    // Prune strips
    if (!strips.has_cond_eff){
        BOR_ISET(rm_fact);
        BOR_ISET(rm_op);
        if (pddlIrrelevanceAnalysis(&strips, &rm_fact, &rm_op, NULL, &err) != 0){
            BOR_INFO2(&err, "Irrelevance analysis failed.");
            fprintf(stderr, "Error: ");
            borErrPrint(&err, 1, stderr);
            return -1;
        }
        if (borISetSize(&rm_fact) > 0 || borISetSize(&rm_op) > 0){
            pddlStripsReduce(&strips, &rm_fact, &rm_op);
            if (borISetSize(&rm_fact) > 0)
                pddlMGroupsReduce(&mgroups, &rm_fact);
        }
        borISetFree(&rm_fact);
        borISetFree(&rm_op);
    }

    // Find fam-groups
    if (opt.fam_groups){
        pddl_famgroup_config_t fam_cfg = PDDL_FAMGROUP_CONFIG_INIT;
        if (pddlFAMGroupsInfer(&mgroups, &strips, &fam_cfg, &err) != 0){
            fprintf(stderr, "Error: ");
            borErrPrint(&err, 1, stderr);
            return -1;
        }
    }

    // Construct FDR
    pddl_fdr_t fdr;
    unsigned fdr_var_flag = PDDL_FDR_VARS_LARGEST_FIRST;
    pddl_mutex_pairs_t mutex;
    pddlMutexPairsInitStrips(&mutex, &strips);
    pddlMutexPairsAddMGroups(&mutex, &mgroups);
    pddlFDRInitFromStrips(&fdr, &strips, &mgroups, &mutex, fdr_var_flag, &err);
    if (pddlPruneFDR(&fdr, &err) != 0){
        BOR_INFO2(&err, "Pruning failed.");
        fprintf(stderr, "Error: ");
        borErrPrint(&err, 1, stderr);
        return -1;
    }

    BOR_INFO(&err, "Number of operators: %d", fdr.op.op_size);
    BOR_INFO(&err, "Number of variables: %d", fdr.var.var_size);
    BOR_INFO(&err, "Number of facts: %d", fdr.var.global_id_size);

    pddl_hpot_t hpot;

    if (pddlHPotInit(&hpot, &fdr, &hpot_cfg, &err) != 0){
        BOR_INFO2(&err, "Cannot find potential heuristic");
        return -1;
    }
    int est = pddlHPotFDRStateEstimate(&hpot, &fdr.var, fdr.init);
    BOR_INFO(&err, "Init state estimate: %d", est);

    // Print out FDR in fast-downward format and potentials
    FILE *fout = stdout;
    if (opt.fdr_out != NULL){
        fout = fopen(opt.fdr_out, "w");
        if (fout == NULL){
            fprintf(stderr, "Error: Could not open '%s'", opt.fdr_out);
            return -1;
        }
    }
    pddlFDRPrintFD(&fdr, NULL, fout);
    printPotentials(&fdr, &hpot, fout);
    if (fout != stdout)
        fclose(fout);

    if (opt.eval_plan != NULL)
        evalPlan(&fdr, &hpot);

    pddlHPotFree(&hpot);

    optsClear();
    pddlFDRFree(&fdr);
    pddlMutexPairsFree(&mutex);
    pddlMGroupsFree(&mgroups);
    pddlStripsFree(&strips);
    pddlLiftedMGroupsFree(&lifted_mgroups);
    pddlFree(&pddl);
    return 0;
}



