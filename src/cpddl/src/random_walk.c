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

#include "pddl/random_walk.h"
#include "assert.h"

static void init(pddl_random_walk_t *rndw,
                 const pddl_fdr_t *fdr,
                 const pddl_fdr_app_op_t *app_op)
{
    bzero(rndw, sizeof(*rndw));
    rndw->fdr = fdr;
    if (app_op == NULL){
        pddl_fdr_app_op_t *app = BOR_ALLOC(pddl_fdr_app_op_t);
        pddlFDRAppOpInit(app, &fdr->var, &fdr->op, &fdr->goal);
        rndw->app = app;
        rndw->owns_app_op = 1;
    }else{
        rndw->app = app_op;
        rndw->owns_app_op = 0;
    }
}

void pddlRandomWalkInit(pddl_random_walk_t *rndw,
                        const pddl_fdr_t *fdr,
                        const pddl_fdr_app_op_t *app_op)
{
    init(rndw, fdr, app_op);
    rndw->rnd = borRandMTNewAuto();
}

void pddlRandomWalkInitSeed(pddl_random_walk_t *rndw,
                            const pddl_fdr_t *fdr,
                            const pddl_fdr_app_op_t *app_op,
                            uint32_t seed)
{
    init(rndw, fdr, app_op);
    rndw->rnd = borRandMTNew(seed);
}

void pddlRandomWalkFree(pddl_random_walk_t *rndw)
{
    if (rndw->owns_app_op && rndw->app != NULL){
        pddl_fdr_app_op_t *app = (pddl_fdr_app_op_t *)rndw->app;
        pddlFDRAppOpFree(app);
        BOR_FREE(app);
    }

    if (rndw->rnd != NULL)
        borRandMTDel(rndw->rnd);
}

int pddlRandomWalkSampleState(pddl_random_walk_t *rndw,
                              const int *start_state,
                              int max_steps,
                              int *resulting_state)
{
    static const double p = 0.5;

    // Length of the walk
    int num_steps = 0;
    for (int i = 0; i < max_steps; ++i){
        if (borRandMT01(rndw->rnd) < p)
            ++num_steps;
    }

    // Perform at most num_steps of a random walk
    int num_performed_steps = 0;
    BOR_ISET(app_ops);
    memcpy(resulting_state, start_state, sizeof(int) * rndw->fdr->var.var_size);
    for (int step = 0; step < num_steps; ++step){
        borISetEmpty(&app_ops);
        pddlFDRAppOpFind(rndw->app, resulting_state, &app_ops);
        if (borISetSize(&app_ops) == 0){
            // No applicable operators -- terminate
            break;
        }else{
            int which_op = borRandMT(rndw->rnd, 0, borISetSize(&app_ops));
            ASSERT(which_op >= 0 && which_op < borISetSize(&app_ops));
            int op_id = borISetGet(&app_ops, which_op);
            const pddl_fdr_op_t *op = rndw->fdr->op.op[op_id];
            pddlFDROpApplyOnState(op, resulting_state);
            ++num_performed_steps;
        }
    }
    borISetFree(&app_ops);

    return num_performed_steps;
}
