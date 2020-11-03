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

#ifndef __PDDL_HPOT_H__
#define __PDDL_HPOT_H__

#include <pddl/pot.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

struct pddl_hpot {
    double **pot; /*!< Potential functions */
    int pot_size; /*!< Number of potential functions */
    int pot_alloc;
    int var_size; /*!< Number of LP variables in the problem */
    double *func; /*!< Auxiliary potential function */
};
typedef struct pddl_hpot pddl_hpot_t;

#define PDDL_HPOT_OBJ_INIT 0x1
#define PDDL_HPOT_OBJ_ALL_STATES 0x2
#define PDDL_HPOT_OBJ_SAMPLES_MAX 0x3
#define PDDL_HPOT_OBJ_SAMPLES_SUM 0x4
#define PDDL_HPOT_OBJ_ALL_STATES_MUTEX 0x5
#define PDDL_HPOT_OBJ_DIVERSE 0x6
#define PDDL_HPOT_OBJ_ALL_STATES_MUTEX_CONDITIONED 0x7
#define PDDL_HPOT_OBJ_ALL_STATES_MUTEX_CONDITIONED_RAND 0x8
#define PDDL_HPOT_OBJ_ALL_STATES_MUTEX_CONDITIONED_RAND2 0x9
#define PDDL_HPOT_OBJ_MAX_INIT_ALL_STATES 0xa

struct pddl_hpot_config {
    int disambiguation; /*!< If true, disambiguation is used */
    int weak_disambiguation; /*!< If true, weak disambiguation is used */
    int obj; /*!< One of PDDL_HPOT_OBJ_*: specifies optimization method */
    int add_init_constr; /*!< Add >= constraint on the initial state */
    double init_constr_coef; /*!< Coeficient used for the initial state
                                  constraint */
    int num_samples; /*!< Number of samples used for sampling based methods */
    int samples_use_mutex; /*!< If true, mutexes are used to filter out
                                unreachable sample states */
    int samples_random_walk; /*!< Uses random walk for sampling */
    int all_states_mutex_size; /*!< Size of sets of facts for
                                    *_ALL_STATES_MUTEX method */
};
typedef struct pddl_hpot_config pddl_hpot_config_t;

#define PDDL_HPOT_CONFIG_INIT { \
        1, /* .disambiguation */ \
        0, /* .weak_disambiguation */ \
        PDDL_HPOT_OBJ_ALL_STATES, /* .obj */ \
        1, /* .add_init_constr */ \
        1., /* .init_constr_coef */ \
        1000, /* .num_samples */ \
        0, /* .samples_use_mutex */ \
        0, /* .samples_random_walk */ \
        0, /* .all_states_mutex_size */ \
    }

int pddlHPotInit(pddl_hpot_t *hpot,
                 const pddl_fdr_t *fdr,
                 const pddl_hpot_config_t *cfg,
                 bor_err_t *err);

void pddlHPotFree(pddl_hpot_t *hpot);

/**
 * Returns heuristic estimate for the given FDR state.
 */
int pddlHPotFDRStateEstimate(const pddl_hpot_t *hpot,
                             const pddl_fdr_vars_t *vars,
                             const int *state);

/**
 * Same as pddlHPotFDRStateEstimate() but no rounding is used.
 */
double pddlHPotFDRStateEstimateDbl(const pddl_hpot_t *hpot,
                                   const pddl_fdr_vars_t *vars,
                                   const int *state);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __PDDL_HPOT_H__ */
