/***
 * cpddl
 * -------
 * Copyright (c)2015 Daniel Fiser <danfis@danfis.cz>,
 * Agent Technology Center, Department of Computer Science,
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

#ifndef __PDDL_CAUSALGRAPH_H__
#define __PDDL_CAUSALGRAPH_H__

#include <boruvka/rbtree_int.h>
#include <pddl/fdr_op.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

// TODO: Rewrite this whole module and get rid of the ugly "build" API

struct pddl_causal_graph_graph {
    int **end_var;  /*!< ID of variable at the end of the edge */
    int **value;    /*!< Value of the edge */
    int *edge_size; /*!< Number edges emanating from i'th variable */
    int var_size;   /*!< Number of variables */
};
typedef struct pddl_causal_graph_graph pddl_causal_graph_graph_t;

struct pddl_causal_graph {
    int var_size;

    /** Graph with edges from precondition vars to effect vars */
    pddl_causal_graph_graph_t successor_graph;
    /** Graph with edges from effect vars to precondition vars */
    pddl_causal_graph_graph_t predecessor_graph;
    /** Bool flag for each variable whether it is important, i.e., if it is
     *  connected through operators with the goal. */
    int *important_var;
    int *var_order; /*!< Ordered array of variable IDs, array ends with -1 */
    int var_order_size;       /*!< Number of elements in .var_order[] */
};
typedef struct pddl_causal_graph pddl_causal_graph_t;

struct pddl_causal_graph_build {
    bor_rbtree_int_t *succ_graph;
    bor_rbtree_int_t *pred_graph;
};
typedef struct pddl_causal_graph_build pddl_causal_graph_build_t;

/**
 * Initializes build structure.
 */
void pddlCausalGraphBuildInit(pddl_causal_graph_build_t *cg_build);

/**
 * Frees allocated resources of build structure.
 */
void pddlCausalGraphBuildFree(pddl_causal_graph_build_t *cg_build);

/**
 * Adds connection from variable in action precondition to a variable in
 * its effect.
 */
void pddlCausalGraphBuildAdd(pddl_causal_graph_build_t *cg_build,
                             int var_pre, int var_eff);

/**
 * Creates a new causal graph of variables
 */
pddl_causal_graph_t *pddlCausalGraphNew(int var_size);

/**
 * Deletes causal graph object.
 */
void pddlCausalGraphDel(pddl_causal_graph_t *cg);

/**
 * Builds causal graph from the build structure.
 */
void pddlCausalGraphBuild(pddl_causal_graph_t *cg,
                          pddl_causal_graph_build_t *cg_build);

/**
 * Fill causal graphs from operators.
 */
void pddlCausalGraphBuildFromOps(pddl_causal_graph_t *cg,
                                 const pddl_fdr_ops_t *ops);

/**
 * Process causal graph and fills .important_var and .var_order* members.
 */
void pddlCausalGraph(pddl_causal_graph_t *cg,
                     const pddl_fdr_part_state_t *goal);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __PDDL_CAUSALGRAPH_H__ */
