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

#include <boruvka/alloc.h>
#include <boruvka/pairheap.h>

#include "pddl/causal_graph.h"

/** One strongly connected component */
struct scc_comp {
    int *var;     /*!< Variables in component */
    int var_size; /*!< Number of variables in component */
};
typedef struct scc_comp scc_comp_t;

/** Strongly connected components */
struct scc {
    scc_comp_t *comp; /*!< List of components */
    int comp_size;    /*!< Number of components */
};
typedef struct scc scc_t;

static void graphInit(pddl_causal_graph_graph_t *g, int var_size);
static void graphFree(pddl_causal_graph_graph_t *g);
/** Copies src graph into dst with unimportant variables */
static void graphCopyImportant(pddl_causal_graph_graph_t *dst,
                               const pddl_causal_graph_graph_t *src,
                               const int *important_var);
/** Removes from graph edges that are not within found strongly connected
 *  components. */
static void graphPruneBySCC(pddl_causal_graph_graph_t *graph,
                            const scc_t *scc);
/** Fills .important_var array with 0/1 signaling whether there is
 *  connection between the variable and a goal. */
static void markImportantVars(pddl_causal_graph_t *cg,
                              const pddl_fdr_part_state_t *goal);
/** Creates an ordering of variables based on the given graph */
static void createOrdering(pddl_causal_graph_t *cg);
/** Determines strongly connected components */
static scc_t *sccNew(const pddl_causal_graph_graph_t *graph, int var_size);
static void sccDel(scc_t *);

pddl_causal_graph_t *pddlCausalGraphNew(int var_size)
{
    pddl_causal_graph_t *cg;

    cg = BOR_ALLOC(pddl_causal_graph_t);
    cg->var_size = var_size;
    graphInit(&cg->successor_graph, var_size);
    graphInit(&cg->predecessor_graph, var_size);
    cg->important_var = BOR_ALLOC_ARR(int, cg->var_size);
    cg->var_order = BOR_ALLOC_ARR(int, cg->var_size + 1);

    return cg;
}

void pddlCausalGraphDel(pddl_causal_graph_t *cg)
{
    graphFree(&cg->successor_graph);
    graphFree(&cg->predecessor_graph);
    if (cg->important_var)
        BOR_FREE(cg->important_var);
    if (cg->var_order)
        BOR_FREE(cg->var_order);
    BOR_FREE(cg);
}

void pddlCausalGraph(pddl_causal_graph_t *cg,
                     const pddl_fdr_part_state_t *goal)
{
    // Set up .important_var[]
    markImportantVars(cg, goal);

    // Compute ordering of variables (.var_order*)
    createOrdering(cg);
}

static void graphInit(pddl_causal_graph_graph_t *g, int var_size)
{
    g->var_size = var_size;
    g->edge_size = BOR_CALLOC_ARR(int, g->var_size);
    g->value     = BOR_CALLOC_ARR(int *, g->var_size);
    g->end_var   = BOR_CALLOC_ARR(int *, g->var_size);
}

static void graphFree(pddl_causal_graph_graph_t *g)
{
    int i;

    for (i = 0; i < g->var_size; ++i){
        if (g->value[i])
            BOR_FREE(g->value[i]);
        if (g->end_var[i])
            BOR_FREE(g->end_var[i]);
    }

    BOR_FREE(g->edge_size);
    BOR_FREE(g->value);
    BOR_FREE(g->end_var);
}

static void graphCopyImportant(pddl_causal_graph_graph_t *dst,
                               const pddl_causal_graph_graph_t *src,
                               const int *important_var)
{
    int var, i, edge_size, *value, *end_var, *src_value, *src_end_var, ins, w;

    for (var = 0; var < src->var_size; ++var){
        // skip unimportant vars
        if (!important_var[var])
            continue;

        edge_size = dst->edge_size[var] = src->edge_size[var];
        value = dst->value[var] = BOR_ALLOC_ARR(int, edge_size);
        end_var = dst->end_var[var] = BOR_ALLOC_ARR(int, edge_size);

        src_value = src->value[var];
        src_end_var = src->end_var[var];

        for (i = 0, ins = 0; i < edge_size; ++i){
            w = src_end_var[i];
            if (!important_var[w])
                continue;

            end_var[ins] = w;
            value[ins++] = src_value[i];
        }
        dst->edge_size[var] = ins;
    }
}

static void graphPruneVarBySCC(pddl_causal_graph_graph_t *graph, int var,
                               const scc_comp_t *comp)
{
    int gi, ci, gv, cv, edge_size, *value, *end_var, ins;

    if (comp->var_size == 1){
        graph->edge_size[var] = 0;
    }

    edge_size = graph->edge_size[var];
    value = graph->value[var];
    end_var = graph->end_var[var];

    // Both arrays in comp and in end_var are sorted
    ins = 0;
    for (gi = ci = 0; gi < edge_size && ci < comp->var_size;){
        gv = end_var[gi];
        cv = comp->var[ci];

        if (gv == cv){
            // Keep this edge
            if (ins != gi){
                value[ins] = value[gi];
                end_var[ins] = end_var[gi];
            }
            ++ins;
            ++gi;
            ++ci;

        }else if (gv < cv){
            ++gi;
        }else{ // (gv > cv)
            ++ci;
        }
    }

    graph->edge_size[var] = ins;
}

static void graphPruneBySCC(pddl_causal_graph_graph_t *graph,
                            const scc_t *scc)
{
    int ci, v;
    const scc_comp_t *comp;

    for (ci = 0; ci < scc->comp_size; ++ci){
        comp = scc->comp + ci;
        for (v = 0; v < comp->var_size; ++v){
            graphPruneVarBySCC(graph, comp->var[v], comp);
        }
    }
}



/** Edge from variable to variable in rbtree-based graph */
struct _graph_edge_t {
    bor_rbtree_int_node_t rbtree;
    int value;
};
typedef struct _graph_edge_t graph_edge_t;

/** Variable node in rbtree-based graph */
struct _graph_node_t {
    bor_rbtree_int_node_t rbtree;
    bor_rbtree_int_t *edge;
    int edge_size;
};
typedef struct _graph_node_t graph_node_t;

/** Converts rbtree-based graph to array based graph */
static void rbtreeToGraph(bor_rbtree_int_t *rb,
                          pddl_causal_graph_graph_t *graph)
{
    bor_rbtree_int_node_t *tree_node, *tree_edge;
    graph_node_t *node;
    graph_edge_t *edge;
    int var, i, edge_size, *value = NULL, *end_var = NULL;

    while (!borRBTreeIntEmpty(rb)){
        tree_node = borRBTreeIntExtractMin(rb);
        node = bor_container_of(tree_node, graph_node_t, rbtree);
        var = borRBTreeIntKey(tree_node);

        edge_size = graph->edge_size[var] = node->edge_size;
        if (edge_size > 0){
            value = graph->value[var] = BOR_ALLOC_ARR(int, edge_size);
            end_var = graph->end_var[var] = BOR_ALLOC_ARR(int, edge_size);
        }

        for (i = 0; !borRBTreeIntEmpty(node->edge); ++i){
            tree_edge = borRBTreeIntExtractMin(node->edge);
            edge = bor_container_of(tree_edge, graph_edge_t, rbtree);

            end_var[i] = borRBTreeIntKey(tree_edge);
            value[i]   = edge->value;

            BOR_FREE(edge);
        }

        borRBTreeIntDel(node->edge);
        BOR_FREE(node);
    }
}

/** Adds one connection to the graph (or increases edge's value by one if
 *  the connection is already there). */
static void graphAddConnection(bor_rbtree_int_t *graph, int from, int to)
{
    bor_rbtree_int_node_t *tree_node, *tree_edge;
    graph_node_t *node;
    graph_edge_t *edge;

    if (from == to)
        return;

    tree_node = borRBTreeIntFind(graph, from);
    if (tree_node != NULL){
        node = bor_container_of(tree_node, graph_node_t, rbtree);
    }else{
        node = BOR_ALLOC(graph_node_t);
        node->edge = borRBTreeIntNew();
        node->edge_size = 0;
        borRBTreeIntInsert(graph, from, &node->rbtree);
    }

    tree_edge = borRBTreeIntFind(node->edge, to);
    if (tree_edge != NULL){
        edge = bor_container_of(tree_edge, graph_edge_t, rbtree);
    }else{
        edge = BOR_ALLOC(graph_edge_t);
        edge->value = 0;
        borRBTreeIntInsert(node->edge, to, &edge->rbtree);
        ++node->edge_size;
    }

    ++edge->value;
}

void pddlCausalGraphBuildInit(pddl_causal_graph_build_t *cg_build)
{
    cg_build->succ_graph = borRBTreeIntNew();
    cg_build->pred_graph = borRBTreeIntNew();
}

void pddlCausalGraphBuildFree(pddl_causal_graph_build_t *cg_build)
{
    borRBTreeIntDel(cg_build->succ_graph);
    borRBTreeIntDel(cg_build->pred_graph);
}

void pddlCausalGraphBuildAdd(pddl_causal_graph_build_t *cg_build,
                             int pre_var, int eff_var)
{
    graphAddConnection(cg_build->succ_graph, pre_var, eff_var);
    graphAddConnection(cg_build->pred_graph, eff_var, pre_var);
}

void pddlCausalGraphBuildFromOps(pddl_causal_graph_t *cg,
                                 const pddl_fdr_ops_t *ops)
{
    pddl_causal_graph_build_t build;

    pddlCausalGraphBuildInit(&build);

    for (int i = 0; i < ops->op_size; ++i){
        const pddl_fdr_op_t *op = ops->op[i];

        for (int prei = 0; prei < op->pre.fact_size; ++prei){
            int pre_var = op->pre.fact[prei].var;
            for (int effi = 0; effi < op->eff.fact_size; ++effi){
                int eff_var = op->eff.fact[effi].var;
                pddlCausalGraphBuildAdd(&build, pre_var, eff_var);
            }
        }
        /*
        // TODO:
        for (int effi = 0; effi < op->eff.fact_size; ++effi){
            int eff_var = op->eff.fact[effi].var;
            for (int effi2 = 0; effi2 < op->eff.fact_size; ++effi2){
                int eff_var2 = op->eff.fact[effi2].var;
                pddlCausalGraphBuildAdd(&build, eff_var, eff_var2);
            }
        }
        */
    }

    pddlCausalGraphBuild(cg, &build);
    pddlCausalGraphBuildFree(&build);
}

void pddlCausalGraphBuild(pddl_causal_graph_t *cg,
                          pddl_causal_graph_build_t *cg_build)
{
    rbtreeToGraph(cg_build->succ_graph, &cg->successor_graph);
    rbtreeToGraph(cg_build->pred_graph, &cg->predecessor_graph);
}


static void markImportantVarsDFS(pddl_causal_graph_t *cg, int var)
{
    int i, len, *other_var, w;

    len = cg->predecessor_graph.edge_size[var];
    other_var = cg->predecessor_graph.end_var[var];
    for (i = 0; i < len; ++i){
        w = other_var[i];
        if (!cg->important_var[w]){
            cg->important_var[w] = 1;
            markImportantVarsDFS(cg, w);
        }
    }
}

static void markImportantVars(pddl_causal_graph_t *cg,
                              const pddl_fdr_part_state_t *goal)
{
    // Set all variables as unimportant
    bzero(cg->important_var, sizeof(int) * cg->var_size);

    for (int fi = 0; fi < goal->fact_size; ++fi){
        int var = goal->fact[fi].var;
        if (!cg->important_var[var]){
            cg->important_var[var] = 1;
            markImportantVarsDFS(cg, var);
        }
    }
}

struct _order_var_t {
    bor_pairheap_node_t heap;
    int var; /*!< ID of the variable */
    int w;   /*!< Incoming weight */
    int ins; /*!< True if inserted into heap */
};
typedef struct _order_var_t order_var_t;

static int heapLT(const bor_pairheap_node_t *a,
                  const bor_pairheap_node_t *b, void *_)
{
    order_var_t *v1 = bor_container_of(a, order_var_t, heap);
    order_var_t *v2 = bor_container_of(b, order_var_t, heap);
    if (v1->w == v2->w)
        return v1->var < v2->var;
    return v1->w < v2->w;
}

/** Computes sum of incoming weights for each variable */
static void computeIncomingWeights(const pddl_causal_graph_graph_t *graph,
                                   order_var_t *var)
{
    int i, j, edge_size, *value, *end_var;

    for (i = 0; i < graph->var_size; ++i){
        edge_size = graph->edge_size[i];
        value = graph->value[i];
        end_var = graph->end_var[i];
        for (j = 0; j < edge_size; ++j){
            var[end_var[j]].w += value[j];
        }
    }
}

/** Updates all var's successors' incoming weights and this their position
 *  in the heap. */
static void orderingUpdateSucc(bor_pairheap_t *heap,
                               order_var_t *order_var,
                               const pddl_causal_graph_graph_t *graph,
                               int var, int w)
{
    int i, edge_size, *end_var, v;

    edge_size = graph->edge_size[var];
    end_var   = graph->end_var[var];
    for (i = 0; i < edge_size; ++i){
        v = end_var[i];
        if (order_var[v].ins){
            order_var[v].w -= w;
            borPairHeapDecreaseKey(heap, &order_var[v].heap);
        }
    }
}

static void updateOrdering(pddl_causal_graph_t *cg,
                           order_var_t *order_var,
                           const pddl_causal_graph_graph_t *graph)
{
    bor_pairheap_t *heap;
    bor_pairheap_node_t *hnode;
    order_var_t *minvar;
    int i;

    // Create and initialize heap
    heap = borPairHeapNew(heapLT, NULL);
    for (i = 0; i < cg->var_size; ++i){
        if (!cg->important_var[i])
            continue;

        borPairHeapAdd(heap, &order_var[i].heap);
        order_var[i].ins = 1;
    }

    // Order variables as DAG
    while (!borPairHeapEmpty(heap)){
        // Remove node with lowest incoming weight
        hnode = borPairHeapExtractMin(heap);
        minvar = bor_container_of(hnode, order_var_t, heap);
        minvar->ins = 0;

        // Add variable to the ordering array
        cg->var_order[cg->var_order_size++] = minvar->var;

        // Update all successor variables
        orderingUpdateSucc(heap, order_var, graph, minvar->var, minvar->w);
    }

    borPairHeapDel(heap);
}

static void reverseVarOrder(int *arr, int size)
{
    int i, len, tmp;

    len = size / 2;
    for (i = 0; i < len; ++i){
        BOR_SWAP(arr[i], arr[size - 1 - i], tmp);
    }
}

static void createOrdering(pddl_causal_graph_t *cg)
{
    pddl_causal_graph_graph_t scc_graph;
    scc_t *scc;
    order_var_t *var;
    int i;

    // Create a copy of successor graph without unimportant variables
    graphInit(&scc_graph, cg->var_size);
    graphCopyImportant(&scc_graph, &cg->successor_graph,
                       cg->important_var);

    // Compute strongly connected components
    scc = sccNew(&scc_graph, cg->var_size);

    // Remove edges outside strongly connected components
    graphPruneBySCC(&scc_graph, scc);

    // Initialize array with variables
    var = BOR_ALLOC_ARR(order_var_t, cg->var_size);
    for (i = 0; i < cg->var_size; ++i){
        var[i].var = i;
        var[i].w = 0;
        var[i].ins = 0;
    }

    // Compute sum of all incoming weights for all variables
    computeIncomingWeights(&scc_graph, var);

    // Create ordering from components
    cg->var_order_size = 0;
    updateOrdering(cg, var, &scc_graph);
    reverseVarOrder(cg->var_order, cg->var_order_size);
    cg->var_order[cg->var_order_size] = -1;

    BOR_FREE(var);
    sccDel(scc);
    graphFree(&scc_graph);
}


/** Context for DFS during computing SCC */
struct scc_dfs {
    const pddl_causal_graph_graph_t *graph;
    int cur_index;
    int *index;
    int *lowlink;
    int *in_stack;
    int *stack;
    int stack_size;
};
typedef struct scc_dfs scc_dfs_t;

static int cmpInt(const void *a, const void *b)
{
    return (*(int *)a) - (*(int *)b);
}

static void sccTarjanStrongconnect(scc_t *scc, scc_dfs_t *dfs, int var)
{
    int i, len, *end_var, w;
    scc_comp_t *comp;

    dfs->index[var] = dfs->lowlink[var] = dfs->cur_index++;
    dfs->stack[dfs->stack_size++] = var;
    dfs->in_stack[var] = 1;

    len     = dfs->graph->edge_size[var];
    end_var = dfs->graph->end_var[var];
    for (i = 0; i < len; ++i){
        w = end_var[i];
        if (dfs->index[w] == -1){
            sccTarjanStrongconnect(scc, dfs, w);
            dfs->lowlink[var] = BOR_MIN(dfs->lowlink[var], dfs->lowlink[w]);
        }else if (dfs->in_stack[w]){
            dfs->lowlink[var] = BOR_MIN(dfs->lowlink[var], dfs->lowlink[w]);
        }
    }

    if (dfs->index[var] == dfs->lowlink[var]){
        // Find how deep unroll stack
        for (i = dfs->stack_size - 1; dfs->stack[i] != var; --i)
            dfs->in_stack[dfs->stack[i]] = 0;
        dfs->in_stack[dfs->stack[i]] = 0;

        // Create new component
        ++scc->comp_size;
        scc->comp = BOR_REALLOC_ARR(scc->comp, scc_comp_t, scc->comp_size);
        comp = scc->comp + scc->comp_size - 1;
        comp->var_size = dfs->stack_size - i;
        comp->var = BOR_ALLOC_ARR(int, comp->var_size);

        // Copy variable IDs from stack to the component
        memcpy(comp->var, dfs->stack + i, sizeof(int) * comp->var_size);
        if (comp->var_size > 1)
            qsort(comp->var, comp->var_size, sizeof(int), cmpInt);

        // Shrink stack
        dfs->stack_size = i;
    }
}

static void sccTarjan(scc_t *scc, const pddl_causal_graph_graph_t *graph,
                      int var_size)
{
    scc_dfs_t dfs;
    int i, var;

    // Initialize structure for Tarjan's algorithm
    dfs.graph = graph;
    dfs.cur_index = 0;
    dfs.index    = BOR_ALLOC_ARR(int, 4 * var_size);
    dfs.lowlink  = dfs.index + var_size;
    dfs.in_stack = dfs.lowlink + var_size;
    dfs.stack    = dfs.in_stack + var_size;
    dfs.stack_size = 0;
    for (i = 0; i < var_size; ++i){
        dfs.index[i] = dfs.lowlink[i] = -1;
        dfs.in_stack[i] = 0;
    }

    for (var = 0; var < var_size; ++var){
        if (dfs.index[var] == -1)
            sccTarjanStrongconnect(scc, &dfs, var);
    }

    BOR_FREE(dfs.index);
}

static scc_t *sccNew(const pddl_causal_graph_graph_t *graph, int var_size)
{
    scc_t *scc;

    scc = BOR_ALLOC(scc_t);
    scc->comp_size = 0;
    scc->comp = NULL;

    // Run Tarjan's algorithm for finding strongly connected components.
    sccTarjan(scc, graph, var_size);

    return scc;
}

static void sccDel(scc_t *scc)
{
    int i;

    for (i = 0; i < scc->comp_size; ++i){
        if (scc->comp[i].var != NULL)
            BOR_FREE(scc->comp[i].var);
    }
    if (scc->comp)
        BOR_FREE(scc->comp);
    BOR_FREE(scc);
}
