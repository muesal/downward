/***
 * cpddl
 * -------
 * Copyright (c)2018 Daniel Fiser <danfis@danfis.cz>,
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
#include <boruvka/sort.h>
#include "pddl/fdr_app_op.h"
#include "pddl/causal_graph.h"
#include "assert.h"

/**
 * Base building structure for a tree node.
 */
struct pddl_fdr_app_op_tree {
    int var; /*!< Decision variable */
    bor_iset_t ops; /*!< List of immediate operators that are returned once
                         this node is reached */
    struct pddl_fdr_app_op_tree **val; /*!< Subtrees indexed by the value of
                                        the decision variable */
    int val_size;
    struct pddl_fdr_app_op_tree *def; /*!< Default subtree containing operators
                                       without precondition on the decision
                                       variable */
};
typedef struct pddl_fdr_app_op_tree pddl_fdr_app_op_tree_t;

/** Creates a new tree node (and recursively all subtrees) */
static pddl_fdr_app_op_tree_t *treeNew(const int *op_ids, int len,
                                   const int *var,
                                   const pddl_fdr_ops_t *ops);

static int opsSortCmp(const void *a, const void *b, void *ud)
{
    const pddl_fdr_app_op_t *s = ud;
    int op_id1 = *(const int *)a;
    int op_id2 = *(const int *)b;
    const pddl_fdr_op_t *opa = s->ops->op[op_id1];
    const pddl_fdr_op_t *opb = s->ops->op[op_id2];
    const pddl_fdr_part_state_t *prea = &opa->pre;
    const pddl_fdr_part_state_t *preb = &opb->pre;

    for (int i = 0; i < s->var_size; ++i){
        int var = s->var_order[i];
        int aset = pddlFDRPartStateIsSet(prea, var);
        int bset = pddlFDRPartStateIsSet(preb, var);

        if (aset && bset){
            int aval = pddlFDRPartStateGet(prea, var);
            int bval = pddlFDRPartStateGet(preb, var);

            if (aval < bval){
                return -1;
            }else if (aval > bval){
                return 1;
            }
        }else if (!aset && bset){
            return -1;
        }else if (aset && !bset){
            return 1;
        }
    }

    // make the sort stable
    if (op_id1 < op_id2)
        return -1;
    return 1;
}

static int *sortedOps(const pddl_fdr_app_op_t *app)
{
    int *op_ids;

    op_ids = BOR_ALLOC_ARR(int, app->ops->op_size);
    for (int i = 0; i < app->ops->op_size; ++i)
        op_ids[i] = i;

    borSort(op_ids, app->ops->op_size, sizeof(int), opsSortCmp, (void *)app);
    return op_ids;
}

static void treeBuildSetOps(pddl_fdr_app_op_tree_t *tree, const int *ops, int len)
{
    borISetEmpty(&tree->ops);
    for (int i = 0; i < len; ++i)
        borISetAdd(&tree->ops, ops[i]);
}

static int treeBuildDef(pddl_fdr_app_op_tree_t *tree,
                        const int *op_ids, int len,
                        const int *var,
                        const pddl_fdr_ops_t *ops)
{
    int size;

    for (size = 1;
         size < len
            && !pddlFDRPartStateIsSet(&ops->op[op_ids[size]]->pre, *var);
         ++size);

    tree->var = *var;
    tree->def = treeNew(op_ids, size, var + 1, ops);

    return size;
}

static void treeBuildPrepareVal(pddl_fdr_app_op_tree_t *tree, int val)
{
    int i;

    tree->val_size = val + 1;
    tree->val = BOR_ALLOC_ARR(pddl_fdr_app_op_tree_t *, tree->val_size);

    for (i = 0; i < tree->val_size; ++i)
        tree->val[i] = NULL;
}

static int treeBuildVal(pddl_fdr_app_op_tree_t *tree,
                        const int *op_ids, int len,
                        const int *var,
                        const pddl_fdr_ops_t *ops)
{
    int val = pddlFDRPartStateGet(&ops->op[op_ids[0]]->pre, *var);

    int size;
    for (size = 1;
         size < len
            && pddlFDRPartStateGet(&ops->op[op_ids[size]]->pre, *var) == val;
         ++size);

    tree->var = *var;
    tree->val[val] = treeNew(op_ids, size, var + 1, ops);

    return size;
}

static pddl_fdr_app_op_tree_t *treeNew(const int *op_ids, int len,
                                   const int *var,
                                   const pddl_fdr_ops_t *ops)
{
    pddl_fdr_app_op_tree_t *tree;
    const pddl_fdr_op_t *last_op = ops->op[op_ids[len - 1]];
    const pddl_fdr_part_state_t *last_pre = &last_op->pre;
    const pddl_fdr_op_t *first_op = ops->op[op_ids[0]];
    const pddl_fdr_part_state_t *first_pre = &first_op->pre;
    int start;

    tree = BOR_ALLOC(pddl_fdr_app_op_tree_t);
    tree->var = -1;
    borISetInit(&tree->ops);
    tree->val = NULL;
    tree->val_size = 0;
    tree->def = NULL;

    if (len == 0)
        return tree;

    // Find first variable that is set for at least one operator.
    // The operators are sorted so that it is enough to check the last
    // operator in the array.
    for (; *var != -1; ++var){
        if (pddlFDRPartStateIsSet(last_pre, *var))
            break;
    }

    if (*var == -1){
        // If there isn't any operator with set value anymore insert all
        // operators as immediate ops and exit.
        treeBuildSetOps(tree, op_ids, len);
        return tree;
    }

    // Now we know that array of operators contain at least one operator
    // with set value of current variable.

    // Prepare val array -- we now that the last operator in array has
    // largest value.
    treeBuildPrepareVal(tree, pddlFDRPartStateGet(last_pre, *var));

    // Initialize index of the first element with current value
    start = 0;

    // First check unset values from the beggining of the array
    if (!pddlFDRPartStateIsSet(first_pre, *var)){
        start = treeBuildDef(tree, op_ids, len, var, ops);
    }

    // Then build subtree for each value
    while (start < len){
        start += treeBuildVal(tree, op_ids + start, len - start, var, ops);
    }

    return tree;
}

static void treeDel(pddl_fdr_app_op_tree_t *tree)
{
    int i;

    borISetFree(&tree->ops);
    if (tree->val){
        for (i = 0; i < tree->val_size; ++i)
            if (tree->val[i])
                treeDel(tree->val[i]);
        BOR_FREE(tree->val);
    }

    if (tree->def)
        treeDel(tree->def);

    BOR_FREE(tree);
}

static void appOpInit(pddl_fdr_app_op_t *app,
                      const pddl_fdr_ops_t *ops,
                      int var_size,
                      const pddl_causal_graph_t *cg)
{
    int *sorted_ops = NULL;

    app->ops = ops;

    // Copy var_order to internal storage
    app->var_size = var_size;
    app->var_order = BOR_ALLOC_ARR(int, var_size + 1);
    memcpy(app->var_order, cg->var_order,
           sizeof(int) * cg->var_order_size);
    for (int ins = cg->var_order_size, i = 0;
            ins < var_size && i < var_size; ++i){
        if (!cg->important_var[i])
            app->var_order[ins++] = i;
    }
    app->var_order[var_size] = -1;

    if (ops->op_size > 0)
        sorted_ops = sortedOps(app);

    app->root = treeNew(sorted_ops, ops->op_size, app->var_order, ops);

    if (sorted_ops)
        BOR_FREE(sorted_ops);
}

void pddlFDRAppOpInit(pddl_fdr_app_op_t *app,
                      const pddl_fdr_vars_t *vars,
                      const pddl_fdr_ops_t *ops,
                      const pddl_fdr_part_state_t *goal)
{
    pddl_causal_graph_t *cg;
    cg = pddlCausalGraphNew(vars->var_size);
    pddlCausalGraphBuildFromOps(cg, ops);
    pddlCausalGraph(cg, goal);
    appOpInit(app, ops, vars->var_size, cg);
    pddlCausalGraphDel(cg);
}

void pddlFDRAppOpFree(pddl_fdr_app_op_t *app)
{
    if (app->root)
        treeDel(app->root);
    if (app->var_order != NULL)
        BOR_FREE(app->var_order);
}

static int treeFind(const pddl_fdr_app_op_tree_t *tree,
                    const int *vals,
                    bor_iset_t *ops)
{
    // insert all immediate operators
    borISetUnion(ops, &tree->ops);
    int found = borISetSize(&tree->ops);

    // check whether this node should check on any variable value
    if (tree->var != -1){
        // get corresponding value from state
        int val = vals[tree->var];

        // and use tree corresponding to the value if present
        if (val != -1
                && val < tree->val_size
                && tree->val[val]){
            found += treeFind(tree->val[val], vals, ops);
        }

        // use default tree if present
        if (tree->def)
            found += treeFind(tree->def, vals, ops);
    }

    return found;
}

int pddlFDRAppOpFind(const pddl_fdr_app_op_t *app,
                     const int *state,
                     bor_iset_t *ops)
{
    if (app->root == NULL)
        return 0;
    return treeFind(app->root, state, ops);
}
