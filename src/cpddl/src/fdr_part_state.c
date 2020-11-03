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
#include "pddl/fdr_part_state.h"
#include "assert.h"

void pddlFDRPartStateInit(pddl_fdr_part_state_t *ps)
{
    bzero(ps, sizeof(*ps));
}

void pddlFDRPartStateInitCopy(pddl_fdr_part_state_t *dst,
                              const pddl_fdr_part_state_t *src)
{
    dst->fact_size = src->fact_size;
    dst->fact_alloc = src->fact_alloc;
    ASSERT(dst->fact_alloc >= dst->fact_size);
    dst->fact = BOR_ALLOC_ARR(pddl_fdr_fact_t, src->fact_alloc);
    memcpy(dst->fact, src->fact, sizeof(*dst->fact) * dst->fact_size);
}

void pddlFDRPartStateFree(pddl_fdr_part_state_t *ps)
{
    if (ps->fact != NULL)
        BOR_FREE(ps->fact);
}

void pddlFDRPartStateSet(pddl_fdr_part_state_t *ps, int var, int val)
{
    for (int i = 0; i < ps->fact_size; ++i){
        if (ps->fact[i].var == var){
            ps->fact[i].val = val;
            return;
        }
    }

    if (ps->fact_alloc == ps->fact_size){
        if (ps->fact_alloc == 0)
            ps->fact_alloc = 1;
        ps->fact_alloc *= 2;
        ps->fact = BOR_REALLOC_ARR(ps->fact, pddl_fdr_fact_t, ps->fact_alloc);
    }

    ps->fact[ps->fact_size].var = var;
    ps->fact[ps->fact_size].val = val;
    ++ps->fact_size;
    for (int i = ps->fact_size - 1; i > 0; --i){
        if (ps->fact[i - 1].var > var){
            ps->fact[i].var = ps->fact[i - 1].var;
            ps->fact[i].val = ps->fact[i - 1].val;
            ps->fact[i - 1].var = var;
            ps->fact[i - 1].val = val;
        }
    }
}

int pddlFDRPartStateGet(const pddl_fdr_part_state_t *ps, int var)
{
    for (int i = 0; i < ps->fact_size && i <= var; ++i){
        if (ps->fact[i].var == var)
            return ps->fact[i].val;
    }
    return -1;
}

int pddlFDRPartStateIsSet(const pddl_fdr_part_state_t *ps, int var)
{
    return pddlFDRPartStateGet(ps, var) >= 0;
}

int pddlFDRPartStateCmp(const pddl_fdr_part_state_t *p1,
                        const pddl_fdr_part_state_t *p2)
{
    int cmp = p1->fact_size - p2->fact_size;
    for (int i = 0; i < p1->fact_size && cmp == 0; ++i)
        cmp = p1->fact[i].var - p2->fact[i].var;
    return cmp;
}

void pddlFDRPartStateRemapFacts(pddl_fdr_part_state_t *ps,
                                const pddl_fdr_vars_remap_t *remap)
{
    int ins = 0;
    for (int i = 0; i < ps->fact_size; ++i){
        pddl_fdr_fact_t *fact = ps->fact + i;
        if (remap->remap[fact->var][fact->val] != NULL){
            const pddl_fdr_val_t *v = remap->remap[fact->var][fact->val];
            ASSERT(v != NULL);
            fact->var = v->var_id;
            fact->val = v->val_id;
            ps->fact[ins++] = *fact;
        }
    }
    ps->fact_size = ins;
}
