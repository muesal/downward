//
// Created by salome on 26/09/2020.
//

#ifndef DOWNWARD_MUTEXES_H
#define DOWNWARD_MUTEXES_H

#include "util.h"

#include "../option_parser.h"
#include "../plugin.h"
#include "../global_state.h"
#include "../utils/rng.h"
#include "../utils/rng_options.h"


#include "../heuristics/hm_heuristic.h"
#include <iostream>
#include <fstream>
#include <zconf.h>


using namespace std;
using namespace hm_heuristic;

namespace options {
    class Options;
}

class MutexTable {

public:
    vector<vector<FactPair>> mutexes;
    VariablesProxy &variables;

    explicit MutexTable(Options &opts, VariablesProxy &variables, State &state);

    static bool unassigned(map<int, int> &state, int variable_id);

    static vector<int> set_minus(const vector<int> &domain, const vector<FactPair> &mutex, int variable_id);

    static vector<FactPair> intersection(const vector<FactPair> &one, const vector<FactPair> &two);

    vector<FactPair> get_mutex_with_state(map<int, int> &state);

    void get_mutex_with_fact(int variable, int value, vector<FactPair> &mf);

    vector<vector<int>>
    multi_fact_disambiguation(map<int, int> &state);

    ~MutexTable() = default;

    VariablesProxy getVariableProxy() const;
};


#endif //DOWNWARD_MUTEXES_H
