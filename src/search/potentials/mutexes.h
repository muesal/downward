#ifndef DOWNWARD_MUTEXES_H
#define DOWNWARD_MUTEXES_H

#include "util.h"

#include "../option_parser.h"
#include "../task_utils/task_properties.h"
#include "../plugin.h"
#include "../global_state.h"
#include "../utils/rng.h"
#include "../utils/rng_options.h"
#include "../heuristic.h" // something here includes utils::log

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <zconf.h>


using namespace std;
using namespace potentials;

namespace options {
    class Options;
}

class MutexTable {

    using Pair = std::pair<FactPair, FactPair>;

private:
    vector<Pair> mutexes;
    VariablesProxy variables;
    TaskProxy task_proxy;
    std::map<Pair, int> hm_table;

    void generate_all_pairs(); // recursively generate all tuples.
    static vector<Pair> generate_all_pairs(vector<FactPair> &tuples); // generate all pairs.

    void init_hm_table(vector<FactPair> &t);

    void update_hm_table();

    bool all_reachable(vector<FactPair> &t);

    static bool unassigned(map<int, int> &state, int variable_id);

    static vector<int> set_minus(const vector<int> &domain, const vector<FactPair> &mutex, int variable_id);

    static vector<FactPair> intersection(const vector<FactPair> &one, const vector<FactPair> &two);

    vector<FactPair> get_mutex_with_state(map<int, int> &state);

    void get_mutex_with_fact(int variable, int value, vector<FactPair> &mf);

public:
    explicit MutexTable(TaskProxy task_proxy);

    vector<vector<int>>
    multi_fact_disambiguation(map<int, int> &state);

    const VariablesProxy *getVariablesProxy() const;

    ~MutexTable() = default;
};

#endif //DOWNWARD_MUTEXES_H
