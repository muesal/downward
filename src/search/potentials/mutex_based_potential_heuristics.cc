//
// Created by salome on 02.08.20.
//
#include "potential_function.h"
#include "potential_heuristic.h"
#include "potential_optimizer.h"
#include "util.h"

 #include "../heuristics/hm_heuristic.h"

#include "../option_parser.h"
#include "../plugin.h"
#include "../global_state.h"
#include "../utils/rng.h"
#include "../utils/rng_options.h"

#include <memory>
#include <vector>

using namespace std;
using namespace hm_heuristic;
using Tuple = vector<FactPair>;

namespace potentials {

    // take the domain and FactPairs which are mutex to the current state, return the remaining possible vlues for this variable.
    vector<int> set_minus(vector<int> domain, vector<FactPair> mutex, int variable_id) {
        vector<int> set_minus;
        for (FactPair m : mutex) {
            for (auto i = domain.begin(); i < domain.end(); i++) {
                if (m.var == variable_id && m.value == *i){
                    domain.erase(i);
                }
            }
        }
        return set_minus;
    }

    // Overload, to have method call wihthout argument mf?
    // or to have one which only takes elements in mf?
    // save all FactPairs which are mmutex to the fact <variable, value>
    static void get_mutex_with_fact(int variable, int value, vector<Tuple> &mutexes, vector<FactPair> &mf) {
        for (Tuple mutex : mutexes) {
            if (mutex[0].var == variable && mutex[0].value == value){
                mf.push_back(mutex[1]);
            } else if (mutex[1].var == variable && mutex[1].value == value){
                mf.push_back(mutex[0]);
            }
        }
    }

    // Save all mutex pairs which are unreachable for this (partial) state
    static vector<FactPair> get_mutex_with_state(State &state, vector<Tuple> &mutexes) {
        vector<FactPair> mf;
        for (Tuple mutex : mutexes) {
            if (state[mutex[0].var].get_value() == mutex[0].value){
                mf.push_back(mutex[1]);
            } else if (state[mutex[1].var].get_value() == mutex[1].value){
                mf.push_back(mutex[0]);
            }
        }
        return mf;
    }

    /**
     * TODO: The following probably have the first 11 lines in common...
     */

    static bool single_fact_diasambiguation(State &state, VariablesProxy &variables, vector<Tuple> &mutexes){
        //get all domains and id's of the variables
        vector<vector<int>> domains;
        vector<int> variable_id;
        for(VariableProxy v : variables) {
            vector<int> dom;
            for (int d = 0; d < v.get_domain_size(); d++){
                dom.push_back(d);
            }
            domains.push_back(dom);
            variable_id.push_back(v.get_id());
        }
        vector<FactPair> mp = get_mutex_with_state(state, mutexes);

        //while state changes check for single fact disambiguations and assign them
        bool changed = true;
        while (changed) {
            changed = false;
            for(size_t v = 0; v < state.size(); v++) {
                if (!state.is_defined(v)){
                    //algorithm 1 l. 4 (D_V <- F_V \ M_p)
                    domains[v] = set_minus(domains[v], mp, variable_id[v]);

                    //algorithm 1 l. 6
                    if (domains[v].size() == 0) {
                        return false;
                    }

                    //algorithm 1 l. 5
                    if (domains[v].size() == 1){
                        get_mutex_with_fact(variable_id[v], domains[v][0], mutexes, mp);
                        changed = true;
                    }
                }
            }
        }
        return true;
    }

    static vector<vector<int>> multi_fact_diasambiguation(State &state, VariablesProxy &variables, vector<Tuple> &mutexes){
        //get all domains and id's of the variables
        vector<vector<int>> domains;
        vector<int> variable_id;
        for(VariableProxy v : variables) {
            vector<int> dom;
            for (int d = 0; d < v.get_domain_size(); d++){
                dom.push_back(d);
            }
            domains.push_back(dom);
            variable_id.push_back(v.get_id());
        }
        vector<FactPair> a = get_mutex_with_state(state, mutexes);

        //while state changes check for single fact disambiguations and assign them
        bool changed = true;
        while (changed) {
            changed = false;
            for(size_t v = 0; v < state.size(); v++) {
                if (!state.is_defined(v)){
                    auto size = domains[v].size();
                    //algorithm 2 l. 7 (D_V <- D_V \ a)
                    domains[v] = set_minus(domains[v], a, variable_id[v]);

                    //algorithm 2 l. 8
                    if (domains[v].size() < size) {
                        vector<FactPair> mf;
                        get_mutex_with_fact(variable_id[v], domains[v][0], mutexes, mf);
                        for (size_t f = 1; f < domains[v].size(); f++) {
                            //TODO: this does append all mutex facts, instead of deleteing the ones, which should not be taken... better solution, then appending to mf?
                            get_mutex_with_fact(variable_id[v], domains[v][f], mutexes, mf);
                            if(mf.size() == 0) {
                                break;
                            }
                        }
                        if(mf.size() > 0) {
                            a.insert(a.begin(), mf.begin(), mf.end());
                        }
                        changed = true;
                    }
                }
            }
        }
        return domains;
    }

    static unique_ptr<PotentialFunction> create_mutex_based_potential_function(
            const Options &opts) {

        PotentialOptimizer optimizer(opts);
        const AbstractTask &task = *opts.get<shared_ptr<AbstractTask>>("transform");
        TaskProxy task_proxy(task);
        VariablesProxy variables = task_proxy.get_variables();
        
        State initial = task_proxy.get_initial_state();

        //get all mutexes
        auto hm = make_shared<HMHeuristic>(opts);
        std::vector<Tuple> mutexes = hm->get_unreachable_tuples(initial);
        hm.reset();

        single_fact_diasambiguation(initial, variables, mutexes);
        multi_fact_diasambiguation(initial, variables, mutexes);

        // optimize the optimizer for the initial state
        optimizer.optimize_for_state(initial);
        return optimizer.get_potential_function();
    }

    static shared_ptr<Heuristic> _parse(OptionParser &parser) {
        parser.document_synopsis(
                "Mutex-based potential heuristics",
                get_admissible_potentials_reference());
        parser.add_option<int>("m", "use h2 heuristic", "2", Bounds("0", "infinity"));
        prepare_parser_for_admissible_potentials(parser);
        Options opts = parser.parse();
        if (parser.dry_run())
            return nullptr;

        return make_shared<PotentialHeuristic>(
                opts, create_mutex_based_potential_function(opts));
    }

static Plugin<Evaluator> _plugin(
            "mutex_based_potential", _parse, "heuristics_potentials");
}