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
    
    // save all FactPairs which are mutex to the fact <variable, value>
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

    static vector<FactPair> intersection(vector<FactPair> one, vector<FactPair> two) {
        vector<FactPair> three;
        for (FactPair fact : one) {
            for (FactPair f : two) {
                if (f == fact) {
                    three.push_back(fact);
                    break;
                }
            }
        }
        return three;
    }

    /**
     * Parameters:
     *  state: partial state
     *  k: amount of remaining variables to assign
     *  f_v: position of the first assigned variable
     *  domains: domains of all variables
     */
    static vector<State> get_all_extensions(State &s, int k, int f_v, vector<vector<int>> &domains){
        vector<State> states;
        return get_all_extensions(s, k, f_v, domains, 0, states);
    }

    /**
     * Parameters:
     *  state: partial state
     *  k: amount of remaining variables to assign
     *  f_v: position of the first assigned variable
     *  domains: domains of all variables
     *  last_assigned_variable: last variable which was assigned
     *  extended: list of all states enough extended states
     */
    static vector<State> get_all_extensions(State &s, int k, int f_v, vector<vector<int>> &domains, int last_assigned_variable, vector<State> &extended){
        if (k == 0) {
            extended.push_back(s);
            return extended;
        }
        vector<State> ex;
        // variables assigned from beginning to end, therefore enough unassigned variables need to remain
        // size - k if only not assigned values left, size - k - 1 otherwise
        int before = f_v < last_assigned_variable? 0 : 1;
        for (int i = last_assigned_variable; i < s.size() - k - before; i++) {
            for (int d : domains[i]) {
                State e = State(s);
                e.add_variable(i, d);
                ex = get_all_extensions(e, k-1, f_v, domains, i + 1, extended); //TODO: nur aufrufe, falls k > 0, sonst e anhängen? für bessere performance...
                extended.insert(extended.end(), ex.begin(), ex.end());
                ex.clear();
            }
        }
        return extended;
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
            // TODO: for assigned values, add all values or only the assigned one?
            if ( !state.is_defined(v.get_id())) {
                for (int d = 0; d < v.get_domain_size(); d++){
                    dom.push_back(d);
                }
            } else {
                dom.push_back(state[v.get_id()].get_value());
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
                        vector<FactPair> mf, mf2;
                        get_mutex_with_fact(variable_id[v], domains[v][0], mutexes, mf);
                        for (size_t f = 1; f < domains[v].size(); f++) {
                            get_mutex_with_fact(variable_id[v], domains[v][f], mutexes, mf2);
                            mf = intersection(mf, mf2);
                            mf2.clear();
                            if(mf2.size() == 0 || mf.size() == 0) {
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

    /**
     * Parameters:
     *  state: partial state
     *  k: size of extended states
     *  f_v: position of the already assigned variable
     *  variables: all variables in state
     *  mutexes: all mutex tuples
     */
    static int c_k_f(vector<State> states, int k, int f_v, VariablesProxy &variables, vector<Tuple> &mutexes){
        int sum = 0;
        int mult;
        //TODO: state mitgeben und id der gesetzten variable? so kann das array für den state (arr = -1 für jeden state) einmal initialisiert, und dann im vector (vec == arr initialisieren) geändert werden.
        // get all extensions kann auch genutzt werden in opt mk, mit leerem state und k = 1... und vollständiger domain...dort ebenfalls unreachables weglassen?
        for (State e : states){
            mult = 1;
            vector<vector<int>> domains = multi_fact_diasambiguation(e, variables, mutexes); // This probably does some redundant work, would it be useful to give more narrowed down domains?
            for (vector<int> d : domains){
                mult *= d.size();
            }
            sum += mult;
        }
        return sum;
    }

    // weights of Eq 12, inner vector for facts, outer for the different variables. ordered the same way as
    // variables in proxy, i suppose....
    static vector<vector<float>> opt_k_m(int k, VariablesProxy &variables, vector<Tuple> &mutexes, AbstractTask &task){
        // variables to return
        vector<State> extended_states;
        vector<vector<float>> weights;

        //needed to calculate the weights
        vector<float> weights_f;
        float w;
        float sum;

        //needed to get all states
        vector<int> values(variables.size(), -1);
        State state = State(task, move(values));
        vector<State> states;
        vector<vector<int>> domains;
        //TODO: generate all states and look for them, instead of generating all multiple times?

        for (int i = 0; i < variables.size(); i++) {
            sum = 0;
            weights_f.clear();
            for (int d = 0; d < variables[i].get_domain_size(); d++) {
                //get and add the new state and all its extensions
                state.add_variable(i, d);
                domains = multi_fact_diasambiguation(state, variables, mutexes);
                states = get_all_extensions(state, k, i, domains);
                extended_states.insert(extended_states.end(), states.begin(), states.end());
                
                //get and add the weight of this fact
                w = (float) c_k_f(states, k, i, variables, mutexes);
                weights_f.push_back(w);
                sum += w;
                
                states.clear();
            }
            state.add_variable(i, -1);
            for (int d = 0; d < variables[i].get_domain_size(); d++) {
                weights_f[d] /= sum;
            }
            weights.push_back(weights_f);
        }
        //TODO: also return the states
        return weights;
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

        optimizer.optimize_for_weighted_samples

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