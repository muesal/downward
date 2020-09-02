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
#include <set>
#include <map>

//For the table:
#include <iostream>
#include <fstream>

using namespace std;
using namespace hm_heuristic;
using Tuple = vector<FactPair>;
using Weight = tuple<int, int, long double>; // Variable id, value and corresponding weight

namespace potentials {

    /**
     * Used in fact_disambiguation's, to subtract all non-possible values (mutexes) from the domain of a variable
     * @param domain domain of the variable
     * @param mutex facts which are mutex with the current state
     * @param variable_id id of the variable
     * @return the remaining domain of the variable
     */ //TODO: better name...
    vector<int> set_minus(const vector<int> &domain, const Tuple &mutex, int variable_id) {
        set<int> dom(domain.begin(), domain.end()); // Sets are unique
        for (FactPair m : mutex) {                  // For each fact in mutex...
            if (m.var == variable_id) {             // which affects this variable...
                dom.erase(m.value);                 // remove its value from the domain.
            }
        }
        return vector<int>(dom.begin(), dom.end());
    }

    /**
     * Used in multi_fact_disambiguation to get elements which are in both mutex sets.
     * @param one mutex-set one
     * @param two mutex-set two
     * @return three facts i both one and two
     */ //TODO: better name...
    static Tuple intersection(const Tuple &one, const Tuple &two) {
        Tuple three;
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
     * save all FactPairs which are mutex to the fact <variable, value>
     * @param variable id of the fact
     * @param value of the fact
     * @param mutexes set of mutexes
     * @param mf tuple in which to save the mutexes
     */
    static void get_mutex_with_fact(int variable, int value, vector<Tuple> &mutexes, Tuple &mf) {
        for (Tuple mutex : mutexes) {
            if (mutex[0].var == variable && mutex[0].value == value) {
                mf.push_back(mutex[1]);
            } else if (mutex[1].var == variable && mutex[1].value == value) {
                mf.push_back(mutex[0]);
            }
        }
    }

    /**
      * Save all mutex pairs which are unreachable for this (partial) state
      * @param variable id of the fact
      * @param value of the fact
      * @param mutexes set of mutexes
      * @return mf tuple in which to save the mutexes
      */ //TODO: löschen, falls multi_fact_disambiguation (state) gelöscht wurde
    static Tuple get_mutex_with_state(State &state, vector<Tuple> &mutexes) {
        Tuple mf;
        for (Tuple mutex : mutexes) {
            if (state[mutex[0].var].get_value() == mutex[0].value) {
                mf.push_back(mutex[1]);
            } else if (state[mutex[1].var].get_value() == mutex[1].value) {
                mf.push_back(mutex[0]);
            }
        }
        return mf;
    }

    /**
      * Save all facts which are unreachable for this (partial) state
      * @param variable id of the fact
      * @param value of the fact
      * @param mutexes set of mutexes
      * @return mf tuple in which to save the mutexes
      */
    static Tuple get_mutex_with_state(map<int, int> &state, vector<Tuple> &mutexes) {
        Tuple mf;
        for (Tuple mutex : mutexes) {
            if (state[mutex[0].var] == mutex[0].value) {
                mf.push_back(mutex[1]);
            } else if (state[mutex[1].var] == mutex[1].value) {
                mf.push_back(mutex[0]);
            }
        }
        return mf;
    }

    /**
      * Extend the partial state to all possible states with k more assigned values
      * @param state vector containing the facts of a partial state
      * @param k amount of remaining variables to assign
      * @param domains domains of all variables
      * @param extended the states which have been extended so far
      * @return a vector of all extended vectors
      */
    static vector<map<int, int>>
    get_all_extensions(map<int, int> &state, int k, int last_assigned, vector<vector<int>> &domains,
                       vector<map<int, int>> &extended) {
        if (k == 0) {
            // TODO: does this realy save a copy of the state?
            extended.push_back(state);  // Add the state to the list, if it is fully extended.
            return extended;
        }

        // variables assigned from beginning to end, therefore enough unassigned variables need to remain free
        // size - k if only not assigned values left, size - k - 1 otherwise
        int before = state.begin()->second < last_assigned ? 0 : 1;
        for (size_t i = last_assigned; i < domains.size() - k - before; i++) { //domains.size(), as there are exactly as many domains as variables in the state
            if (state.find(i) == state.end()) { // if the variable is unassigned
                for (int d : domains[i]) {
                    state[i] = d;               // assign it to all domains and get all extensions
                    get_all_extensions(state, k - 1, i + 1, domains, extended); //TODO: does this change extended?
                }
                state.erase(i);                 // unassign the variable
            }
        }
        return extended;
    }

    /**
      * Extend the partial state to all possible states with k more assigned values
      * @param state vector containing the facts of a partial state
      * @param task
      * @param k amount of remaining variables to assign
      * @param domains domains of all variables
      * @return a vector of all extended states as vectors of facts
      */
    static vector<map<int, int>>
    get_all_extensions(map<int, int> &state, int k, vector<vector<int>> &domains) {
        vector<map<int, int>> states;
        return get_all_extensions(state, k, 0, domains, states);
    }

    /**
     * TODO: The following probably have the first 11 lines in common... AND is single_fact needed? can it be used somewhere to simplify things, eg. et the beginning of mutli fact?
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
        Tuple mp = get_mutex_with_state(state, mutexes);

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
    */

    /**
     * Smaller all the domains of all variables by pruning unreachable facts
     * @param state the state of interest
     * @param variables all variables of this task
     * @param mutexes the set of mutexes for this task
     * @return the smalled domains
     */
    //TODO: löschen?
    static vector<vector<int>>
    multi_fact_diasambiguation(State state, VariablesProxy &variables, vector<Tuple> &mutexes) {
        //get all domains and id's of the variables
        vector<vector<int>> domains;
        vector<int> variable_id;
        for (VariableProxy v : variables) {
            vector<int> dom;
            // TODO: for assigned values, add all values or only the assigned one?
            if (!state.is_defined(v.get_id())) {
                for (int d = 0; d < v.get_domain_size(); d++) {
                    dom.push_back(d);
                }
            } else {
                dom.push_back(state[v.get_id()].get_value());
            }
            domains.push_back(dom);
            variable_id.push_back(v.get_id());
        }
        Tuple a = get_mutex_with_state(state, mutexes);

        //while state changes check for single fact disambiguations and assign them
        bool changed = true;
        while (changed) {
            changed = false;
            for (size_t v = 0; v < state.size(); v++) {
                if (!state.is_defined(v)) {
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
                            if (mf2.empty() || mf.empty()) {
                                break;
                            }
                        }
                        if (!mf.empty()) {
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
     * Smaller all the domains of all variables by pruning unreachable facts
     * @param state the state of interest as vector of facts (pairs<int, int>)
     * @param variables all variables of this task  TODO: it could also just get all domains as vector<vector<int>>?
     * @param mutexes the set of mutexes for this task
     * @return the smalled domains
     */
    static vector<vector<int>>
    multi_fact_disambiguation(map<int, int> state, VariablesProxy &variables, vector<Tuple> &mutexes) {
        // get all domains and id's of the variables (D_V<- F_V for every V in Variables)
        vector<vector<int>> domains;
        for (VariableProxy v : variables) {
            vector<int> dom;
            // For already assigned values add only this fact, for the rest the whole domain TODO: is this correct?
            auto i = state.find(v.get_id());
            if (i != state.end()) {
                dom.push_back(i->second);
            } else {
                for (int d = 0; d < v.get_domain_size(); d++) {
                    dom.push_back(d);
                }
            }
            domains.push_back(dom);
        }
        // A <- M_p
        Tuple a = get_mutex_with_state(state, mutexes);

        //while state changes check for single fact disambiguations and assign them
        bool changed = true;
        while (changed) {
            changed = false;
            for (VariableProxy v : variables) {
                if (state.find(v.get_id()) == state.end()) {    // Variable is not assigned
                    auto size = domains[v.get_id()].size();
                    //algorithm 2 l. 7 (D_V <- D_V \ a)
                    domains[v.get_id()] = set_minus(domains[v.get_id()], a, v.get_id());

                    //algorithm 2 l. 8 (A <- A U [intersection over all facts in D_V of M_(p U {f})] )
                    if (domains[v.get_id()].size() < size) {
                        vector<FactPair> mf, mf2;
                        get_mutex_with_fact(v.get_id(), domains[v.get_id()][0], mutexes, mf);
                        for (size_t f = 1; f < domains[v.get_id()].size(); f++) {
                            get_mutex_with_fact(v.get_id(), domains[v.get_id()][f], mutexes, mf2);
                            mf = intersection(mf, mf2); // Get all facts which are not part of any reachable state
                            if (mf.empty()) break;
                            mf2.clear();
                        }
                        if (!mf.empty()) {
                            a.insert(a.begin(), mf.begin(), mf.end());  // add states to a
                        }
                        changed = true;
                    }
                }
            }
        }
        return domains;
    }


    /**
     * Calculate C^k_f(M) for one state, as it handles a list of states it can also be used for K_f
     * it is the sum [ of the products [ of the number of reachable facts ] over all variables ] over all (k-extended) states
     * @param states the extensions of the state of interest (P^{f}_k resp. P^(tU{f})_(|t|+k))
     * @param variables the variables of this task
     * @param mutexes the set of mutexes for this task
     * @return the sum over the product of all reachable facts of each state
     */
    static long double c_k_f(const vector<map<int, int>> &states, VariablesProxy &variables, vector<Tuple>
    &mutexes) {
        long double sum = 0;
        long double mult;
        for (const map<int, int> &e : states) {
            mult = 1;
            vector<vector<int>> domains = multi_fact_disambiguation(e, variables,
                                                                    mutexes); // This probably does some redundant work, would it be useful to give more narrowed down domains?
            for (const vector<int> &d : domains) {
                mult *= d.size();
                if (mult == 0) break;
            }
            sum += mult;
        }
        return sum;
    }

/**
  * Calculates the weights of all facts for all states (Eq. 12)
  * @param k size extended states should have
  * @param variables all variables
  * @param mutexes all mutex facts
  * @param task the task
  * @return a vector containing for each state (outer vector) the weights of each fact (inner vectors)
  */
    static vector<Weight> opt_k_m(int k, VariablesProxy &variables, vector<Tuple> &mutexes) {
        // TODO: weights_f could be replaced by a pointer to the last entry of the last variable in facts...
        vector<Weight> facts;               // vector containing facts and their corresponding weight.
        vector<Weight> weights_f;           // to temporarily store the c_k_f values of one variable
        map<int, int> assigned_variables;   // to temporarily store one 'state'
        vector<map<int, int >> states;      // to temporarily store the extended states of a fact
        vector<vector<int>> domains;        // temporarily contain the multi_fact_disambiguated domain with one fact
        long double w, sum;

        for (size_t i = 0; i < variables.size(); i++) {         // for all variables and each fact_V
            sum = 0;
            weights_f.clear();
            for (int d = 0; d < variables[i].get_domain_size(); d++) {
                assigned_variables.clear();                     // empty the 'state'
                assigned_variables[i] = d;                      // add this fact
                domains = multi_fact_disambiguation(assigned_variables, variables,
                                                    mutexes);// get all non-mutex domains
                states = get_all_extensions(assigned_variables, k,
                                            domains);           // get all extended states for this fact
                w =  c_k_f(states, variables,
                                  mutexes);                  // Get the non-normalized weight of this fact ... TODO: hand over domains instead of variables?))
                states.clear();
                weights_f.emplace_back(i, d, w);                // and add it to the list.
                sum += w;
            }

            // normalize all weights
            for (int d = 0; d < variables[i].get_domain_size(); d++) {
                get<2>(weights_f[d]) = get<2>(weights_f[d]) / sum;
                facts.push_back(weights_f[d]);
            }
        }

        return facts;
    }

    static unique_ptr<PotentialFunction> create_mutex_based_potential_function(
            const Options &opts) {

        PotentialOptimizer optimizer(opts);
        const AbstractTask &task = *opts.get<shared_ptr<AbstractTask>>("transform");
        TaskProxy task_proxy(task);
        VariablesProxy variables = task_proxy.get_variables();

        State initial = task_proxy.get_initial_state();

        //get all mutexes
        //This creates a file containing the mutex table. If one already exists, it will be used instead.
        std::vector<Tuple> mutexes;
        string filename = "mutex_table.txt";
        ifstream table;
        table.open(filename);
        if (!table) {
            cout << "Perform hm-heuristic and store mutexes to " << filename << endl;
            auto hm = make_shared<HMHeuristic>(opts);
            mutexes = hm->get_unreachable_tuples(initial);
            hm.reset();
            //write to file
            ofstream create(filename);
            if (create.is_open()) {
                for (Tuple mutex : mutexes) {
                    create << mutex[0].var << " " << mutex[0].value << " ";
                    create << mutex[1].var << " " << mutex[1].value << "\n";
                }
                create.close();
            } else {
                cout << "Unable to store mutex table";
            }
        } else {
            //read from file
            mutexes.clear();
            string line;
            string number;
            int var, value;
            while (getline(table, line)) {
                istringstream ss(line);
                ss >> number;
                var = stoi(number);
                ss >> number;
                value = stoi(number);
                FactPair fact1(var, value);
                ss >> number;
                var = stoi(number);
                ss >> number;
                value = stoi(number);
                FactPair fact2(var, value);

                mutexes.push_back({fact1, fact2});
            }
            table.close();
        }

        int k = 1; //TODO: get from options
        vector<Weight> weights = opt_k_m(k, variables, mutexes);
        optimizer.optimize_for_weighted_samples(weights);
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