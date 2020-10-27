#include "potential_function.h"
#include "potential_heuristic.h"
#include "potential_max_heuristic.h"
#include "potential_optimizer.h"
#include "mutexes.h"
#include "util.h"

#include <memory>
#include <vector>
#include <set>
#include <map>

using namespace std;
using Tuple = vector<FactPair>;
using Weight = tuple<int, int, long double>; // Variable id, value and corresponding weight

namespace potentials {

    /**
     * @param state the state
     * @param variable_id the id of the variable
     * @return whether this variables is assigned
     */
    static bool unassigned(map<int, int> &state, int variable_id) {
        return state.find(variable_id) == state.end();
    }

    /**
      * Extend the partial state to all possible states with k more assigned values
      * @param state vector containing the facts of a partial state
      * @param k amount of remaining variables to assign
      * @param assigned vector of variable ids which were assigned from the beginning
      * @param domains domains of all variables
      * @param extended the states which have been extended so far
      * @return a vector of all extended vectors
      */
    static vector<map<int, int>>
    get_all_extensions(map<int, int> &state, int k, int last_assigned, vector<int> &assigned,
                       vector<vector<int>> &domains, vector<map<int, int>> &extended) {
        if (k == 0) {
            extended.push_back(state);  // Add the state to the list, if it is fully extended.
            return extended;
        }

        // variables assigned from beginning to end, therefore enough unassigned variables need to remain free
        int before = assigned.size();
        while (before > 0 && last_assigned < assigned[before - 1]) {
            before--;
        }
        before = (int) assigned.size() - before - (k - 1);
        for (size_t i = last_assigned; i < domains.size() -
                                           before; i++) { //domains.size(), as there are exactly as many domains as variables in the state
            if (unassigned(state, i)) { // if the variable is unassigned
                for (int d : domains[i]) {
                    state[i] = d;               // assign it to all domains and get all extensions
                    get_all_extensions(state, k - 1, (int) i + 1, assigned, domains, extended);
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
        if (k <= 0) {
            states.push_back(state);
            return states;
        }
        states.reserve(min(domains.size() * domains[0].size() * k, states.max_size()));
        vector<int> variables;
        variables.reserve(state.size());
        for (auto s : state) {
            variables.push_back(s.first);
        }
        return get_all_extensions(state, k, 0, variables, domains, states);
    }

    /**
     * Calculate C^k_f(M) for one state, as it handles a list of states it can also be used for K_f
     * it is the sum [ of the products [ of the number of reachable facts ] over all variables ] over all (k-extended) states
     * @param states the extensions of the state of interest (P^{f}_k resp. P^(tU{f})_(|t|+k))
     * @param variables the variables of this task
     * @param mutexes the set of mutexes for this task
     * @return the sum over the product of all reachable facts of each state
     */
    static long double c_k_f(vector<map<int, int>> &states, MutexTable &table) {
        long double sum = 0;
        long double mult;
        for (map<int, int> &e : states) {
            mult = 1;
            vector<vector<int>> domains = table.multi_fact_disambiguation(e);
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
      * @param k amount of variables to extend
      * @param variables all variables
      * @param mutexes all mutex facts
      * @param task the task
      * @return a vector containing for each state (outer vector) the weights of each fact (inner vectors)
      */
    static vector<Weight>
    opt_k_m(int k, map<int, int> &assigned_variables, MutexTable &table) {
        const VariablesProxy *variables = table.getVariablesProxy();
        vector<Weight> facts;               // vector containing facts and their corresponding weight.
        facts.reserve(min(variables->size() * variables->operator[](0).get_domain_size(), facts.max_size()));
        vector<Weight> weights_f;           // to temporarily store the c_k_f values of one variable
        vector<map<int, int >> states;      // to temporarily store the extended states of a fact
        vector<vector<int>> domains;        // temporarily contain the multi_fact_disambiguated domain with one fact
        long double w, sum;

        for (size_t i = 0; i < variables->size(); i++) {         // for all all unassigned variables V and each fact_V
            if (unassigned(assigned_variables, (int) i)) {
                sum = 0;
                weights_f.clear();
                for (int d = 0; d < variables->operator[](i).get_domain_size(); d++) {
                    assigned_variables[i] = d;                          // add this fact
                    domains = table.multi_fact_disambiguation(assigned_variables);    // get all non-mutex domains
                    if (domains[0].empty()) continue;                   // partial state is a dead end
                    states = get_all_extensions(assigned_variables, k - 1,
                                                domains);           // get all extended states for this fact
                    // (k - 1, as one additional variables is already assigned)
                    w = c_k_f(states, table);                    // Get the non-normalized weight of this fact ...
                    states.clear();
                    weights_f.emplace_back(i, d, w);                    // and add it to the list.
                    sum += w;
                }

                // normalize all weights
                for (auto & d : weights_f) {
                    get<2>(d) = get<2>(d) / sum;
                    facts.push_back(d);
                }

                // remove state[i]
                assigned_variables.erase(assigned_variables.find(i));
            }
        }
        return facts;
    }

    /**
      * Calculates the weights of all facts for all states (Eq. 12)
      * @param k amount of variables to extend
      * @param variables all variables
      * @param mutexes all mutex facts
      * @param task the task
      * @return a vector containing for each state (outer vector) the weights of each fact (inner vectors)
      */
    static vector<Weight> opt_k_m(int k, MutexTable &table) {
        utils::g_log << "Using Mutex Based Potential Heuristics." << endl;
        map<int, int> assigned_variables;   // to temporarily store one 'state'
        return opt_k_m(k, assigned_variables, table);
    }

    /**
     * Calculate n optimization functions, using opt_t_k_m (eq. 14)
     * @param t size of partial state (t uniformly randomly chosen facts)
     * @param k
     * @param n
     * @param variables
     * @param mutexes
     * @param optimizer
     * @return
     */
    static vector<unique_ptr<PotentialFunction>>
    opt_t_k_m(int t, int k, int n, MutexTable table, PotentialOptimizer &optimizer) {
        utils::g_log << "Using Mutex Based Potential Ensemble Heuristics." << endl;
        const VariablesProxy *variables = table.getVariablesProxy();
        random_device rd;
        mt19937 gen(rd()); // apparently this is necessary for uniform int distribution...
        uniform_int_distribution<> dist_v(0, variables->size() - 1);

        // generate n optimization functions for one randomly sampled state
        vector<unique_ptr<PotentialFunction>> functions;
        functions.resize(n);
        vector<Weight> weights;
        map<int, int> state;
        for (int i = 0; i < n; i++) {
            state.clear();
            for (int t_i = 0; t_i < t; t_i++) {
                int v = dist_v(gen);
                if (unassigned(state, v)) {
                    uniform_int_distribution<> dist_f(0, variables->operator[](v).get_domain_size() - 1);
                    state[v] = dist_f(gen);
                } else t_i--;
            }
            weights = opt_k_m(k - t, state, table);
            if (weights.empty()) {
                n--; // State was a dead end;
            } else {
                optimizer.optimize_for_weighted_samples(weights);
                functions[i] = optimizer.get_potential_function();
                utils::g_log << "Calculated " << i << "potential functions." << endl;
            }
        }
        return functions;
    }

    static unique_ptr<PotentialFunction> create_mutex_based_potential_function(
            Options &opts) {
        PotentialOptimizer optimizer(opts);
        const AbstractTask &task = *opts.get<shared_ptr<AbstractTask>>("transform");
        TaskProxy task_proxy(task);
        VariablesProxy variables = task_proxy.get_variables();

        MutexTable *table = optimizer.get_mutex_table();
        if (table == nullptr) {
            table = new MutexTable(task_proxy);
        }

        int k = opts.get<int>("k");
        assert(k < (int) variables.size()); // k may not be bigger than the size of one state
        vector<Weight> weights = opt_k_m(k, *table);
        optimizer.optimize_for_weighted_samples(weights);
        return optimizer.get_potential_function();
    }

    static vector<unique_ptr<PotentialFunction>> create_mutex_based_ensemble_potential_function(
            const Options &opts) {
        PotentialOptimizer optimizer(opts);
        const AbstractTask &task = *opts.get<shared_ptr<AbstractTask>>("transform");
        TaskProxy task_proxy(task);
        VariablesProxy variables = task_proxy.get_variables();

        MutexTable *table = optimizer.get_mutex_table();
        if (table == nullptr) {
            table = new MutexTable(task_proxy);
        }

        int k = opts.get<int>("k");
        assert(k < (int) variables.size()); // k may not be bigger than the size of one state
        int t = opts.get<int>("t");
        assert(t <= k);                     // t may not be bigger than k
        int n = opts.get<int>("n");

        return opt_t_k_m(t, k, n, *table, optimizer);
    }

    static shared_ptr<Heuristic> _parse_single(OptionParser &parser) {
        parser.document_synopsis(
                "Mutex-based potential heuristics",
                get_admissible_potentials_reference());
        parser.add_option<int>(
                "mutex",
                "Use mutexes in potential optimizer",
                "1",
                Bounds("0", "1"));
        parser.add_option<int>(
                "init-const",
                "Use the addition constraint on the initial state in potential optimizer",
                "0",
                Bounds("0", "1"));
        parser.add_option<int>(
                "m",
                "use h2 heuristic",
                "2",
                Bounds("0", "infinity"));
        parser.add_option<int>(
                "k",
                "size of extended state",
                "1",
                Bounds("0", "infinity"));
        prepare_parser_for_admissible_potentials(parser);
        Options opts = parser.parse();
        if (parser.dry_run())
            return nullptr;

        return make_shared<PotentialHeuristic>(
                opts, create_mutex_based_potential_function(opts));
    }

    static shared_ptr<Heuristic> _parse_ensemble(OptionParser &parser) {
        parser.document_synopsis(
                "Mutex-based ensemble potential heuristics",
                get_admissible_potentials_reference());
        parser.add_option<int>(
                "mutex",
                "Use mutexes in potential optimizer",
                "1",
                Bounds("0", "1"));
        parser.add_option<int>(
                "init-const",
                "Use the addition constraint on the initial state in potential optimizer",
                "0",
                Bounds("0", "1"));
        parser.add_option<int>(
                "m",
                "use h2 heuristic",
                "2",
                Bounds("0", "infinity"));
        parser.add_option<int>(
                "k",
                "size of extended state",
                "2",
                Bounds("0", "infinity"));
        parser.add_option<int>(
                "t",
                "amount of uniformly randomly chosen facts (t <= k)",
                "1",
                Bounds("0", "infinity"));
        parser.add_option<int>(
                "n",
                "Number of states to sample",
                "50",
                Bounds("0", "infinity"));
        prepare_parser_for_admissible_potentials(parser);
        Options opts = parser.parse();
        utils::add_rng_options(parser);
        if (parser.dry_run())
            return nullptr;

        return make_shared<PotentialMaxHeuristic>(
                opts, create_mutex_based_ensemble_potential_function(opts));
    }

    static Plugin<Evaluator> _plugin_single(
            "mutex_based_potential", _parse_single, "heuristics_potentials");

    static Plugin<Evaluator> _plugin_ensemble(
            "mutex_based_ensemble_potential", _parse_ensemble, "heuristics_potentials");
}