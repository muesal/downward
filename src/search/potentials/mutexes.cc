#include "mutexes.h"

MutexTable::MutexTable(TaskProxy task_proxy)
        : variables(task_proxy.get_variables()),
          task_proxy(task_proxy) {
    utils::g_log << "Start building mutex table." << endl;
    vector<FactPair> s_tup = task_properties::get_fact_pairs(task_proxy.get_initial_state());
    init_hm_table(s_tup);
    update_hm_table();

    for (pair<Pair, int> pair: this->hm_table) {
        if (pair.second == 1) {
            mutexes.push_back(pair.first);
        }
    }
    utils::g_log << "Built mutex table." << endl;
}

bool MutexTable::unassigned(map<int, int> &state, int variable_id) {
    return state.find(variable_id) == state.end();
}

/**
 * Used in fact_disambiguation's, to subtract all non-possible values (mutexes) from the domain of a variable
 * @param domain domain of the variable
 * @param mutex facts which are mutex with the current state
 * @param variable_id id of the variable
 * @return the remaining domain of the variable
 */
vector<int> MutexTable::set_minus(const vector<int> &domain, const vector<FactPair> &mutex, int variable_id) {
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
 */
vector<FactPair> MutexTable::intersection(const vector<FactPair> &one, const vector<FactPair> &two) {
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
  * Save all facts which are unreachable for this (partial) state
  * @param variable id of the fact
  * @param value of the fact
  * @param mutexes set of mutexes
  * @return mf tuple in which to save the mutexes
  */
vector<FactPair> MutexTable::get_mutex_with_state(map<int, int> &state) {
    vector<FactPair> mf;
    for (Pair mutex : mutexes) {
        if (!unassigned(state, mutex.first.var) && state.find(mutex.first.var)->second == mutex.first.value) {
            mf.push_back(mutex.second);
        } else if (!unassigned(state, mutex.second.var) && state.find(mutex.second.var)->second == mutex.second.value) {
            mf.push_back(mutex.first);
        }
    }
    return mf;
}

/**
   * save all FactPairs which are mutex to the fact <variable, value>
   * @param variable id of the fact
   * @param value of the fact
   * @param mutexes set of mutexes
   * @param mf tuple in which to save the mutexes
   */
void MutexTable::get_mutex_with_fact(int variable, int value, vector<FactPair> &mf) {
    for (Pair mutex : mutexes) {
        if (mutex.first.var == variable && mutex.first.value == value) {
            mf.push_back(mutex.second);
        } else if (mutex.second.var == variable && mutex.second.value == value) {
            mf.push_back(mutex.first);
        }
    }
}

/**
* Smaller all the domains of all variables by pruning unreachable facts
* @param state the state of interest as vector of facts (pairs<int, int>)
* @param variables all variables of this task
* @param mutexes the set of mutexes for this task
* @return the smalled domains
*/
vector<vector<int>>
MutexTable::multi_fact_disambiguation(map<int, int> &state) {
    // get all domains and id's of the variables (D_V<- F_V for every V in Variables)
    vector<vector<int>> domains;
    domains.reserve(variables.size());
    for (VariableProxy v : variables) {
        vector<int> dom;
        // For already assigned values add only this fact, for the rest the whole domain
        if (!unassigned(state, v.get_id())) {
            dom.push_back(state[v.get_id()]);
        } else {
            for (int d = 0; d < v.get_domain_size(); d++) {
                dom.push_back(d);
            }
        }
        domains.push_back(dom);
    }

    // A <- M_p
    vector<FactPair> a = get_mutex_with_state(state);
    a.reserve(mutexes.size());

    //while state changes check for single fact disambiguations and assign them
    bool changed = true;
    while (changed) {
        changed = false;
        for (VariableProxy v : variables) {
            int id = v.get_id();
            if (unassigned(state, id)) {    // Variable is not assigned
                auto size = domains[id].size();
                //algorithm 2 l. 7 (D_V <- D_V \ a)
                domains[id] = set_minus(domains[id], a, id);

                //algorithm 2 l. 8 (A <- A U [intersection over all facts in D_V of M_(p U {f})] )
                if (domains[id].size() < size) {
                    if (domains[id].empty()) {  // The state is a dead end.
                        for (auto &domain : domains) {
                            domain.clear();
                        }
                        return domains; // return empty domains, as none are leading
                    }
                    vector<FactPair> mf, mf2;
                    get_mutex_with_fact(id, domains[id][0], mf);
                    for (size_t f = 1; f < domains[id].size(); f++) {
                        get_mutex_with_fact(id, domains[id][f], mf2);
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

const VariablesProxy *MutexTable::getVariablesProxy() const {
    auto var = &variables;
    return var;
}

/**
 * Generate the hm_table.
 * Expand and store all pairs of facts with the unreachable value (1).
 */
void MutexTable::generate_all_pairs() {
    int num_variables = task_proxy.get_variables().size();
    for (int i = 0; i < num_variables; ++i) {
        for (int j = 0; j < task_proxy.get_variables()[i].get_domain_size(); ++j) {
            Pair pair({i, j}, {i, j}); // single facts must be in the table as well
            hm_table[pair] = 0;
            for (int i2 = i + 1; i2 < num_variables; ++i2) {
                for (int j2 = 0; j2 < task_proxy.get_variables()[i2].get_domain_size(); ++j2) {
                    Pair pair2({i, j}, {i2, j2});
                    hm_table[pair2] = 1;
                }
            }
        }
    }
}

/**
 * Generate all pairs from the given facts.
 * @param t list f facts
 * @param base generated pairs
 */
vector<MutexTable::Pair> MutexTable::generate_all_pairs(vector<FactPair> &t) {
    vector<Pair> base;
    for (size_t i = 0; i < t.size(); ++i) {
        base.emplace_back(t[i], t[i]); // for reachability of single facts
        for (size_t j = i + 1; j < t.size(); ++j) {
            base.emplace_back(t[i], t[j]);
        }
    }
    return base;
}

/**
 * All pairs which are in t are reachable. Set their values in the hm-table to zero.
 * @param t initial state
 */
void MutexTable::init_hm_table(vector<FactPair> &t) {
    generate_all_pairs();
    vector<Pair> pairs = generate_all_pairs(t);
    for (Pair pair : pairs) {
        hm_table[pair] = 0;
    }
}

/**
 * While something changes, go through all operators. For the operators where all preconditions are reachable, set all
 * effects to reachable, and delete the operator from the list, as there is no need to look at it again. Stop, when
 * nothing changes anymore.
 */
void MutexTable::update_hm_table() {
    was_updated = true;
    while (was_updated) {
        was_updated = false;
        for (OperatorProxy op : task_proxy.get_operators()) {
            vector<FactPair> pre = task_properties::get_fact_pairs(op.get_preconditions());
            sort(pre.begin(), pre.end());

            bool all = all_reachable(pre);
            if (all) { // if all preconditions are reachable
                vector<FactPair> eff;
                set<int> op_vars;
                for (EffectProxy effect : op.get_effects()) {
                    eff.push_back(effect.get_fact().get_pair());
                    op_vars.insert(effect.get_fact().get_pair().var);
                }
                sort(eff.begin(), eff.end());

                vector<Pair> effs = generate_all_pairs(eff);
                for (Pair pair : effs) {
                    if (hm_table[pair] == 1) {
                        hm_table[pair] = 0;
                        was_updated = true;
                    }
                }

                for (FactPair f : eff) {
                    extend_fact(f, pre, op_vars);
                }
            }
        }
    }
}

/**
 * Go through all hm_entries for the fact pairs of t, if any of them is not reachable return false.
 * @param t vector of facts
 * @return true if all facts in t are reachable
 */
bool MutexTable::all_reachable(vector<FactPair> &t) {
    vector<Pair> pairs = generate_all_pairs(t);
    for (Pair pair : pairs) {
        if (hm_table[pair] == 1) {
            return false;
        }
    }
    return true;
}

/**
 * Retrn True whenthe fact is reachable with all preconditions.
 * Same as all_reachable, with the assumption that preconditions are reachable themselves.
 * @param preconditions
 * @param fact
 * @return true when all reachable
 */
bool MutexTable::reachable_from(const vector<FactPair> & preconditions, FactPair fact) {
    int var = fact.var;
    for (FactPair pre : preconditions) {
        if ( var == pre.var && fact.value != pre.value) {
            return false; // contradicts the precondition
        } else if (var < pre.var ) {
            Pair pair(fact, pre);
            if (hm_table[pair] == 1) { // not reachable
                return false;
            }
        } else {
            Pair pair(pre, fact);
            if (hm_table[pair] == 1) { // not reachable
                return false;
            }
        }
    }
    return true;
}

/**
 * Go over all facts f2 of variables which are neither the fact, in the precondition nor the effect (i.e., not in op_vars)
 * and add <f,f2> as reachable, if precondition u {f2} are reachable.
 * @param fact the fact
 * @param pre all FactPairs of the precondition
 * @param eff_vars list of all variables in the effeect and the preconditions of the operator.
 */
void MutexTable::extend_fact(const FactPair &fact, const vector<FactPair> &pre, const set<int> &eff_vars) {
    for (int i = 0; i < fact.var; ++i) { // all variables
        if (eff_vars.count(i) == 0) { // which are on in the effect
            for (int j = 0; j < task_proxy.get_variables()[i].get_domain_size(); ++j) {
                // if {preconditions u f2} is reachable set {f,f2} to reachable
                FactPair f2(i, j);
                Pair pair(f2, fact);
                if (hm_table[pair] == 1) {
                    if (reachable_from(pre, f2)) {
                        hm_table[pair] = 0;
                        was_updated = true;
                    }
                }
            }
        }
    }
    for (int i = fact.var + 1; i < variables.size(); ++i) {
        if (eff_vars.count(i) == 0) { // which are on in the effect
            for (int j = 0; j < task_proxy.get_variables()[i].get_domain_size(); ++j) {
                // if {preconditions u f2} is reachable set {f,f2} to reachable
                FactPair f2(i, j);
                Pair pair(fact, f2);
                if (hm_table[pair] == 1) {
                    if (reachable_from(pre, f2)) {
                        hm_table[pair] = 0;
                        was_updated = true;
                    }
                }
            }
        }
    }
}
