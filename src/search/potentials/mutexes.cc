//
// Created by salome on 26/09/2020.
//

#include "mutexes.h"

MutexTable::MutexTable(const Options &opts, VariablesProxy variables, State &state)
        : variables(variables) {
    auto hm = make_shared<HMHeuristic>(opts);
    mutexes = hm->get_unreachable_tuples(state);
    hm.reset();
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
    for (vector<FactPair> mutex : mutexes) {
        if (!unassigned(state, mutex[0].var) && state.find(mutex[0].var)->second == mutex[0].value) {
            mf.push_back(mutex[1]);
        } else if (!unassigned(state, mutex[1].var) && state.find(mutex[1].var)->second == mutex[1].value) {
            mf.push_back(mutex[0]);
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
    for (vector<FactPair> mutex : mutexes) {
        if (mutex[0].var == variable && mutex[0].value == value) {
            mf.push_back(mutex[1]);
        } else if (mutex[1].var == variable && mutex[1].value == value) {
            mf.push_back(mutex[0]);
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
