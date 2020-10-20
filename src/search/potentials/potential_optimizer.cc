#include "potential_optimizer.h"

#include "potential_function.h"

#include "../option_parser.h"

#include "../task_utils/task_properties.h"
#include "../utils/collections.h"
#include "../utils/memory.h"
#include "../utils/system.h"

#include <limits>
#include <unordered_map>

using namespace std;
using utils::ExitCode;

namespace potentials {
static int get_undefined_value(VariableProxy var) {
    return var.get_domain_size();
}

static int get_undefined_value_for_operator(VariableProxy var, OperatorProxy o) {
    return var.get_domain_size() + o.get_id() + 1;
}

PotentialOptimizer::PotentialOptimizer(const Options &opts)
    : task(opts.get<shared_ptr<AbstractTask>>("transform")),
      task_proxy(*task),
      lp_solver(opts.get<lp::LPSolverType>("lpsolver")),
      max_potential(opts.get<double>("max_potential")),
      num_lp_vars(0),
      use_mutexes(opts.get<int>("mutex") == 1),
      initial_constraint(opts.get<int>("init-const") == 1) {
    task_properties::verify_no_axioms(task_proxy);
    task_properties::verify_no_conditional_effects(task_proxy);
    if (use_mutexes) {
        State initial = task_proxy.get_initial_state();
        VariablesProxy vars = task_proxy.get_variables();
        table = new MutexTable(opts, vars, initial);
    }
    initialize();
}

void PotentialOptimizer::initialize() {
    VariablesProxy vars = task_proxy.get_variables();
    std::size_t num_operators = use_mutexes ? task_proxy.get_operators().size() : 0;
    lp_var_ids.resize(vars.size() * ( 1 + num_operators));
    fact_potentials.resize(vars.size());
    for (VariableProxy var : vars) {
        lp_var_ids[var.get_id()].resize(var.get_domain_size() + 1 + num_operators);
        for (int val = 0; val < var.get_domain_size() + 1; ++val) {
            lp_var_ids[var.get_id()][val] = num_lp_vars++;
        }
        if (use_mutexes) {
            for (std::size_t val = 1; val < task_proxy.get_operators().size() + 1; val++) {
                lp_var_ids[var.get_id()][val + var.get_domain_size()] = num_lp_vars++;
            }
        }
        fact_potentials[var.get_id()].resize(var.get_domain_size());
    }
    if (use_mutexes) {
        construct_mutex_lp();
    } else {
        construct_lp();
    }
}

bool PotentialOptimizer::has_optimal_solution() const {
    return lp_solver.has_optimal_solution();
}

void PotentialOptimizer::optimize_for_state(const State &state) {
    optimize_for_samples({state}
                         );
}

int PotentialOptimizer::get_lp_var_id(const FactProxy &fact) const {
    int var_id = fact.get_variable().get_id();
    int value = fact.get_value();
    assert(utils::in_bounds(var_id, lp_var_ids));
    assert(utils::in_bounds(value, lp_var_ids[var_id]));
    return lp_var_ids[var_id][value];
}

void PotentialOptimizer::optimize_for_all_states() {
    if (!potentials_are_bounded()) {
        cerr << "Potentials must be bounded for all-states LP." << endl;
        utils::exit_with(ExitCode::SEARCH_INPUT_ERROR);
    }
    vector<double> coefficients(num_lp_vars, 0.0);
    for (FactProxy fact : task_proxy.get_variables().get_facts()) {
        coefficients[get_lp_var_id(fact)] = 1.0 / fact.get_variable().get_domain_size();
    }
    lp_solver.set_objective_coefficients(coefficients);
    solve_and_extract();
    if (!has_optimal_solution()) {
        ABORT("all-states LP unbounded even though potentials are bounded.");
    }
}

void PotentialOptimizer::optimize_for_samples(const vector<State> &samples) {
    vector<double> coefficients(num_lp_vars, 0.0);
    for (const State &state : samples) {
        for (FactProxy fact : state) {
            coefficients[get_lp_var_id(fact)] += 1.0;
        }
    }
    lp_solver.set_objective_coefficients(coefficients);
    solve_and_extract();
}

// this method is used for mutex based potential heuristics.
// as partial states are used, no states are delivered, but tuples of variable id and value and their weight.
void PotentialOptimizer::   optimize_for_weighted_samples(
    const std::vector<tuple<int, int, long double>> &weights) {
    vector<double> coefficients(num_lp_vars, 0.0);
    // TODO assert that id and value are in range (like in line 61 and 62)
    for (const tuple<int, int, float> fact : weights) {
        coefficients[lp_var_ids[get<0>(fact)][get<1>(fact)]] += get<2>(fact);
    }
    lp_solver.set_objective_coefficients(coefficients);
    solve_and_extract();
}

const shared_ptr<AbstractTask> PotentialOptimizer::get_task() const {
    return task;
}

MutexTable * PotentialOptimizer::get_mutex_table() {
    return table;
}

bool PotentialOptimizer::potentials_are_bounded() const {
    return max_potential != numeric_limits<double>::infinity();
}

void PotentialOptimizer::construct_lp() {
    double upper_bound = (potentials_are_bounded() ? max_potential :
                          lp_solver.get_infinity());

    vector<lp::LPVariable> lp_variables;
    lp_variables.reserve(num_lp_vars);
    for (int lp_var_id = 0; lp_var_id < num_lp_vars; ++lp_var_id) {
        // Use dummy coefficient for now. Adapt coefficient later.
        lp_variables.emplace_back(-lp_solver.get_infinity(), upper_bound, 1.0);
    }

    vector<lp::LPConstraint> lp_constraints;
    for (OperatorProxy op : task_proxy.get_operators()) {
        // Create constraint:
        // Sum_{V in vars(eff(o))} (P_{V=pre(o)[V]} - P_{V=eff(o)[V]}) <= cost(o)
        unordered_map<int, int> var_to_precondition;
        for (FactProxy pre : op.get_preconditions()) {
            var_to_precondition[pre.get_variable().get_id()] = pre.get_value();
        }
        lp::LPConstraint constraint(-lp_solver.get_infinity(), op.get_cost());
        vector<pair<int, int>> coefficients;
        for (EffectProxy effect : op.get_effects()) {
            VariableProxy var = effect.get_fact().get_variable();
            int var_id = var.get_id();

            // Set pre to pre(op) if defined, otherwise to u = |dom(var)|.
            int pre = -1;
            auto it = var_to_precondition.find(var_id);
            if (it == var_to_precondition.end()) {
                pre = get_undefined_value(var);
            } else {
                pre = it->second;
            }

            int post = effect.get_fact().get_value();
            int pre_lp = lp_var_ids[var_id][pre];
            int post_lp = lp_var_ids[var_id][post];
            assert(pre_lp != post_lp);
            coefficients.emplace_back(pre_lp, 1);
            coefficients.emplace_back(post_lp, -1);
        }
        sort(coefficients.begin(), coefficients.end());
        for (const auto &coeff : coefficients)
            constraint.insert(coeff.first, coeff.second);
        lp_constraints.push_back(constraint);
    }

    /* Create full goal state. Use value |dom(V)| as "undefined" value
       for variables V undefined in the goal. */
    vector<int> goal(task_proxy.get_variables().size(), -1);
    for (FactProxy fact : task_proxy.get_goals()) {
        goal[fact.get_variable().get_id()] = fact.get_value();
    }
    for (VariableProxy var : task_proxy.get_variables()) {
        if (goal[var.get_id()] == -1)
            goal[var.get_id()] = get_undefined_value(var);
    }

    for (VariableProxy var : task_proxy.get_variables()) {
        /*
          Create constraint (using variable bounds): P_{V=goal[V]} = 0
          When each variable has a goal value (including the
          "undefined" value), this is equivalent to the goal-awareness
          constraint \sum_{fact in goal} P_fact <= 0. We can't set the
          potential of one goal fact to +2 and another to -2, but if
          all variables have goal values, this is not beneficial
          anyway.
        */
        int var_id = var.get_id();
        lp::LPVariable &lp_var = lp_variables[lp_var_ids[var_id][goal[var_id]]];
        lp_var.lower_bound = 0;
        lp_var.upper_bound = 0;

        int undef_val_lp = lp_var_ids[var_id][get_undefined_value(var)];
        for (int val = 0; val < var.get_domain_size(); ++val) {
            int val_lp = lp_var_ids[var_id][val];
            // Create constraint: P_{V=v} <= P_{V=u}
            // Note that we could eliminate variables P_{V=u} if V is
            // undefined in the goal.
            lp::LPConstraint constraint(-lp_solver.get_infinity(), 0);
            constraint.insert(val_lp, 1);
            constraint.insert(undef_val_lp, -1);
            lp_constraints.push_back(constraint);
        }
    }
    lp_solver.load_problem(lp::LPObjectiveSense::MAXIMIZE, lp_variables, lp_constraints);

    // add the constraint \sum_{f\in I} P(f) = h^P_I(I).
    if (initial_constraint) {
        utils::g_log << "Adding initial constraint." << endl;
        State init = task_proxy.get_initial_state();
        optimize_for_samples({init});
        float heuristic_value = 0;
        for (FactProxy fact : init) {
            int var_id = fact.get_variable().get_id();
            int value = fact.get_value();
            assert(utils::in_bounds(var_id, fact_potentials));
            assert(utils::in_bounds(value, fact_potentials[var_id]));
            heuristic_value += fact_potentials[var_id][value];
        }

        lp::LPConstraint constraint(heuristic_value, heuristic_value);
        for (FactProxy fact : init) {
            constraint.insert(lp_var_ids[fact.get_variable().get_id()][fact.get_value()], 1);
        }

        lp_solver.add_temporary_constraints({constraint});
    }
}

void PotentialOptimizer::construct_mutex_lp() {
    double upper_bound = (potentials_are_bounded() ? max_potential :
                          lp_solver.get_infinity());

    vector<lp::LPVariable> lp_variables;
    lp_variables.reserve(num_lp_vars);
    for (int lp_var_id = 0; lp_var_id < num_lp_vars; ++lp_var_id) {
        // Use dummy coefficient for now. Adapt coefficient later.
        lp_variables.emplace_back(-lp_solver.get_infinity(), upper_bound, 1.0);
    }

    vector<lp::LPConstraint> lp_constraints;
    for (OperatorProxy op : task_proxy.get_operators()) {
        // Create constraint:
        // Sum_{V in vars(eff(o))} (P_{V=pre(o)[V]} - P_{V=eff(o)[V]}) <= cost(o)
        unordered_map<int, int> var_to_precondition;
        for (FactProxy pre : op.get_preconditions()) {
            var_to_precondition[pre.get_variable().get_id()] = pre.get_value();
        }
        lp::LPConstraint constraint(-lp_solver.get_infinity(), op.get_cost());
        vector<pair<int, int>> coefficients;
        for (EffectProxy effect : op.get_effects()) {
            VariableProxy var = effect.get_fact().get_variable();
            int var_id = var.get_id();

            // Set pre to pre(op) if defined, otherwise to U_E^o_V = |dom(var)| + op.id.
            int pre = -1;
            auto it = var_to_precondition.find(var_id);
            if (it == var_to_precondition.end()) {
                pre = get_undefined_value_for_operator(var, op);
            } else {
                pre = it->second;
            }

            int post = effect.get_fact().get_value();
            int pre_lp = lp_var_ids[var_id][pre];
            int post_lp = lp_var_ids[var_id][post];
            assert(pre_lp != post_lp);
            coefficients.emplace_back(pre_lp, 1);
            coefficients.emplace_back(post_lp, -1);
        }
        sort(coefficients.begin(), coefficients.end());
        for (const auto &coeff : coefficients)
            constraint.insert(coeff.first, coeff.second);
        lp_constraints.push_back(constraint);
    }

    /* Create full goal state. Use value |dom(V)| as "undefined" value
       for variables V undefined in the goal. */
    vector<int> goal(task_proxy.get_variables().size(), -1);
    map<int, int> goal_map;
    for (FactProxy fact : task_proxy.get_goals()) {
        goal[fact.get_variable().get_id()] = fact.get_value();
        goal_map[fact.get_variable().get_id()] = fact.get_value();
    }
    vector<vector<int>> domains = table->multi_fact_disambiguation(goal_map);
    if (domains[0].empty()) {
        // problem unsolvable, if goal-state it contains a mutex.
        utils::exit_with(ExitCode::SEARCH_UNSOLVABLE);
    }
    for (VariableProxy var : task_proxy.get_variables()) {
        if (goal[var.get_id()] == -1) {
                goal[var.get_id()] = get_undefined_value(var);
        }
    }

    // add constraints on U_G_V
    for (VariableProxy var : task_proxy.get_variables()) {
        int var_id = var.get_id();
        lp::LPVariable &lp_var = lp_variables[lp_var_ids[var_id][goal[var_id]]];
        lp_var.lower_bound = 0;
        lp_var.upper_bound = 0;

        int undef_val_lp = lp_var_ids[var_id][get_undefined_value(var)];
        for (int val : domains[var_id]) {
            int val_lp = lp_var_ids[var_id][val];
            /*
             Create constraint (using variable bounds): P_{V=goal[V]} = 0
             When each variable has a goal value (including the
             "undefined goal" value), this is equivalent to the goal-awareness
             constraint \sum_{fact in goal} P_fact <= 0. We can't set the
             potential of one goal fact to +2 and another to -2, but if
             all variables have goal values, this is not beneficial
             anyway.
           */
            lp::LPConstraint constraint(-lp_solver.get_infinity(), 0);
            constraint.insert(val_lp, 1);
            constraint.insert(undef_val_lp, -1);
            lp_constraints.push_back(constraint);
        }
    }

    // add constraints for U_E^o_V
    for (OperatorProxy o : task_proxy.get_operators()) {
        map<int, int> pre;
        for(FactProxy fact : o.get_preconditions()) {
            pre[fact.get_variable().get_id()] = fact.get_value();
        }
        domains = table->multi_fact_disambiguation(pre);
        if (domains[0].empty()) {
            // if a precondition contains a mutex, this operation will never be on any reachable path.
            continue;
        }
        for (VariableProxy var : task_proxy.get_variables()) {
            int var_id = var.get_id();

            int undef_val_lp = lp_var_ids[var_id][get_undefined_value_for_operator(var, o)];
            for (int val : domains[var_id]) {
                int val_lp = lp_var_ids[var_id][val];
                // Create constraint: P_{V=v} <= P_{V=u}
                // Note that we could eliminate variables P_{V=u} if V is
                // undefined in the goal.
                lp::LPConstraint constraint(-lp_solver.get_infinity(), 0);
                constraint.insert(val_lp, 1);
                constraint.insert(undef_val_lp, -1);
                lp_constraints.push_back(constraint);
            }
        }
    }
    lp_solver.load_problem(lp::LPObjectiveSense::MAXIMIZE, lp_variables, lp_constraints);

    // add the constraint \sum_{f\in I} P(f) = h^P_I(I).
    if (initial_constraint) {
        utils::g_log << "Adding initial constraint." << endl;
        State init = task_proxy.get_initial_state();
        optimize_for_samples({init});
        float heuristic_value = 0;
        for (FactProxy fact : init) {
            int var_id = fact.get_variable().get_id();
            int value = fact.get_value();
            assert(utils::in_bounds(var_id, fact_potentials));
            assert(utils::in_bounds(value, fact_potentials[var_id]));
            heuristic_value += fact_potentials[var_id][value];
        }

        lp::LPConstraint constraint(heuristic_value, heuristic_value);
        for (FactProxy fact : init) {
            constraint.insert(lp_var_ids[fact.get_variable().get_id()][fact.get_value()], 1);
        }

        lp_solver.add_temporary_constraints({constraint});
    }
}

void PotentialOptimizer::solve_and_extract() {
    lp_solver.solve();
    if (has_optimal_solution()) {
        extract_lp_solution();
    }
}

void PotentialOptimizer::extract_lp_solution() {
    assert(has_optimal_solution());
    const vector<double> solution = lp_solver.extract_solution();
    for (FactProxy fact : task_proxy.get_variables().get_facts()) {
        fact_potentials[fact.get_variable().get_id()][fact.get_value()] =
            solution[get_lp_var_id(fact)];
    }
}

unique_ptr<PotentialFunction> PotentialOptimizer::get_potential_function() const {
    assert(has_optimal_solution());
    return utils::make_unique_ptr<PotentialFunction>(fact_potentials);
}
}
