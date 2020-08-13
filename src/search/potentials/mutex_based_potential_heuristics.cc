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
using Tuple = std::vector<FactPair>;

namespace potentials {

    /* TODO: these need a class partial_state...*/

    static void single_fact_diasambiguation(State &state){
        //while state changes check for single fact disambiguations and assign them
        state.get_task(); //placeholder
    }

    static void multi_fact_diasambiguation(State &state){
        //while the set changes, check for disambiguations.
        state.get_task(); //placeholder
    } 

    static unique_ptr<PotentialFunction> create_mutex_based_potential_function(
            const Options &opts) {
        PotentialOptimizer optimizer(opts);
        const AbstractTask &task = *opts.get<shared_ptr<AbstractTask>>("transform");
        TaskProxy task_proxy(task);

        //get all pairs which are mutex
        auto hm = make_shared<HMHeuristic>(opts);
        State initial = task_proxy.get_initial_state();
        //TODO: get the table
        //GlobalState gs_initial;
        //hm->compute_heuristic(task.);
        std::map<Tuple, int> mutexes;
        for (pair<Tuple, int> pair: hm->hm_table) {
            if (pair.second == numeric_limits<int>::max())
                mutexes.insert(pair);
        }
        hm.reset();

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