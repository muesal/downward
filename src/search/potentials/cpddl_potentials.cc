#include "potential_function.h"
#include "potential_heuristic.h"
#include "potential_optimizer.h"
#include "util.h"

#include "../option_parser.h"
#include "../plugin.h"

#include "../utils/system.h"

using namespace std;

namespace potentials {

static unique_ptr<PotentialFunction> create_potential_function(Options &opts) {
    const AbstractTask &task = *opts.get<shared_ptr<AbstractTask>>("transform");
    TaskProxy task_proxy(task);
    std::vector<std::vector<double>> fact_potentials;
    int num_facts = 0;
    fact_potentials.resize(task_proxy.get_variables().size());
    for (int i = 0; i < fact_potentials.size(); ++i) {
        fact_potentials[i].resize(task_proxy.get_variables()[i].get_domain_size());
        num_facts += fact_potentials[i].size();
    }

    int var,val;
    double weight;
    string line;
    while (line != "begin_potentials") {
        getline(cin,line);
    }
    for (int i = 0; i < num_facts; ++i) {
        getline(cin,line);
        std::stringstream ss(line);
        ss >> var;
        ss >> val;
        ss >> weight;
        fact_potentials[var][val] = weight;
    }

    return utils::make_unique_ptr<PotentialFunction>(fact_potentials);
}

static shared_ptr<Heuristic> _parse(OptionParser &parser) {
    prepare_parser_for_admissible_potentials(parser);
    Options opts = parser.parse();
    if (parser.dry_run())
        return nullptr;

    return make_shared<PotentialHeuristic>(
        opts, create_potential_function(opts));
}

static Plugin<Evaluator> _plugin_cpddl_potentials(
    "cpddl_potentials", _parse, "heuristics_potentials");
}
