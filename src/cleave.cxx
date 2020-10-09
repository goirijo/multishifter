#include "./cleave.hpp"
#include "./misc.hpp"
#include <casmutils/xtal/structure_tools.hpp>
#include "./chain.hpp"
#include "./common_options.hpp"
#include <casmutils/mush/shift.hpp>
#include <multishift/slice_settings.hpp>

template<>
mush::fs::path mush::make_target_directory<mush::SUBCOMMAND::CLEAVE>(const mush::MultiRecord& record)
{
    return mush::make_cleave_dirname(record.cleavage);
}

void setup_subcommand_cleave(CLI::App& app)
{
    auto input_path_ptr = std::make_shared<mush::fs::path>();
    auto output_path_ptr = std::make_shared<mush::fs::path>();
    auto celavages_ptr = std::make_shared<std::vector<double>>();

    CLI::App* chain_sub = app.add_subcommand("cleave", "Create slab structures separated by a range of specified values in Angstrom.");

    populate_subcommand_input_option(chain_sub, input_path_ptr.get());
    populate_subcommand_output_option(chain_sub, output_path_ptr.get());

    chain_sub->add_option("-v,--values", *celavages_ptr, "List of cleavage values to insert between slabs.")->required();

    chain_sub->callback([=]() { run_subcommand_chain<mush::SUBCOMMAND::CLEAVE>(*input_path_ptr, *output_path_ptr, *celavages_ptr, {1,1}, std::cout); });
}

