#include "./shift.hpp"
#include "./common_options.hpp"
#include "./chain.hpp"
#include "multishift/shifter.hpp"
#include <casmutils/xtal/structure_tools.hpp>
#include <multishift/slice_settings.hpp>

template<>
mush::fs::path mush::make_target_directory<mush::SUBCOMMAND::SHIFT>(const mush::MultiRecord& record)
{
    return mush::make_shift_dirname(record.a_index,record.b_index);
}

void setup_subcommand_shift(CLI::App& app)
{
    auto input_path_ptr = std::make_shared<mush::fs::path>();
    auto output_path_ptr = std::make_shared<mush::fs::path>();
    auto grid_dims_ptr = std::make_shared<std::vector<int>>();

    CLI::App* shift_sub = app.add_subcommand("shift", "Shift slabs parallel to each other at regular intervals.");

    populate_subcommand_input_option(shift_sub, input_path_ptr.get());
    populate_subcommand_output_option(shift_sub, output_path_ptr.get());

        shift_sub->add_option("-g,--grid",
                     *grid_dims_ptr,
                     "Grid dimensions to divide the ab-lane of the slab into. Each grid point will correspond to a particular shift "
                     "vector. Periodic image not counted in grid.")
        ->expected(2)
        ->required();

    shift_sub->callback([=]() { run_subcommand_chain<mush::SUBCOMMAND::SHIFT>(*input_path_ptr, *output_path_ptr, {0.0}, *grid_dims_ptr, std::cout); });
}
