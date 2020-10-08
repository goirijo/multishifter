#include "./shift.hpp"
#include "multishift/shifter.hpp"
#include <casmutils/xtal/structure_tools.hpp>
#include <multishift/slice_settings.hpp>

void setup_subcommand_shift(CLI::App& app)
{
    auto settings_path_ptr = std::make_shared<mush::fs::path>();

    CLI::App* shift_sub = app.add_subcommand(
        "shift", "Shift unit cell along interface to generate stacking faults and structures for gamma surface calculations.");
    shift_sub->add_option("-s,--settings", *settings_path_ptr, "Settings file with path to slab unit, slab thickness, and grid density.")
        ->required();

    shift_sub->callback([settings_path_ptr]() { run_subcommand_shift(*settings_path_ptr, std::cout); });
}

void run_subcommand_shift(const mush::fs::path& settings_path, std::ostream& log)
{
}
