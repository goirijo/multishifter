#include "./cleave.hpp"
#include "./misc.hpp"
#include <casmutils/xtal/structure_tools.hpp>
#include <casmutils/mush/shift.hpp>
#include <multishift/slice_settings.hpp>


void setup_subcommand_cleave(CLI::App& app)
{
    auto settings_path_ptr = std::make_shared<mush::fs::path>();

    CLI::App* cleave_sub = app.add_subcommand("cleave", "Cleave unit cell to insert/remove empty space at the interface layer.");
    cleave_sub
        ->add_option("-s,--settings", *settings_path_ptr, "Settings file with path to slab unit, slab thickness, and cleavage values.")
        ->required();

    cleave_sub->callback([settings_path_ptr]() { run_subcommand_cleave(*settings_path_ptr, std::cout); });
}

void run_subcommand_cleave(const mush::fs::path& settings_path, std::ostream& log)
{
}

