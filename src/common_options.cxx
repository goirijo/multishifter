#include "./common_options.hpp"

void populate_subcommand_output_option(CLI::App* sub, mush::fs::path* out)
{
    sub->add_option("-o,--output",*out,"Target output file or directory.")->required();//->check(CLI::NonexistentPath);
}

void populate_subcommand_input_option(CLI::App* sub, mush::fs::path* in)
{
    sub->add_option("-i,--input",*in,"Source input structure or slab file.")->required();//->check(CLI::ExistingFile);
}

void populate_subcommand_fractional(CLI::App* sub, bool* frac_ptr, CLI::Option* needed)
{
    sub->add_flag("--fractional", *frac_ptr, "Specifies that the parameters passed to "+needed->get_name()+" are in fractional coordinates, not Cartesian.")->needs(needed);
}
