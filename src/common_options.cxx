#include "./common_options.hpp"

void populate_subcommand_output_option(CLI::App* sub, mush::fs::path* out, bool require)
{
    sub->add_option("-o,--output",*out,"Target output file")->required(require);
}

void populate_subcommand_input_option(CLI::App* sub, mush::fs::path* in, bool require)
{
    sub->add_option("-i,--input",*in,"Source input file")->required(require);
}

void populate_subcommand_fractional(CLI::App* sub, bool* frac_ptr, CLI::Option* needed)
{
    sub->add_flag("--fractional", *frac_ptr, "Specifies that the parameters passed to "+needed->get_name()+" are in fractional coordinates, not Cartesian")->needs(needed);
}
