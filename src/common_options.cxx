#include "./common_options.hpp"

void populate_subcommand_output_option(CLI::App* sub, mush::fs::path* out, bool require)
{
    sub->add_option("-o,--oputput",*out,"Target output file")->required(require);
}

void populate_subcommand_input_option(CLI::App* sub, mush::fs::path* in, bool require)
{
    sub->add_option("-i,--input",*in,"Source input file")->required(require);
}
