#include <CLI/CLI.hpp>
#include <ostream>
#include "./misc.hpp"
#include "./slice.hpp"
#include "./cleave.hpp"
#include "./shift.hpp"
#include "./chain.hpp"
#include "./fourier.hpp"
#include "./twist.hpp"
#include "./stack.hpp"
#include "./mutate.hpp"
#include "./translate.hpp"
#include "./align.hpp"

int main(int argc, char** argv)
{
    mush::fs::path settings_path;

    CLI::App app{"Collection of tools that generate structures for crystal surface/interface calculaions."};

    setup_subcommand_slice(app);
    setup_subcommand_stack(app);
    setup_subcommand_translate(app);
    setup_subcommand_align(app);
    setup_subcommand_mutate(app);
    setup_subcommand_chain(app);
    /* setup_subcommand_cleave(app); */
    /* setup_subcommand_shift(app); */
    /* setup_subcommand_fourier(app); */
    /* setup_subcommand_twist(app); */

    app.require_subcommand();

    CLI11_PARSE(app, argc, argv);

    return 0;
}
