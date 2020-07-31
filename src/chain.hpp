#ifndef CHAIN_SUBCOMMAND_HH
#define CHAIN_SUBCOMMAND_HH

#include <CLI/CLI.hpp>
#include <multishift/definitions.hpp>

void setup_subcommand_chain(CLI::App& app);
void run_subcommand_chain(const mush::fs::path& settings_path, std::ostream& log);

#endif
