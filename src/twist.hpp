#ifndef TWIST_SUBCOMMAND_HH
#define TWIST_SUBCOMMAND_HH

#include <CLI/CLI.hpp>
#include <multishift/definitions.hpp>

void setup_subcommand_twist(CLI::App& app);
void run_subcommand_twist(const mush::fs::path& settings_path, std::ostream& log);

#endif
