#ifndef CLEAVE_SUBCOMMAND_HH
#define CLEAVE_SUBCOMMAND_HH

#include "misc.hpp"
#include <multishift/definitions.hpp>
#include <casmutils/xtal/structure.hpp>
#include <CLI/CLI.hpp>


void setup_subcommand_cleave(CLI::App& app);
void run_subcommand_cleave(const mush::fs::path& settings_path, std::ostream& log);
#endif
