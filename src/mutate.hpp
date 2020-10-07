#ifndef MUTATE_SUBCOMMAND_HH
#define MUTATE_SUBCOMMAND_HH

#include <CLI/CLI.hpp>
#include <filesystem>
#include <multishift/definitions.hpp>

void setup_subcommand_mutate(CLI::App& app);
void run_subcommand_mutate(const mush::fs::path& input_path, const mush::fs::path& output_path, const std::vector<double>& mutation, bool frac, std::ostream& log);

#endif

