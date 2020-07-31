#ifndef STACK_SUBCOMMAND_HH
#define STACK_SUBCOMMAND_HH

#include <CLI/CLI.hpp>
#include <filesystem>
#include <multishift/definitions.hpp>

void setup_subcommand_stack(CLI::App& app);
void run_subcommand_stack(const std::vector<mush::fs::path>& input_paths, const mush::fs::path& output_path, std::ostream& log);

#endif
