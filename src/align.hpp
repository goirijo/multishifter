#ifndef ALIGN_COMMAND_HH
#define ALIGN_COMMAND_HH

#include <CLI/CLI.hpp>
#include <filesystem>
#include <multishift/definitions.hpp>

void setup_subcommand_align(CLI::App& app);
void run_subcommand_align(const mush::fs::path& input_paths, const mush::fs::path& output_path, bool prismatic, std::ostream& log);

#endif
