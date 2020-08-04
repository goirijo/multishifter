#ifndef TRANSLATE_COMMAND_HH
#define TRANSLATE_COMMAND_HH

#include <CLI/CLI.hpp>
#include <filesystem>
#include <multishift/definitions.hpp>

void setup_subcommand_translate(CLI::App& app);
void run_subcommand_translate(const mush::fs::path& settings_path, const CLI::App& sub, std::ostream& log);

#endif
