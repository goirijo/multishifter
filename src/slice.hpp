#ifndef MAIN_SLICE_HH
#define MAIN_SLICE_HH

#include "./misc.hpp"
#include "casmutils/xtal/structure.hpp"
#include <casmutils/mush/shift.hpp>
#include "multishift/slicer.hpp"
#include <filesystem>
#include <ostream>
#include <unordered_map>
#include <vector>
#include <CLI/CLI.hpp>

using Structure = casmutils::xtal::Structure;

/* void write_slicer_structures(const mush::Slicer& slicer, const mush::fs::path& slices_path, std::ostream& log); */

void setup_subcommand_slice(CLI::App& app);
void run_subcommand_slice(const mush::fs::path& input_path, const mush::fs::path& output_path, const std::vector<int>& millers, bool align, std::ostream& log);

#endif
