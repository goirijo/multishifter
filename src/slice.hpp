#ifndef MAIN_SLICE_HH
#define MAIN_SLICE_HH

#include "./misc.hpp"
#include "casmutils/xtal/structure.hpp"
#include "multishift/slicer.hpp"
#include <filesystem>
#include <ostream>
#include <unordered_map>
#include <vector>

using Structure = casmutils::xtal::Structure;

void write_slicer_structures(const mush::Slicer& slicer, const mush::fs::path& slices_path, std::ostream& log);

std::unordered_map<std::string, MultiRecord> write_cleaver_structures(const std::vector<Structure>& cleaved_structures,
                                                                         const std::vector<double>& cleavage_values,
                                                                         const mush::fs::path& cleaved_root,
                                                                         std::ostream& log,
                                                                         const MultiRecord& starting_state);

#endif
