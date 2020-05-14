#ifndef MAIN_SLICE_HH
#define MAIN_SLICE_HH

#include "./misc.hpp"
#include "multishift/slicer.hpp"
#include <ostream>

void write_slicer_structures(const mush::Slicer& slicer, const mush::fs::path& slices_path, std::ostream& log);

#endif
