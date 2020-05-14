#ifndef MAIN_MISC_HH
#define MAIN_MISC_HH

#include <multishift/definitions.hpp>
#include <multishift/slab.hpp>
#include <multishift/slicer.hpp>
#include <casmutils/xtal/structure.hpp>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>

using json = nlohmann::json;
using Structure=mush::cu::xtal::Structure;

json load_json(const mush::fs::path& json_path);
void write_json(const json& json, const mush::fs::path& target);

///If directory already exists, throw exception, otherwise continue normally
void cautious_create_directory(const mush::fs::path new_dir);

std::string extract_name_from_settings(const json& settings);
json extract_slice_subsettings(const json& settings);

#endif
