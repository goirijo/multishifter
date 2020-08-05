#ifndef MAIN_MISC_HH
#define MAIN_MISC_HH

#include <multishift/definitions.hpp>
#include <multishift/slicer.hpp>
#include <casmutils/xtal/structure.hpp>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using json = nlohmann::json;
using Structure = mush::cu::xtal::Structure;

enum class COMMAND
{
    SHIFT,
    CLEAVE,
    TWIST
};

/**
 * Keeps all the info for a particular structure, cleavage,
 * shifts, and twists.
 */

struct MultiRecord
{
    double cleavage=0.0;
    int a_index=0;
    int b_index=0;
    /* double a_frac=0; */
    /* double b_frac=0; */
    /* double a_cart=0; */
    /* double b_cart=0; */
    /* double angle=0.0; */
    std::string id() const;
    std::vector<std::string> equivalent_structures;
};


json load_json(const mush::fs::path& json_path);
void write_json(const json& json, const mush::fs::path& target);

///If directory already exists, throw exception, otherwise continue normally
void cautious_create_directory(const mush::fs::path new_dir);

std::string extract_name_from_settings(const json& settings);
/* json extract_slice_subsettings(const json& settings); */

///Converts map from path to MultiRecord into the final json format, which
///might include more information than what's in the MultiRecord objects
json record_to_json(const std::unordered_map<std::string,MultiRecord>& record);

#endif
