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

namespace mush
{

/**
 * Keeps all the info for a particular structure, cleavage
 */

struct MultiRecord
{
    double cleavage=0.0;
    int a_index=0;
    int b_index=0;
    std::string id() const;
    std::vector<std::string> equivalent_structures;
};

enum class SUBCOMMAND {CLEAVE,SHIFT,CHAIN};

std::string make_cleave_dirname(double cleavage);
std::string make_shift_dirname(int a, int b);
template<SUBCOMMAND>
mush::fs::path make_target_directory(const mush::MultiRecord& record);


json load_json(const mush::fs::path& json_path);
void write_json(const json& json, const mush::fs::path& target);

///If directory already exists, throw exception, otherwise continue normally
void cautious_create_directory(const mush::fs::path new_dir);

}
#endif
