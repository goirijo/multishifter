#include "./misc.hpp"
#include <fstream>
#include <stdexcept>
#include <string>

namespace mush
{
std::string MultiRecord::id() const
{
    std::string id = std::to_string(a_index) + ":" + std::to_string(b_index) + ":";

    std::stringstream cleavestream;
    cleavestream << std::fixed << std::setprecision(6) << cleavage;

    id+=cleavestream.str();
    return id;
}

std::string make_cleave_dirname(double cleavage)
{
    std::stringstream cleavestream;
    cleavestream << std::fixed << std::setprecision(6) << cleavage;
    return "cleave__" + cleavestream.str();
}

std::string make_shift_dirname(int a, int b) { return "shift__" + std::to_string(a) + "." + std::to_string(b); }


json load_json(const mush::fs::path& json_path)
{
    std::ifstream settings_stream(json_path);
    json j;
    settings_stream >> j;
    return j;
}

void write_json(const json& json, const mush::fs::path& target)
{
    std::ofstream json_stream(target);
    json_stream << json.dump(4);
    return;
}

void cautious_create_directory(const mush::fs::path new_dir)
{
    if (mush::fs::exists(new_dir))
    {
        throw std::runtime_error("Will not continue because " + new_dir.string() + " already exists.");
    }

    mush::fs::create_directory(new_dir);
    return;
}

}

