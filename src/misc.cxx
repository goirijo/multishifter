#include "./misc.hpp"
#include <fstream>

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
    json_stream<<json;
    return;
}

void cautious_create_directory(const mush::fs::path new_dir)
{
    if(mush::fs::exists(new_dir))
    {
        throw std::runtime_error("Will not continue because "+new_dir.string()+" already exists.");
    }

    mush::fs::create_directory(new_dir);
    return;
}

std::string extract_name_from_settings(const json& settings) { return settings["name"]; }
json extract_slice_subsettings(const json& settings) { return settings["slice"];}

