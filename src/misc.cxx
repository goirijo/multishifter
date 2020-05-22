#include "./misc.hpp"
#include <fstream>
#include <stdexcept>
#include <string>

std::string MultiRecord::id() const
{
    std::string id = std::to_string(a_index) + ":" + std::to_string(b_index) + ":";

    std::stringstream cleavestream;
    cleavestream << std::fixed << std::setprecision(6) << cleavage;

    id+=cleavestream.str()+":";

    std::stringstream twiststream;
    twiststream << std::fixed << std::setprecision(6) << angle;

    id+=twiststream.str();
    return id;
}

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
    json_stream << json;
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

std::string extract_name_from_settings(const json& settings)
{
    if (settings.count("name") == 0)
    {
        throw std::runtime_error("Could not find entry 'name' in settings file.");
    }
    return settings["name"];
}
/* json extract_slice_subsettings(const json& settings) { return settings["slice"];} */

json record_to_json(const std::unordered_map<std::string, MultiRecord>& record)
{
    //TODO: Maybe it's better to have the path as a value, and the id as a key?
    json j;
    for (const auto& [path, mr] : record)
    {
        json mr_json;
        mr_json["cleavage"] = mr.cleavage;
        mr_json["a_index"] = mr.a_index;
        mr_json["b_index"] = mr.b_index;
        mr_json["angle"] = mr.angle;
        mr_json["id"] = mr.id();
        mr_json["equivalent_structures"] = mr.equivalent_structures;

        if(mr.equivalent_structures.size()==0)
        {
            mr_json["equivalent_structures"].push_back(mr_json["id"]);
        }

        j[path] = mr_json;
    }
    return j;
}

