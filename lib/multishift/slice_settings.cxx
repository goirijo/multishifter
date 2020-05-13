#include "./slice_settings.hpp"
#include <filesystem>

namespace mush
{
        SliceSettings SliceSettings::from_json(const json& input_settings)
        {
            fs::path prim_path(input_settings["prim"]);
            std::vector<int> millers(input_settings["miller_indexes"]);

            return SliceSettings(prim_path,Eigen::Vector3i(millers[0],millers[1],millers[2]));
        }
}
