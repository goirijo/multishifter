#include "./definitions.hpp"
#include "./slab.hpp"
#include <nlohmann/json.hpp>

namespace mush
{
    /**
     * Simple structure to hold settings required to slice a given structure
     * through a particular set of Miller indexes. The class can be used to generate
     * an approriately oriented sliced structure, which is the starting point
     * for generating shifts and twists.
     */

    struct SliceSettings
    {
        SliceSettings(const fs::path& prim_path, const Eigen::Vector3i& miller_indexes): prim_path(prim_path),miller_indexes(miller_indexes){}

        static SliceSettings from_json(const json& input_settings);

        fs::path prim_path;
        Eigen::Vector3i miller_indexes;
    };
};
