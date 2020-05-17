#ifndef SLICE_SETTINGS_HH
#define SLICE_SETTINGS_HH

#include "./definitions.hpp"
#include "./slab.hpp"
#include <nlohmann/json.hpp>
#include <vector>

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

    /**
     * Simple structure to hold settings for creating a slab.
     * It's just the path to a starting slab unit (i.e. prim that
     * has already been sliced), and the number of times that
     * unit should be repeated along the c direction.
     */

    struct SlabSettings
    {
        SlabSettings(const fs::path& slab_unit_path, int stacks): slab_unit_path(slab_unit_path), stacks(stacks){};

        static SlabSettings from_json(const json& input_settings);

        fs::path slab_unit_path;
        int stacks;
    };

    /**
     * Holds the settings that define the cleavege values
     */

    struct CleavageSettings
    {
        CleavageSettings(const std::vector<double>& cleavage_values):cleavage_values(cleavage_values){}

        static CleavageSettings from_json(const json& input_settings);

        std::vector<double> cleavage_values;
    };

    /**
     * Holds settings that define the grid density in the a-b plane.
     * The grid must be uniform, so only two integers are needed: the
     * density along the a and b vectors.
     */

    struct ShiftSettings
    {
        ShiftSettings(int a_density, int b_density):a_density(a_density),b_density(b_density){}

        static ShiftSettings from_json(const json& input_settings);

        int a_density;
        int b_density;
    };
};

#endif
