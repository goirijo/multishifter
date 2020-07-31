#include "./slice_settings.hpp"
#include <filesystem>
#include <vector>

namespace mush
{
        SliceSettings SliceSettings::from_json(const json& input_settings)
        {
            fs::path prim_path(input_settings["prim"]);
            std::vector<int> millers(input_settings["miller_indexes"]);

            return SliceSettings(prim_path,Eigen::Vector3i(millers[0],millers[1],millers[2]));
        }

        SlabSettings SlabSettings::from_json(const json& input_settings)
        {
            return SlabSettings(input_settings["slab_unit"],input_settings["stacks"]);
        }

        CleavageSettings CleavageSettings::from_json(const json& input_settings)
        {
            return CleavageSettings(input_settings["cleave"]);
        }

        ShiftSettings ShiftSettings::from_json(const json& input_settings)
        {
            std::vector<int> ab=input_settings["shift_grid"];
            return ShiftSettings(ab[0],ab[1]);
        }

        FourierSettings FourierSettings::from_json(const json& input_settings)
        {
            fs::path data_path(input_settings["data"]);
            double resolution=0.0005;
            if(input_settings.contains("k_resolution"))
            {
            resolution=input_settings["k_resolution"];
            }
            return FourierSettings(data_path,resolution);
        }

        TwisterSettings TwisterSettings::from_json(const json& input_settings)
        {
            std::vector<double> angles=input_settings["angles"];
            return TwisterSettings(angles);
        }
}
