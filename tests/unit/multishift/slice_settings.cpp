#include "../../autotools.hh"
#include <filesystem>
#include <fstream>
#include <multishift/slice_settings.hpp>

#include <gtest/gtest.h>
#include <memory>

using namespace mush;

class LoadSettingsFromJson : public testing::Test
{
    protected:
        std::unique_ptr<SliceSettings> slice_settings_ptr;
        std::unique_ptr<SlabSettings> slab_settings_ptr;
        std::unique_ptr<CleavageSettings> cleavage_settings_ptr;
        std::unique_ptr<ShiftSettings> shift_settings_ptr;

        fs::path settings_path;
        virtual void SetUp() override
        {
            settings_path=autotools::input_filesdir/"complete_settings.json";
            std::ifstream settings_stream(settings_path);
            json j;
            settings_stream>>j;

            slice_settings_ptr.reset(new SliceSettings(SliceSettings::from_json(j)));
            slab_settings_ptr.reset(new SlabSettings(SlabSettings::from_json(j)));
            cleavage_settings_ptr.reset(new CleavageSettings(CleavageSettings::from_json(j)));
            shift_settings_ptr.reset(new ShiftSettings(ShiftSettings::from_json(j)));
        }

};

TEST_F(LoadSettingsFromJson, SlicePrimPath)
{
    EXPECT_EQ("./b2.vasp",slice_settings_ptr->prim_path);
}

TEST_F(LoadSettingsFromJson, SliceMillerIndexes)
{
    EXPECT_EQ(slice_settings_ptr->miller_indexes,Eigen::Vector3i(1,-1,2));
}

TEST_F(LoadSettingsFromJson, SlabUnitPath)
{
    EXPECT_EQ("./b2_floor0.vasp",slab_settings_ptr->slab_unit_path);
}

TEST_F(LoadSettingsFromJson, SlabStacks)
{
    EXPECT_EQ(6,slab_settings_ptr->stacks);
}

TEST_F(LoadSettingsFromJson, ShiftGrid)
{
    EXPECT_EQ(shift_settings_ptr->a_density,6);
    EXPECT_EQ(shift_settings_ptr->b_density,8);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
