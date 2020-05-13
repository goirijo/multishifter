#include "../../autotools.hh"
#include <filesystem>
#include <fstream>
#include <multishift/slice_settings.hpp>

#include <gtest/gtest.h>
#include <memory>

using namespace mush;

class LoadSliceSettingsFromJson : public testing::Test
{
    protected:
        std::unique_ptr<SliceSettings> settings_ptr;
        fs::path settings_path;
        virtual void SetUp() override
        {
            settings_path=autotools::input_filesdir/"slice_settings.json";
            std::ifstream settings_stream(settings_path);
            json j;
            settings_stream>>j;
            settings_ptr.reset(new SliceSettings(SliceSettings::from_json(j)));
        }

};

TEST_F(LoadSliceSettingsFromJson, PrimPath)
{
    EXPECT_EQ("./b2.vasp",settings_ptr->prim_path);
}

TEST_F(LoadSliceSettingsFromJson, MillerIndexes)
{
    EXPECT_EQ(settings_ptr->miller_indexes,Eigen::Vector3i(1,-1,2));
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
