#include "gtest/gtest.h"
#include "casm/casm_io/jsonParser.hh"
#include "multishift/shift.hpp"

CASM::jsonParser generate_required_shift_settings_json()
{
    CASM::jsonParser shift_settings;
    shift_settings["slab"]="./this/is/not/real.vasp";

    return shift_settings;
}

CASM::jsonParser generate_shift_settings_json()
{
    CASM::jsonParser shift_settings;
    shift_settings["slab"]="./this/is/not/real.vasp";
    shift_settings["a"]=19;
    shift_settings["b"]=21;
    shift_settings["cleavage"]=std::vector<double>{-0.3,-0.1,0.0,0.1,0.3,0.5,1.0,2.0};

    return shift_settings;
}

TEST(ShiftSettingsTest, jsonRead)
{
    auto shift_json=generate_shift_settings_json();
    auto shift_settings=mush::ShiftSettings::from_json(shift_json);

    ASSERT_EQ(shift_settings.slab_path(),shift_json["slab"].get<std::string>());
    ASSERT_EQ(shift_settings.a_points(),shift_json["a"].get<int>());
    ASSERT_EQ(shift_settings.b_points(),shift_json["b"].get<int>());
    ASSERT_EQ(shift_settings.cleavage_values(),shift_json["cleavage"].get<std::vector<double>>());
}

TEST(ShiftSettingsTest, jsonReadDefault)
{
    auto shift_json=generate_required_shift_settings_json();
    auto shift_settings=mush::ShiftSettings::from_json(shift_json);

    ASSERT_EQ(shift_settings.slab_path(),shift_json["slab"].get<std::string>());
    ASSERT_EQ(shift_settings.a_points(),1);
    ASSERT_EQ(shift_settings.b_points(),1);
    ASSERT_EQ(shift_settings.cleavage_values(),std::vector<double>{0.0});
}

TEST(ShiftSettingsTest, jsonSerialize)
{
    auto shift_json=generate_shift_settings_json();
    auto shift_settings=mush::ShiftSettings::from_json(shift_json);
    auto rejson=shift_settings.to_json();

    ASSERT_EQ(shift_json, rejson);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
