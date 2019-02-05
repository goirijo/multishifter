#include "gtest/gtest.h"
#include "casm/casm_io/jsonParser.hh"
#include "multishift/base.hpp"

CASM::jsonParser generate_required_base_settings_json()
{
    CASM::jsonParser base_settings;
    base_settings["prim"]="./this/is/not/real.vasp";
    base_settings["millers"]=std::vector<int>{2,3,4};

    return base_settings;
}

CASM::jsonParser generate_base_settings_json()
{
    CASM::jsonParser base_settings;
    base_settings["prim"]="./this/is/not/real.vasp";
    base_settings["millers"]=std::vector<int>{2,3,4};
    base_settings["floor_slab_index"]=5;
    base_settings["stacks"]=10;

    return base_settings;
}

TEST(BaseSettingsTest, jsonRead)
{
    auto base_json=generate_base_settings_json();
    auto base_settings=mush::BaseSettings::from_json(base_json);

    ASSERT_EQ(base_settings.prim_path(),base_json["prim"].get<std::string>());
    ASSERT_EQ(base_settings.millers(),base_json["millers"].get<Eigen::Vector3i>());
    ASSERT_EQ(base_settings.floor_slab_atom_index(),base_json["floor_slab_index"].get<int>());
    ASSERT_EQ(base_settings.stacks(),base_json["stacks"].get<int>());
}

TEST(BaseSettingsTest, jsonReadDefault)
{
    auto base_json=generate_required_base_settings_json();
    auto base_settings=mush::BaseSettings::from_json(base_json);

    ASSERT_EQ(base_settings.prim_path(),base_json["prim"].get<std::string>());
    ASSERT_EQ(base_settings.millers(),base_json["millers"].get<Eigen::Vector3i>());
    ASSERT_EQ(base_settings.floor_slab_atom_index(),0);
    ASSERT_EQ(base_settings.stacks(),1);
}

TEST(BaseSettingsTest, jsonSerialize)
{
    auto base_json=generate_base_settings_json();
    auto base_settings=mush::BaseSettings::from_json(base_json);
    auto rejson=base_settings.to_json();

    ASSERT_EQ(base_json, rejson);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
