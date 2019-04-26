#include "gtest/gtest.h"
#include "casm/casm_io/jsonParser.hh"
#include "casm/casm_io/json_io/container.hh"
#include "multishift/fourier.hpp"

CASM::jsonParser generate_fourier_settings_json()
{
    CASM::jsonParser fourier_settings;

    fourier_settings["data"]="./this/is/not/real.vasp";
    fourier_settings["lattice"]="./this/is/also/not/real.vasp";
    fourier_settings["values"]=std::vector<std::string>{"value1","value2","value3"};

    return fourier_settings;
}

CASM::jsonParser generate_analytical_settings()
{
    CASM::jsonParser analytical_settings;

    analytical_settings["x"] = "X";
    analytical_settings["y"] = "Y";
    analytical_settings["library"] = "numpy";
    analytical_settings["minimum_magnitude"] = 1e-8;

    return analytical_settings;
}

TEST(FourierSettingsTest, FourierjsonRead)
{
    auto fourier_json=generate_fourier_settings_json();
    auto fourier_settings=mush::FourierSettings::from_json(fourier_json);

    ASSERT_EQ(fourier_settings.data_path(),fourier_json["data"].get<std::string>());
    ASSERT_EQ(fourier_settings.lattice_path(),fourier_json["lattice"].get<std::string>());
    ASSERT_EQ(fourier_settings.values_to_interpolate(),fourier_json["values"].get<std::vector<std::string>>());
}

TEST(FourierSettingsTest, AnalyticjsonRead)
{
    auto analytic_json=generate_analytical_settings();
    auto analytic_settings=mush::FourierAnalyticalSettings::from_json(analytic_json);

    ASSERT_EQ(analytic_settings.x(),analytic_json["x"].get<std::string>());
    ASSERT_EQ(analytic_settings.y(),analytic_json["y"].get<std::string>());
    ASSERT_EQ(analytic_settings.math_library(),analytic_json["library"].get<std::string>());
    ASSERT_EQ(analytic_settings.min_magnitude(),analytic_json["minimum_magnitude"].get<double>());
}

TEST(FourierSettingsTest, AnaliticjsonReadDefault)
{
    CASM::jsonParser empty_json;
    auto analytic_settings=mush::FourierAnalyticalSettings::from_json(empty_json);

    ASSERT_EQ(analytic_settings.x(),"xx");
    ASSERT_EQ(analytic_settings.y(),"yy");
    ASSERT_EQ(analytic_settings.math_library(),"np");
    ASSERT_EQ(analytic_settings.min_magnitude(),1e-12);
}

TEST(FourierSettingsTest, FourierjsonSerialize)
{
    auto fourier_json=generate_fourier_settings_json();
    auto fourier_settings=mush::FourierSettings::from_json(fourier_json);
    auto rejson=fourier_settings.to_json();

    ASSERT_EQ(fourier_json, rejson);
}

TEST(FourierSettingsTest, AnalyticjsonSerialize)
{
    auto analytic_json=generate_analytical_settings();
    auto analytic_settings=mush::FourierAnalyticalSettings::from_json(analytic_json);
    auto rejson=analytic_settings.to_json();

    ASSERT_EQ(analytic_json, rejson);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
