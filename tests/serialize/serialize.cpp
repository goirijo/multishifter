#include "gtest/gtest.h"
#include "casm/casm_io/jsonParser.hh"
#include "casm/crystallography/Lattice.hh"
#include "multishift/surf.hpp"
#include "multishift/misc.hpp"
#include "multishift/fourier.hpp"
#include <complex>

mush::Interpolator fake_interpolator()
{
    int a_div=13;
    int b_div=14;

    std::vector<double> af,bf,vs;
    for(int a=0; a<a_div; ++a)
    {
        for(int b=0; b<b_div; ++b)
        {
            af.push_back(a*1.0/a_div);
            bf.push_back(b*1.0/b_div);
            vs.push_back(std::sqrt(a*b));
        }
    }

    Eigen::Matrix3d fake_lat_mat;
    fake_lat_mat<<1.59,-2.76,0.0,-1.59,-2.76,5.18,-15.96,0.0,0.0;

    return mush::Interpolator(CASM::Lattice(fake_lat_mat),af,bf,vs);
}

TEST(Serialization, InterPoint)
{
    mush::InterPoint ip_raw(0.5,0.25,std::complex<double>(1.3,-5.5));
    auto ip_serial=ip_raw.serialize();
    auto ip_rebuild=mush::InterPoint::deserialize(ip_serial);

    ASSERT_TRUE(lazy::almost_equal<double>(ip_raw.a_frac,ip_rebuild.a_frac));
    ASSERT_TRUE(lazy::almost_equal<double>(ip_raw.b_frac,ip_rebuild.b_frac));
    ASSERT_TRUE(lazy::almost_equal<double>(ip_raw.value.real(),ip_rebuild.value.real()));
    ASSERT_TRUE(lazy::almost_equal<double>(ip_raw.value.imag(),ip_rebuild.value.imag()));
}

TEST(Serialization, InterpolatorLattices)
{
    auto ipolator=fake_interpolator();
    auto serialized=ipolator.serialize();
    auto reconstructed=mush::Interpolator::deserialize(serialized);

    ASSERT_TRUE(lazy::almost_equal<double>(ipolator.real_lattice()[0](0),reconstructed.real_lattice()[0](0)));
    ASSERT_TRUE(lazy::almost_equal<double>(ipolator.real_lattice()[0](1),reconstructed.real_lattice()[0](1)));
    ASSERT_TRUE(lazy::almost_equal<double>(ipolator.real_lattice()[0](2),reconstructed.real_lattice()[0](2)));

    ASSERT_TRUE(lazy::almost_equal<double>(ipolator.real_lattice()[1](0),reconstructed.real_lattice()[1](0)));
    ASSERT_TRUE(lazy::almost_equal<double>(ipolator.real_lattice()[1](1),reconstructed.real_lattice()[1](1)));
    ASSERT_TRUE(lazy::almost_equal<double>(ipolator.real_lattice()[1](2),reconstructed.real_lattice()[1](2)));

    ASSERT_TRUE(lazy::almost_equal<double>(ipolator.real_lattice()[2](0),reconstructed.real_lattice()[2](0)));
    ASSERT_TRUE(lazy::almost_equal<double>(ipolator.real_lattice()[2](1),reconstructed.real_lattice()[2](1)));
    ASSERT_TRUE(lazy::almost_equal<double>(ipolator.real_lattice()[2](2),reconstructed.real_lattice()[2](2)));


    ASSERT_TRUE(lazy::almost_equal<double>(ipolator.reciprocal_lattice()[0](0),reconstructed.reciprocal_lattice()[0](0)));
    ASSERT_TRUE(lazy::almost_equal<double>(ipolator.reciprocal_lattice()[0](1),reconstructed.reciprocal_lattice()[0](1)));
    ASSERT_TRUE(lazy::almost_equal<double>(ipolator.reciprocal_lattice()[0](2),reconstructed.reciprocal_lattice()[0](2)));

    ASSERT_TRUE(lazy::almost_equal<double>(ipolator.reciprocal_lattice()[1](0),reconstructed.reciprocal_lattice()[1](0)));
    ASSERT_TRUE(lazy::almost_equal<double>(ipolator.reciprocal_lattice()[1](1),reconstructed.reciprocal_lattice()[1](1)));
    ASSERT_TRUE(lazy::almost_equal<double>(ipolator.reciprocal_lattice()[1](2),reconstructed.reciprocal_lattice()[1](2)));

    ASSERT_TRUE(lazy::almost_equal<double>(ipolator.reciprocal_lattice()[2](0),reconstructed.reciprocal_lattice()[2](0)));
    ASSERT_TRUE(lazy::almost_equal<double>(ipolator.reciprocal_lattice()[2](1),reconstructed.reciprocal_lattice()[2](1)));
    ASSERT_TRUE(lazy::almost_equal<double>(ipolator.reciprocal_lattice()[2](2),reconstructed.reciprocal_lattice()[2](2)));
}

TEST(Serialization, InterpolatorPoints)
{
    auto ipolator=fake_interpolator();
    auto serialized=ipolator.serialize();
    auto reconstructed=mush::Interpolator::deserialize(serialized);

    ASSERT_EQ(ipolator.dims().first, reconstructed.dims().first);
    ASSERT_EQ(ipolator.dims().second, reconstructed.dims().second);

    for(int a=0; a<ipolator.dims().first; ++a)
    {
        for(int b=0; b<ipolator.dims().second; ++b)
        {
            const auto& r_orig=ipolator.sampled_values()[a][b];
            const auto& k_orig=ipolator.k_values()[a][b];

            const auto& r_reco=reconstructed.sampled_values()[a][b];
            const auto& k_reco=reconstructed.k_values()[a][b];

            ASSERT_TRUE(lazy::almost_equal<double>(r_orig.a_frac,r_reco.a_frac));
            ASSERT_TRUE(lazy::almost_equal<double>(r_orig.b_frac,r_reco.b_frac));
            ASSERT_TRUE(lazy::almost_equal<double>(r_orig.value.real(),r_reco.value.real()));
            ASSERT_TRUE(lazy::almost_equal<double>(r_orig.value.imag(),r_reco.value.imag()));

            ASSERT_TRUE(lazy::almost_equal<double>(k_orig.a_frac,k_reco.a_frac));
            ASSERT_TRUE(lazy::almost_equal<double>(k_orig.b_frac,k_reco.b_frac));
            ASSERT_TRUE(lazy::almost_equal<double>(k_orig.value.real(),k_reco.value.real()));
            ASSERT_TRUE(lazy::almost_equal<double>(k_orig.value.imag(),k_reco.value.imag()));
        }
    }
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
