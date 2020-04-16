#include "casmutils/xtal/structure_tools.hpp"
#include <algorithm>
#include <casmutils/misc.hpp>
#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/xtal/site.hpp>
#include <casmutils/xtal/structure.hpp>
#include <gtest/gtest.h>
#include <memory>
#include <multishift/slab.hpp>
#include <tuple>

namespace cu = casmutils;

class SlabTest : public testing::Test
{
protected:
    static cu::xtal::Structure make_hcp_structure()
    {
        Eigen::Matrix3d hcp_col_lat_mat;
        // clang-format off
        // pasted from a vasp file
        hcp_col_lat_mat << 1.5960945300000000, 2.7645168299999998, 0.0000000000000000,
                          -1.5960945300000000, 2.7645168299999998, 0.0000000000000000,
                           0.0000000000000000, 0.0000000000000000, 5.1840195700000002;
        // clang-format on

        cu::xtal::Lattice hcp_lat(hcp_col_lat_mat.transpose());
        cu::xtal::Site site0(cu::xtal::Coordinate::from_fractional(2.0 / 3, 2.0 / 3, 3.0 / 4, hcp_lat), "A");
        cu::xtal::Site site1(cu::xtal::Coordinate::from_fractional(1.0 / 3, 1.0 / 3, 1.0 / 4, hcp_lat), "A");

        cu::xtal::Structure hcp_struc(hcp_lat, std::vector<cu::xtal::Site>{site0, site1});
        return cu::xtal::make_niggli(hcp_struc);
    }

    static cu::xtal::Structure make_fcc_structure()
    {
        Eigen::Matrix3d fcc_col_lat_mat;
        // clang-format off
        // pasted from a vasp file
        fcc_col_lat_mat << 0.000000000000, 1.754750223661, 1.754750223661,
                           1.754750223661, 0.000000000000, 1.754750223661,
                           1.754750223661, 1.754750223661, 0.000000000000;
        // clang-format on

        cu::xtal::Lattice fcc_lat(fcc_col_lat_mat.transpose());
        cu::xtal::Site site0(cu::xtal::Coordinate(0, 0, 0), "A");

        return cu::xtal::Structure(fcc_lat, std::vector<cu::xtal::Site>{site0});
    }

    static cu::xtal::Structure make_b2_structure()
    {
        Eigen::Matrix3d b2_col_lat_mat;
        // clang-format off
        // pasted from a vasp file
        b2_col_lat_mat << 2.878939000000, 0.000000000000, 0.000000000000,
                          0.000000000000, 2.878939000000, 0.000000000000,
                          0.000000000000, 0.000000000000, 2.878939000000;
        // clang-format on

        cu::xtal::Lattice b2_lat(b2_col_lat_mat.transpose());
        cu::xtal::Site site0(cu::xtal::Coordinate::from_fractional(0, 0, 0, b2_lat), "A");
        cu::xtal::Site site1(cu::xtal::Coordinate::from_fractional(0.5, 0.5, 0.5, b2_lat), "B");

        return cu::xtal::Structure(b2_lat, std::vector<cu::xtal::Site>{site0, site1});
    }

    void SetUp() override
    {
        hcp_ptr.reset(new cu::xtal::Structure(SlabTest::make_hcp_structure()));
        fcc_ptr.reset(new cu::xtal::Structure(SlabTest::make_fcc_structure()));
        b2_ptr.reset(new cu::xtal::Structure(SlabTest::make_b2_structure()));

        millers_001 << 0, 0, 1;
        millers_010 << 0, 1, 0;
        millers_101 << 1, 0, 1;
        millers_211 << 2, 1, 1;
        millers_111 << 1, 1, 1;
        millers_11bar1 << 1, -1, 1;
    }

    auto make_sliced_structures(const Eigen::Vector3i& miller_indexes)
    {
        return std::make_tuple(cu::xtal::make_sliced_structure(*hcp_ptr, miller_indexes),
                               cu::xtal::make_sliced_structure(*fcc_ptr, miller_indexes),
                               cu::xtal::make_sliced_structure(*b2_ptr, miller_indexes));
    }

    auto make_sliced_lattices(const Eigen::Vector3i& miller_indexes)
    {
        return std::make_tuple(cu::xtal::make_sliced_lattice(hcp_ptr->lattice(), miller_indexes),
                               cu::xtal::make_sliced_lattice(fcc_ptr->lattice(), miller_indexes),
                               cu::xtal::make_sliced_lattice(b2_ptr->lattice(), miller_indexes));
    }

    std::unique_ptr<cu::xtal::Structure> hcp_ptr;
    std::unique_ptr<cu::xtal::Structure> fcc_ptr;
    std::unique_ptr<cu::xtal::Structure> b2_ptr;

    Eigen::Vector3i millers_001;
    Eigen::Vector3i millers_010;
    Eigen::Vector3i millers_101;
    Eigen::Vector3i millers_211;
    Eigen::Vector3i millers_111;
    Eigen::Vector3i millers_11bar1;
};

TEST_F(SlabTest, SliceLatticeConsistency)
{
    auto miller_set = {millers_001, millers_010, millers_101, millers_211, millers_111, millers_11bar1};
    for (const Eigen::Vector3i& millers : miller_set)
    {
        auto [hcp_lat, fcc_lat, b2_lat] = make_sliced_lattices(millers);
        auto [hcp_slice, fcc_slice, b2_slice] = make_sliced_structures(millers);

        EXPECT_TRUE(cu::is_equal<cu::xtal::LatticeEquals_f>(hcp_lat, hcp_slice.lattice(), 1e-5));
        EXPECT_TRUE(cu::is_equal<cu::xtal::LatticeEquals_f>(fcc_lat, fcc_slice.lattice(), 1e-5));
        EXPECT_TRUE(cu::is_equal<cu::xtal::LatticeEquals_f>(b2_lat, b2_slice.lattice(), 1e-5));
    }
}

TEST_F(SlabTest, LatticeSlice001)
{
    auto [hcp_lat, fcc_lat, b2_lat] = make_sliced_lattices(millers_001);
    EXPECT_TRUE(cu::is_equal<cu::xtal::LatticeEquals_f>(hcp_lat, hcp_ptr->lattice(), 1e-5));
    EXPECT_TRUE(cu::is_equal<cu::xtal::LatticeEquals_f>(fcc_lat, fcc_ptr->lattice(), 1e-5));
    EXPECT_TRUE(cu::is_equal<cu::xtal::LatticeEquals_f>(b2_lat, b2_ptr->lattice(), 1e-5));
}

#include <casm/crystallography/Superlattice.hh>

TEST_F(SlabTest, LatticeSlice_fcc_111)
{
    cu::xtal::Lattice fcc_slice_lat = cu::xtal::make_sliced_lattice(fcc_ptr->lattice(), millers_111);
    Eigen::Matrix3i prim_to_3layer_fcc_mat;
    prim_to_3layer_fcc_mat << -1,0,1,1,-1,1,0,1,1;

    cu::xtal::Structure sliced_struc = cu::xtal::make_sliced_structure(*fcc_ptr, millers_111);
    cu::xtal::print_poscar(sliced_struc, std::cout);
    cu::xtal::print_poscar(*this->fcc_ptr, std::cout);

    auto tmp=CASM::xtal::Superlattice::smooth_prim(fcc_ptr->lattice().__get(), sliced_struc.lattice().__get());
    std::cout<<tmp.transformation_matrix()<<std::endl;

    cu::xtal::Lattice fcc_3layer = cu::xtal::make_superlattice(fcc_ptr->lattice(), prim_to_3layer_fcc_mat);
    /* std::cout<<fcc_slice_lat.column_vector_matrix()<<std::endl<<std::endl<<conventional_fcc.column_vector_matrix()<<std::endl<<std::endl;
     */
    EXPECT_TRUE(cu::is_equal<cu::xtal::LatticeEquals_f>(fcc_slice_lat, fcc_3layer, 1e-5));
}

TEST_F(SlabTest, StructureSlice_b2_101)
{
    cu::xtal::Structure b2_slice = cu::xtal::make_sliced_structure(*b2_ptr, millers_101);

    EXPECT_EQ(b2_slice.basis_sites().size(), 4);

    cu::xtal::Site site0(cu::xtal::Coordinate::from_fractional(0, 0, 0, b2_slice.lattice()), "A");
    cu::xtal::Site site1(cu::xtal::Coordinate::from_fractional(0.5, 0.5, 0, b2_slice.lattice()), "B");

    cu::xtal::SiteEquals_f is_equal_site0(site0, 1e-5);
    EXPECT_TRUE(std::find_if(b2_slice.basis_sites().begin(), b2_slice.basis_sites().end(), is_equal_site0) != b2_slice.basis_sites().end());

    cu::xtal::SiteEquals_f is_equal_site1(site1, 1e-5);
    EXPECT_TRUE(std::find_if(b2_slice.basis_sites().begin(), b2_slice.basis_sites().end(), is_equal_site1) != b2_slice.basis_sites().end());
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
