#include "../../autotools.hh"
#include <algorithm>
#include <casmutils/misc.hpp>
#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/xtal/site.hpp>
#include <casmutils/xtal/structure.hpp>
#include <casmutils/xtal/structure_tools.hpp>
#include <gtest/gtest.h>
#include <memory>
#include <multishift/slab.hpp>
#include <tuple>

namespace cu = casmutils;

class SlicingTest : public testing::Test
{
protected:
    using Structure=cu::xtal::Structure;
    void SetUp() override
    {
        hcp_ptr.reset(new Structure(Structure::from_poscar(mush::autotools::input_filesdir / "hcp.vasp")));
        fcc_ptr.reset(new Structure(Structure::from_poscar(mush::autotools::input_filesdir / "fcc.vasp")));
        b2_ptr.reset(new Structure(Structure::from_poscar(mush::autotools::input_filesdir / "b2.vasp")));

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

    std::unique_ptr<Structure> hcp_ptr;
    std::unique_ptr<Structure> fcc_ptr;
    std::unique_ptr<Structure> b2_ptr;

    Eigen::Vector3i millers_001;
    Eigen::Vector3i millers_010;
    Eigen::Vector3i millers_101;
    Eigen::Vector3i millers_211;
    Eigen::Vector3i millers_111;
    Eigen::Vector3i millers_11bar1;

    double tol = 1e-5;
};

TEST_F(SlicingTest, SliceLatticeConsistency)
{
    auto miller_set = {millers_001, millers_010, millers_101, millers_211, millers_111, millers_11bar1};
    for (const Eigen::Vector3i& millers : miller_set)
    {
        auto [hcp_lat, fcc_lat, b2_lat] = make_sliced_lattices(millers);
        auto [hcp_slice, fcc_slice, b2_slice] = make_sliced_structures(millers);

        EXPECT_TRUE(cu::is_equal<cu::xtal::LatticeEquals_f>(hcp_lat, hcp_slice.lattice(), tol));
        EXPECT_TRUE(cu::is_equal<cu::xtal::LatticeEquals_f>(fcc_lat, fcc_slice.lattice(), tol));
        EXPECT_TRUE(cu::is_equal<cu::xtal::LatticeEquals_f>(b2_lat, b2_slice.lattice(), tol));
    }
}

TEST_F(SlicingTest, LatticeSlice001)
{
    auto [hcp_lat, fcc_lat, b2_lat] = make_sliced_lattices(millers_001);
    EXPECT_TRUE(cu::is_equal<cu::xtal::LatticeEquals_f>(hcp_lat, hcp_ptr->lattice(), tol));
    EXPECT_TRUE(cu::is_equal<cu::xtal::LatticeEquals_f>(fcc_lat, fcc_ptr->lattice(), tol));
    EXPECT_TRUE(cu::is_equal<cu::xtal::LatticeEquals_f>(b2_lat, b2_ptr->lattice(), tol));
}

TEST_F(SlicingTest, LatticeSlice_fcc_111)
{
    cu::xtal::Lattice fcc_slice_lat = cu::xtal::make_sliced_lattice(fcc_ptr->lattice(), millers_111);
    Eigen::Matrix3i prim_to_3layer_fcc_mat;
    prim_to_3layer_fcc_mat << -1, 0, 1, 1, -1, 1, 0, 1, 1;

    Structure sliced_struc = cu::xtal::make_sliced_structure(*fcc_ptr, millers_111);

    cu::xtal::Lattice fcc_3layer = cu::xtal::make_superlattice(fcc_ptr->lattice(), prim_to_3layer_fcc_mat);
    EXPECT_TRUE(cu::is_equal<cu::xtal::LatticeEquals_f>(fcc_slice_lat, fcc_3layer, tol));
}

TEST_F(SlicingTest, StructureSlice_b2_101)
{
    Structure b2_slice = cu::xtal::make_sliced_structure(*b2_ptr, millers_101);

    EXPECT_EQ(b2_slice.basis_sites().size(), 4);

    cu::xtal::Site site0(cu::xtal::Coordinate::from_fractional(0, 0, 0, b2_slice.lattice()), "A");
    cu::xtal::Site site1(cu::xtal::Coordinate::from_fractional(0.5, 0.5, 0, b2_slice.lattice()), "B");

    cu::xtal::SiteEquals_f is_equal_site0(site0, tol);
    EXPECT_TRUE(std::find_if(b2_slice.basis_sites().begin(), b2_slice.basis_sites().end(), is_equal_site0) != b2_slice.basis_sites().end());

    cu::xtal::SiteEquals_f is_equal_site1(site1, tol);
    EXPECT_TRUE(std::find_if(b2_slice.basis_sites().begin(), b2_slice.basis_sites().end(), is_equal_site1) != b2_slice.basis_sites().end());
}

//******************************************************************************//

class SlabTest : public testing::Test
{
protected:
    using Structure=cu::xtal::Structure;
    void SetUp() override
    {
        b2_101_ptr.reset(new Structure(Structure::from_poscar(mush::autotools::input_filesdir / "b2_101.vasp")));
        b2_101_stack5_ptr.reset(new Structure(mush::make_stacked_slab(*b2_101_ptr, 5)));

        cu::xtal::print_poscar(*b2_101_ptr, std::cout);
        cu::xtal::print_poscar(*b2_101_stack5_ptr, std::cout);
    }

    std::unique_ptr<Structure> b2_101_ptr;
    std::unique_ptr<Structure> b2_101_stack5_ptr;
};

TEST_F(SlabTest, Consistent_AB_Vectors)
{

    EXPECT_EQ(b2_101_ptr->lattice().a(), b2_101_stack5_ptr->lattice().a());
    EXPECT_EQ(b2_101_ptr->lattice().b(), b2_101_stack5_ptr->lattice().b());
}

TEST_F(SlabTest, Consistent_C_Vector) { EXPECT_EQ(b2_101_ptr->lattice().c() * 5, b2_101_stack5_ptr->lattice().c()); }

TEST_F(SlabTest, FlooredStructure)
{
    auto coord_is_origin = [](const cu::xtal::Site& s) { return s.cart().isZero(); };
    const auto& b2_101_basis = this->b2_101_ptr->basis_sites();
    int ix = 0;
    for (; ix < b2_101_basis.size(); ++ix)
    {
        if (coord_is_origin(b2_101_basis[ix]))
        {
            break;
        }
    }

    std::string original_origin_label = b2_101_basis[ix].label();
    int floor_ix = 0;
    for (; floor_ix < b2_101_basis.size(); ++floor_ix)
    {
        if (b2_101_basis[floor_ix].label() != original_origin_label)
        {
            break;
        }
    }

    Structure floored_structure = mush::make_floored_structure(*(this->b2_101_ptr), floor_ix);
    cu::xtal::print_poscar(floored_structure,std::cout);
    const auto& floored_basis=floored_structure.basis_sites();
    auto origin_site_it = std::find_if(floored_basis.begin(), floored_basis.end(), coord_is_origin);
    EXPECT_TRUE(origin_site_it->label() != original_origin_label);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
