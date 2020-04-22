#include "../../autotools.hh"
#include "casmutils/xtal/structure.hpp"
#include <gtest/gtest.h>
#include <memory>
#include <multishift/shift.hpp>
#include <vector>

namespace cu = casmutils;

class CleavingTest : public testing::Test
{
protected:
    using Structure=cu::xtal::Structure;
    void SetUp() override
    {
        b2_ptr.reset(new Structure(Structure::from_poscar(mush::autotools::input_filesdir / "b2.vasp")));
        cleaved_structures=mush::make_cleaved_structures(*b2_ptr,cleavage_values);
    }

    std::unique_ptr<Structure> b2_ptr;
    std::vector<double> cleavage_values{-0.1,0,5};
    std::vector<Structure> cleaved_structures;
};

TEST_F(CleavingTest, Consistent_AB_Vectors)
{
    for(const Structure& cleaved_struc : cleaved_structures)
    {
        EXPECT_EQ(b2_ptr->lattice().a(),cleaved_struc.lattice().a());
        EXPECT_EQ(b2_ptr->lattice().b(),cleaved_struc.lattice().b());
    }
}

TEST_F(CleavingTest, C_VectorLengths)
{
    for(int i=0; i<cleavage_values.size(); ++i)
    {
        if(cleavage_values[i]<0)
        {
            EXPECT_TRUE(b2_ptr->lattice().c().norm()>cleaved_structures[i].lattice().c().norm());
        }
        else if(cleavage_values[i]==0)
        {
            EXPECT_EQ(b2_ptr->lattice().c().norm(),cleaved_structures[i].lattice().c().norm());
        }
        else
        {
            EXPECT_TRUE(b2_ptr->lattice().c().norm()<cleaved_structures[i].lattice().c().norm());
        }
    }
}

TEST_F(CleavingTest, SitesHaventMoved)
{
    for(const Structure& cleaved_struc : cleaved_structures)
    {
        EXPECT_TRUE(cleaved_struc.basis_sites().size()==b2_ptr->basis_sites().size());
        for(int i=0; i<cleaved_struc.basis_sites().size(); ++i)
        {
            const cu::xtal::Site& cl_b=cleaved_struc.basis_sites()[i];
            const cu::xtal::Site& b=b2_ptr->basis_sites()[i];
            EXPECT_EQ(cl_b.cart(),b.cart());
            EXPECT_EQ(cl_b.label(),b.label());
        }
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
