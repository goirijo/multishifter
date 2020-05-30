#include "casmutils/xtal/structure.hpp"
#include "gtest/gtest.h"
#include <multishift/twist.hpp>
#include <multishift/slab.hpp>
#include <multishift/definitions.hpp>
#include <casmutils/xtal/structure_tools.hpp>

TEST(CategoryTest, SpecificTest)
{
    ASSERT_EQ(0, 0);
}

TEST(OtherCategoryTest, SpecificTest)
{
    ASSERT_EQ(0, 0);
}

/* TEST(FakeCategory, RealignLattice) */
/* { */
/*     auto slab=mush::cu::xtal::Structure::from_poscar("/home/mesto/programming/multishifter/tests/projects/Mg-mush/slab.vasp"); */
/*     auto aligned_lat=mush::make_aligned_lattice(slab.lattice()); */
/*     slab.set_lattice(aligned_lat,mush::cu::xtal::FRAC); */
/*     casmutils::xtal::write_poscar(slab,"/home/mesto/programming/multishifter/tests/projects/Mg-mush/aligned_slab.vasp"); */
/* } */

/* TEST(CategoryTest, PrintMessage) */
/* { */
/*     EXPECT_EQ(0, 0) << "I didn't actually expect this to be equal."; */
/* } */

/* TEST(CategoryTest, ExpectedFail) */
/* { */
/*     EXPECT_EQ(1, 0); */
/* } */

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
