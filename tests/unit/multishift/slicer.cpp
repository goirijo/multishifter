#include "../../autotools.hh"
#include "casmutils/xtal/structure.hpp"
#include <filesystem>
#include <fstream>
#include <multishift/slicer.hpp>

#include <gtest/gtest.h>
#include <memory>

using namespace mush;

class SlicerSimpleCounting : public testing::Test
{
    protected:
        std::unique_ptr<Slicer> slicer_ptr;
        virtual void SetUp() override
        {
            cu::xtal::Structure prim=cu::xtal::Structure::from_poscar(autotools::input_filesdir/"b2.vasp");
            Eigen::Vector3i miller_indexes(1,1,1);

            slicer_ptr.reset(new Slicer(prim,miller_indexes));
        }
};

TEST_F(SlicerSimpleCounting, SlicedIsBigger)
{
    EXPECT_TRUE(slicer_ptr->prim.basis_sites().size()<slicer_ptr->sliced_prim.basis_sites().size());
}

TEST_F(SlicerSimpleCounting, MultipleFloors)
{
    EXPECT_TRUE(slicer_ptr->floored_sliced_prims.size()==2);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
