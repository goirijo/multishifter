#include "../../autotools.hh"
#include <multishift/shifter.hpp>

#include <gtest/gtest.h>
#include <memory>

using namespace mush;

class ShifterSimpleCounting : public testing::Test
{
protected:
    std::unique_ptr<Shifter> shifter_ptr;
    int a_max = 6;
    int b_max = 6;
    virtual void SetUp() override
    {
        cu::xtal::Structure slab = cu::xtal::Structure::from_poscar(autotools::input_filesdir / "mg_stack.vasp");
        shifter_ptr.reset(new Shifter(slab, a_max, b_max));
    }
};

TEST_F(ShifterSimpleCounting, CountCategories)
{
    EXPECT_EQ(shifter_ptr->shifted_structures.size(), a_max * b_max);
    EXPECT_EQ(shifter_ptr->shift_records.size(), a_max * b_max);
    EXPECT_EQ(shifter_ptr->equivalence_map.size(), a_max * b_max);
    EXPECT_EQ(shifter_ptr->wigner_seitz_shifted_structures.size(), a_max * b_max);
}

TEST_F(ShifterSimpleCounting, WignerSeitzSanity)
{
    assert(shifter_ptr->equivalence_map == categorize_equivalently_shifted_structures(shifter_ptr->wigner_seitz_shifted_structures));
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
