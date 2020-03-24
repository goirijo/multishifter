#include "gtest/gtest.h"
#include <casmutils/xtal/structure.hpp>
#include <multishift/slab.hpp>
#include <memory>

namespace cu=casmutils;

class SlabTest : public testing::Test
{
    void SetUp() override
    {
        return;
    }

    static cu::xtal::Structure make_hcp_structure();
    static cu::xtal::Structure make_fcc_structure();
    static cu::xtal::Structure make_b2_structure();
    static cu::xtal::Structure make_cubic_structure();

    std::unique_ptr<cu::xtal::Structure> hcp_ptr;
    std::unique_ptr<cu::xtal::Structure> fcc_ptr;
    std::unique_ptr<cu::xtal::Structure> b2_ptr;
    std::unique_ptr<cu::xtal::Structure> cubic_ptr;
};

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
