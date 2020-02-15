#include "gtest/gtest.h"

TEST(CategoryTest, SpecificTest)
{
    ASSERT_EQ(0, 0);
}

TEST(OtherCategoryTest, SpecificTest)
{
    ASSERT_EQ(0, 0);
}

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
