#include "../../autotools.hh"
#include <cmath>
#include <gtest/gtest.h>

#include <multishift/definitions.hpp>
#include <multishift/twist.hpp>
#include <vector>

namespace casmutils::xtal
{
xtal::Lattice make_sliced_lattice(const xtal::Lattice& unit_lattice, const Eigen::Vector3i& miller_indexes);
}

using namespace mush;

class TwistTest : public testing::Test
{
protected:
    void SetUp()
    {
        Eigen::Matrix3d col_lat_mat;
        col_lat_mat << 3, 0, 0, 0, 4, 0, 0, 0, 8;
        Lattice ortho_lat(col_lat_mat.col(0), col_lat_mat.col(1), col_lat_mat.col(2));

        sliced_lattices.emplace_back(cu::xtal::make_sliced_lattice(ortho_lat, Eigen::Vector3i(0, 0, 1)));
        return;

        for (int i : {-1, 0, 2})
        {
            for (int j : {-2, 1, 3})
            {
                for (int k : {1, -4, 2})
                {
                    Eigen::Vector3i millers(i, j, k);
                    sliced_lattices.emplace_back(cu::xtal::make_sliced_lattice(ortho_lat, millers));
                }
            }
        }
        return;
    }

    using Lattice = cu::xtal::Lattice;
    std::vector<Lattice> sliced_lattices;
    double tol = 1e-10;

    //TODO: Move this somewhere useful
    static double degrees_between_vectors(const Eigen::Vector3d& a1, const Eigen::Vector3d& a2)
    {
        return 2 * std::atan2((a2.norm() * a1 - a1.norm() * a2).norm(), (a2.norm() * a1 + a1.norm() * a2).norm()) * 180 / M_PI;
    }

    static double signed_degrees_between_vectors(const Eigen::Vector3d& a1, const Eigen::Vector3d& a2, const Eigen::Vector3d& normal)
    {
        int sign = 1;
        if (normal.dot(a1.cross(a2)) < 0)
        {
            sign = -1;
        }
        return sign * degrees_between_vectors(a1, a2);
    }
};

TEST_F(TwistTest, TwistMatrixIsUnitary)
{
    for (const Lattice& start_lat : sliced_lattices)
    {
        for (double degrees : {-5, 0, 30})
        {

            Eigen::Matrix3d R = make_twist_rotation_matrix(start_lat, degrees);
            EXPECT_TRUE(almost_equal(R.determinant(), 1.0, tol));

            // Unitary if R.T*R is identity
            Eigen::Matrix3d RtransposeR = R.transpose() * R;
            Eigen::Matrix3d I(Eigen::Matrix3d::Identity());

            EXPECT_TRUE(almost_equal(RtransposeR, I, tol));
        }
    }
}

TEST_F(TwistTest, TwistedLatticeAngles)
{
    for (const Lattice& start_lat : sliced_lattices)
    {
        Eigen::Vector3d normal = start_lat.a().cross(start_lat.b());
        for (double angle : {-33.0, -15.0, -1.0, -0.1, 0.0, 0.00000001, 0.2, 30.0, 90.0, 180.0})
        {
            Lattice rotated_lat = make_twisted_lattice(start_lat, angle);

            EXPECT_TRUE(almost_equal(normal, rotated_lat.a().cross(rotated_lat.b()),tol));
            EXPECT_TRUE(almost_equal(start_lat.c().dot(normal), rotated_lat.c().dot(normal), tol));

            double a_deg = signed_degrees_between_vectors(start_lat.a(), rotated_lat.a(), normal);
            double b_deg = signed_degrees_between_vectors(start_lat.b(), rotated_lat.b(), normal);

            EXPECT_TRUE(almost_equal(a_deg, b_deg, tol)) << a_deg << " vs " << b_deg;
            EXPECT_TRUE(almost_equal(a_deg, angle, tol)) << a_deg << " vs " << angle;
        }
    }
}

TEST_F(TwistTest, ObviousTwistMatrix)
{
    Eigen::Matrix3d ortho_mat;
    ortho_mat << 4, 0, 0, 0, 6, 0, 0, 0, 9;

    cu::xtal::Lattice ortho_lat(ortho_mat);

    Eigen::Matrix3d rotate_0 = make_twist_rotation_matrix(ortho_lat, 0);
    Eigen::Matrix3d I(Eigen::Matrix3d::Identity());
    EXPECT_TRUE(almost_equal(rotate_0, I, tol));

    Eigen::Matrix3d rotate_90 = make_twist_rotation_matrix(ortho_lat, 90);
    Eigen::Matrix3d exptected_rotate_90;
    exptected_rotate_90 << 0, -1, 0, 1, 0, 0, 0, 0, 1;
    EXPECT_TRUE(almost_equal(rotate_90, exptected_rotate_90, 1e-10));
}

TEST_F(TwistTest, LessObviousTwistMatrix)
{
    Eigen::Matrix3d ortho_mat;
    ortho_mat << 0, 0, 9, 4, 0, 0, 0, 6, 0;

    cu::xtal::Lattice ortho_lat(ortho_mat);

    Eigen::Matrix3d rotate_0 = make_twist_rotation_matrix(ortho_lat, 0);
    Eigen::Matrix3d I(Eigen::Matrix3d::Identity());
    EXPECT_TRUE(almost_equal(rotate_0, I, 1e-10));

    Eigen::Matrix3d rotate_90 = make_twist_rotation_matrix(ortho_lat, 90);
    Eigen::Matrix3d exptected_rotate_90;
    exptected_rotate_90 << 1, 0, 0, 0, 0, -1, 0, 1, 0;

    EXPECT_TRUE(almost_equal(rotate_90, exptected_rotate_90, tol));
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
