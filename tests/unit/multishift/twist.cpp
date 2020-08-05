#include "../../autotools.hh"
#include "casmutils/xtal/lattice.hpp"
#include "casmutils/xtal/site.hpp"
#include "casmutils/xtal/structure.hpp"
#include "casmutils/xtal/structure_tools.hpp"
#include "multishift/slab.hpp"
#include <cmath>
#include <fstream>
#include <gtest/gtest.h>

#include <memory>
#include <multishift/definitions.hpp>
#include <multishift/twist.hpp>
#include <string>
#include <tuple>
#include <vector>

using namespace mush;

class TwistTest : public testing::Test
{
protected:
    void SetUp()
    {
        Eigen::Matrix3d col_lat_mat;
        col_lat_mat << 3, 0, 0, 0, 4, 0, 0, 0, 8;
        Lattice ortho_lat(col_lat_mat.col(0), col_lat_mat.col(1), col_lat_mat.col(2));

        for (int i : {-1, 0, 2})
        {
            for (int j : {-2, 1, 3})
            {
                for (int k : {1, -4, 2})
                {
                    Eigen::Vector3i millers(i, j, k);
                    sliced_lattices.emplace_back(cu::xtal::slice_along_plane(ortho_lat, millers));
                }
            }
        }
        return;
    }

    std::vector<Lattice> sliced_lattices;
    double tol = 1e-10;

    // TODO: Move this somewhere useful
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

            EXPECT_TRUE(almost_equal(normal, rotated_lat.a().cross(rotated_lat.b()), tol));
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

TEST_F(TwistTest, AlginLattice)
{
    for (const Lattice& start_lat : sliced_lattices)
    {
        Lattice aligned_lat = make_aligned_lattice(start_lat);
        EXPECT_TRUE(almost_equal(aligned_lat.a()(1), 0.0, tol));
        EXPECT_TRUE(almost_equal(aligned_lat.a()(2), 0.0, tol));
        EXPECT_TRUE(almost_equal(aligned_lat.b()(2), 0.0, tol));

        EXPECT_TRUE(almost_equal(start_lat.a().norm(), aligned_lat.a().norm(), tol));
        EXPECT_TRUE(almost_equal(start_lat.b().norm(), aligned_lat.b().norm(), tol));
        EXPECT_TRUE(almost_equal(start_lat.c().norm(), aligned_lat.c().norm(), tol));
    }
}

TEST_F(TwistTest, ApproximantMoireLattice)
{
    for (const Lattice& start_lat : sliced_lattices)
    {
        for (double angle : {-3.00, -2.0, -1.0, -0.50, 0.50, 2.00, 5.00, 22.0})
        {
            using LATTICE=MoireApproximant::LATTICE;
            MoireLattice moire(start_lat, angle);
            MoireApproximant approx_moire(moire.aligned_moire_lattice,moire.aligned_lattice,moire.rotated_lattice);
            const auto& moire_lat=approx_moire.approximate_moire_lattice;
            const auto& aligned_lat=approx_moire.approximate_lattices.at(LATTICE::ALIGNED);
            const auto& rot_lat=approx_moire.approximate_lattices.at(LATTICE::ROTATED);

            Eigen::Matrix3d aligned_to_moire_transform = aligned_lat.column_vector_matrix().inverse() * moire_lat.column_vector_matrix();
            Eigen::Matrix3d aligned_to_moire_transform_diff =
                aligned_to_moire_transform - aligned_to_moire_transform.unaryExpr([](double x) { return std::round(x); });

            Eigen::Matrix3d rot_to_moire_transform = rot_lat.column_vector_matrix().inverse() * moire_lat.column_vector_matrix();
            Eigen::Matrix3d rot_to_moire_transform_diff =
                rot_to_moire_transform - rot_to_moire_transform.unaryExpr([](double x) { return std::round(x); });

            EXPECT_TRUE(almost_zero(aligned_to_moire_transform_diff.block<2, 2>(0, 0), 1e-8));
            EXPECT_TRUE(almost_zero(rot_to_moire_transform_diff.block<2, 2>(0, 0), 1e-8));
        }
    }
}

TEST_F(TwistTest, MoireApproximantStrainMatch)
{
    for (const Lattice& start_lat : sliced_lattices)
    {
        for (double angle : {-3.00, -2.0, -1.0, -0.50, 0.50, 2.00, 5.00, 22.0})
        {
            using LATTICE=MoireApproximant::LATTICE;
            MoireLattice moire(start_lat, angle);
            const Lattice& moire_lat = moire.aligned_moire_lattice;
            const Lattice& aligned_lat = moire.aligned_lattice;
            const Lattice& rotated_lat = moire.rotated_lattice;

            MoireApproximant approx_moire(moire_lat, aligned_lat, rotated_lat);

            for(LATTICE lat : {LATTICE::ALIGNED,LATTICE::ROTATED})
            {
                const Eigen::Matrix3d& L=moire.real(lat).column_vector_matrix();
                const Eigen::Matrix3d& L_prime=approx_moire.approximate_lattices.at(lat).column_vector_matrix();

                Eigen::Matrix3d F=L_prime*L.inverse();
                EXPECT_TRUE(almost_equal(F,approx_moire.approximation_deformations.at(lat)));
            }
        }
    }
}

TEST_F(TwistTest, MoireApproximantVectorMatch)
{
    for (const Lattice& start_lat : sliced_lattices)
    {
        for (double angle : {-3.00, -2.0, -1.0, -0.50, 0.50, 2.00, 5.00, 22.0})
        {
            using LATTICE=MoireApproximant::LATTICE;
            MoireLattice moire(start_lat, angle);
            const Lattice& moire_lat = moire.aligned_moire_lattice;
            const Lattice& aligned_lat = moire.aligned_lattice;
            const Lattice& rotated_lat = moire.rotated_lattice;

            MoireApproximant approx_moire(moire_lat, aligned_lat, rotated_lat);

            Lattice aligned_super =
                cu::xtal::make_superlattice(approx_moire.approximate_lattices.at(LATTICE::ALIGNED),
                                            approx_moire.approximate_moire_integer_transformations.at(LATTICE::ALIGNED).cast<int>());

            Lattice rot_super =
                cu::xtal::make_superlattice(approx_moire.approximate_lattices.at(LATTICE::ROTATED),
                                            approx_moire.approximate_moire_integer_transformations.at(LATTICE::ROTATED).cast<int>());

            EXPECT_TRUE(almost_equal(aligned_super.a(), rot_super.a()));
            EXPECT_TRUE(almost_equal(aligned_super.b(), rot_super.b()));
        }
    }
}

TEST_F(TwistTest, MoireGeneratorVectorMatch)
{
    for (const Lattice& start_lat : sliced_lattices)
    {
        for (double angle : {-3.00, -2.0, -1.0, -0.50, 0.50, 2.00, 5.00, 22.0})
        {
            MoireGenerator approx_moire(start_lat, angle);
            using ZONE = MoireGenerator::ZONE;
            using LATTICE = MoireGenerator::LATTICE;

            Lattice aligned_super =
                cu::xtal::make_superlattice(approx_moire.approximate_lattice(ZONE::ALIGNED,LATTICE::ALIGNED),
                                            approx_moire.approximate_moire_integer_transformation(ZONE::ALIGNED,LATTICE::ALIGNED).cast<int>());

            Lattice rot_super =
                cu::xtal::make_superlattice(approx_moire.approximate_lattice(ZONE::ALIGNED,LATTICE::ROTATED),
                                            approx_moire.approximate_moire_integer_transformation(ZONE::ALIGNED,LATTICE::ROTATED).cast<int>());

            EXPECT_TRUE(almost_equal(aligned_super.a(), rot_super.a()));
            EXPECT_TRUE(almost_equal(aligned_super.b(), rot_super.b()));
        }
    }
}

TEST_F(TwistTest, MoireApproximantAlignedRotatedMatch)
{
    for (const Lattice& start_lat : sliced_lattices)
    {
        for (double angle : {-3.00, -2.0, -1.0, -0.50, 0.50, 2.00, 5.00, 22.0})
        {
            // TODO: Make sure that the aligned and rotated lattices
            //(both perfect and approximated) are exclusively related
            // by a rotation
        }
    }
}

TEST_F(TwistTest, ApproximantPrismaticMoireDeformation)
{
    for (const Lattice& start_lat : sliced_lattices)
    {
        for (double angle : {-3.00, -2.0, -1.0, -0.50, 0.50, 2.00, 5.00, 22.0})
        {
            using LATTICE=MoireLattice::LATTICE;
            MoireLattice moire(start_lat, angle);
            MoireApproximant approx_moire(moire.aligned_moire_lattice,moire.aligned_lattice,moire.rotated_lattice);
            const auto& moire_lat=approx_moire.approximate_moire_lattice;
            const auto& aligned_lat=approx_moire.approximate_lattices.at(LATTICE::ALIGNED);
            const auto& rot_lat=approx_moire.approximate_lattices.at(LATTICE::ROTATED);

            Eigen::Matrix3d aligned_to_moire_transform = aligned_lat.column_vector_matrix().inverse() * moire_lat.column_vector_matrix();
            Eigen::Matrix3d aligned_to_moire_transform_diff =
                aligned_to_moire_transform - aligned_to_moire_transform.unaryExpr([](double x) { return std::round(x); });

            Eigen::Matrix3d rot_to_moire_transform = rot_lat.column_vector_matrix().inverse() * moire_lat.column_vector_matrix();
            Eigen::Matrix3d rot_to_moire_transform_diff =
                rot_to_moire_transform - rot_to_moire_transform.unaryExpr([](double x) { return std::round(x); });

            EXPECT_TRUE(almost_zero(aligned_to_moire_transform_diff.block<2, 2>(0, 0), 1e-8));
            EXPECT_TRUE(almost_zero(rot_to_moire_transform_diff.block<2, 2>(0, 0), 1e-8));
        }
    }
}

TEST_F(TwistTest, PrismaticLattice)
{
    for (const Lattice& start_lat : sliced_lattices)
    {
        for (double angle : {-3.00, -2.0, -1.0, -0.50, 0.50, 2.00, 5.00, 22.0})
        {
            Lattice prismatic_lat = make_prismatic_lattice(start_lat);
            EXPECT_EQ(start_lat.a(), prismatic_lat.a());
            EXPECT_EQ(start_lat.b(), prismatic_lat.b());

            EXPECT_TRUE(almost_equal(prismatic_lat.a().dot(prismatic_lat.c()), 0.0));
            EXPECT_TRUE(almost_equal(prismatic_lat.b().dot(prismatic_lat.c()), 0.0));
        }
    }
}

TEST(HexagonalTwistTest, BrillouinOverlapMoiresRelatedByIntegerTransform)
{
    int passed = 0;
    int missed = 0;
    Eigen::Matrix3d row_lat_mat;
    row_lat_mat << 2.4684159756, 0.0000000000, 0.0000000000, -1.2342079878, 2.1377109420, 0.0000000000, 0.0000000000, 0.0000000000,
        9.9990577698;

    Lattice triangular(row_lat_mat.row(0), row_lat_mat.row(1), row_lat_mat.row(2));

    for (double degrees = 0.5; degrees < 361; degrees += 1.0)
    {
        bool full_brillouin_overlap = true;

        MoireLattice moire(triangular, degrees);
        for (auto  lat : {MoireLattice::LATTICE::ALIGNED, MoireLattice::LATTICE::ROTATED})
        {
            for (int i = 0; i < 2; ++i)
            {
                if (!moire.is_within_brillouin_zone_overlap[lat][i])
                {
                    full_brillouin_overlap = false;
                }
            }
        }

        EXPECT_TRUE(full_brillouin_overlap);

        // TODO: Expose transformation matrix stuff (Superlattice?) in cu
        const auto& M = moire.aligned_moire_lattice;
        const auto& Mt = moire.rotated_moire_lattice;
        Eigen::Matrix3d m_to_m = M.column_vector_matrix().inverse() * Mt.column_vector_matrix();

        Eigen::Matrix3i im_to_m = CASM::iround(m_to_m);
        Eigen::Matrix3d error = m_to_m - im_to_m.cast<double>();

        EXPECT_TRUE(almost_equal(error, Eigen::Matrix3d::Zero()));

        ++passed;
    }
}

class ReducedAngleGrapheneTwistTest : public testing::Test
{
protected:
    void SetUp()
    {
        Eigen::Matrix3d row_lat_mat;
        row_lat_mat << 2.4684159756, 0.0000000000, 0.0000000000, -1.2342079878, 2.1377109420, 0.0000000000, 0.0000000000, 0.0000000000,
            9.9990577698;

        graphene_lat_ptr.reset(new Lattice(row_lat_mat.row(0), row_lat_mat.row(1), row_lat_mat.row(2)));
    }

    std::unique_ptr<Lattice> graphene_lat_ptr;
};

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
