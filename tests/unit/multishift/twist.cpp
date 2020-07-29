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

namespace casmutils::xtal
{
xtal::Lattice make_sliced_lattice(const xtal::Lattice& unit_lattice, const Eigen::Vector3i& miller_indexes);

namespace frankenstein
{
xtal::Structure stack(const std::vector<xtal::Structure>& sub_strucs);
}
} // namespace casmutils::xtal

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
                    sliced_lattices.emplace_back(cu::xtal::make_sliced_lattice(ortho_lat, millers));
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

TEST_F(TwistTest, MoireLattice)
{
    std::ofstream dumpstream("./moirelatticeout.txt");
    for (const Lattice& start_lat : sliced_lattices)
    {
        for (double angle : {-3.00, -2.0, -1.0, -0.50, 0.50, 2.00, 5.00, 25.0})
        {
            auto [moire_lat, aligned_lat, rot_lat] = make_aligned_moire_lattice(start_lat, angle);

            Eigen::Matrix3d aligned_to_moire_transform = aligned_lat.column_vector_matrix().inverse() * moire_lat.column_vector_matrix();
            Eigen::Matrix3d aligned_to_moire_transform_diff =
                aligned_to_moire_transform - aligned_to_moire_transform.unaryExpr([](double x) { return std::round(x); });

            Eigen::Matrix3d rot_to_moire_transform = rot_lat.column_vector_matrix().inverse() * moire_lat.column_vector_matrix();
            Eigen::Matrix3d rot_to_moire_transform_diff =
                rot_to_moire_transform - rot_to_moire_transform.unaryExpr([](double x) { return std::round(x); });

            // TODO: How to test the moire lattice? Plotting the two lattices seems
            // to give the correct results. I THINK this check makes sense
            EXPECT_TRUE(rot_to_moire_transform_diff.isApprox(rot_to_moire_transform_diff));

            // This block is for using a dumb python script for visualization
            dumpstream << "****************************\n\n";
            dumpstream << aligned_lat.column_vector_matrix() << "\n\n";
            dumpstream << rot_lat.column_vector_matrix() << "\n\n";
            dumpstream << moire_lat.column_vector_matrix() << "\n\n";
        }
    }
    dumpstream.close();
}

TEST_F(TwistTest, ApproximantMoireLattice)
{
    for (const Lattice& start_lat : sliced_lattices)
    {
        for (double angle : {-3.00, -2.0, -1.0, -0.50, 0.50, 2.00, 5.00, 22.0})
        {
            MoireLattice moire(start_lat, angle);
            MoireApproximant approx_moire(moire.aligned_moire_lattice,moire.aligned_lattice,moire.rotated_lattice);
            const auto& moire_lat=approx_moire.approximate_moire_lattice;
            const auto& aligned_lat=approx_moire.approximate_lattices.at(&moire.aligned_lattice);
            const auto& rot_lat=approx_moire.approximate_lattices.at(&moire.rotated_lattice);

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
            MoireLattice moire(start_lat, angle);
            const Lattice& moire_lat = moire.aligned_moire_lattice;
            const Lattice& aligned_lat = moire.aligned_lattice;
            const Lattice& rotated_lat = moire.rotated_lattice;

            MoireApproximant approx_moire(moire_lat, aligned_lat, rotated_lat);

            for(const Lattice* lat : {&aligned_lat,&rotated_lat})
            {
                const Eigen::Matrix3d& L=lat->column_vector_matrix();
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
            MoireLattice moire(start_lat, angle);
            const Lattice& moire_lat = moire.aligned_moire_lattice;
            const Lattice& aligned_lat = moire.aligned_lattice;
            const Lattice& rotated_lat = moire.rotated_lattice;

            MoireApproximant approx_moire(moire_lat, aligned_lat, rotated_lat);

            Lattice aligned_super =
                cu::xtal::make_superlattice(approx_moire.approximate_lattices.at(&aligned_lat),
                                            approx_moire.approximate_moire_integer_transformations.at(&aligned_lat).cast<int>());

            Lattice rot_super =
                cu::xtal::make_superlattice(approx_moire.approximate_lattices.at(&rotated_lat),
                                            approx_moire.approximate_moire_integer_transformations.at(&rotated_lat).cast<int>());

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
            MoireLattice moire(start_lat, angle);
            MoireApproximant approx_moire(moire.aligned_moire_lattice,moire.aligned_lattice,moire.rotated_lattice);
            const auto& moire_lat=approx_moire.approximate_moire_lattice;
            const auto& aligned_lat=approx_moire.approximate_lattices.at(&moire.aligned_lattice);
            const auto& rot_lat=approx_moire.approximate_lattices.at(&moire.rotated_lattice);

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
        for (const Lattice* lat : {&moire.aligned_moire_lattice, &moire.rotated_moire_lattice})
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

namespace dirty
{
double degrees_to_radians(double angle) { return angle * M_PI / 180.0; }

double radians_to_degrees(double rad) { return rad * 108 / M_PI; }
Eigen::Matrix3d make_rotation_matrix(double angle)
{
    double rad = degrees_to_radians(angle);
    Eigen::Matrix3d rotation = Eigen::Matrix3d::Identity();

    double c = std::cos(rad);
    double s = std::sin(rad);

    rotation(0, 0) = c;
    rotation(0, 1) = -s;
    rotation(1, 0) = s;
    rotation(1, 1) = c;

    return rotation;
}

cu::xtal::Lattice make_hexagonal_lattice()
{
    Eigen::Vector3d a(0, 3, 0);
    Eigen::Matrix3d R = make_rotation_matrix(60);
    Eigen::Vector3d b = R * a;
    Eigen::Vector3d c(0, 0, 1);

    return cu::xtal::Lattice(a, b, c);
}

cu::xtal::Lattice make_rectangular_lattice()
{
    Eigen::Vector3d a(5, 0, 0);
    Eigen::Vector3d b(0, 2, 0);
    Eigen::Vector3d c(0, 0, 1);
    return cu::xtal::Lattice(a, b, c);
}

cu::xtal::Lattice make_square_lattice()
{
    Eigen::Vector3d a(4, 0, 0);
    Eigen::Vector3d b(0, 4, 0);
    Eigen::Vector3d c(0, 0, 1);
    return cu::xtal::Lattice(a, b, c);
}

cu::xtal::Lattice make_ugly_lattice()
{
    Eigen::Vector3d a(3, 0, 0);
    Eigen::Vector3d b(4, 0, 0);
    Eigen::Matrix3d R = make_rotation_matrix(20);
    b = R * b;
    Eigen::Vector3d c(0, 0, 1);
    return cu::xtal::Lattice(a, b, c);
}

cu::xtal::Lattice make_almost_ugly_lattice()
{
    Eigen::Vector3d a(4, 0, 0);
    Eigen::Vector3d b(4, 0, 0);
    Eigen::Matrix3d R = make_rotation_matrix(20);
    b = R * b;
    Eigen::Vector3d c(0, 0, 1);
    return cu::xtal::Lattice(a, b, c);
}

cu::xtal::Lattice read_lattice()
{
    std::ifstream failstream("./FAIL/target.txt");
    Eigen::Matrix3d failmat;
    double val;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            failstream >> val;
            failmat(i, j) = val;
        }
    }
    return cu::xtal::Lattice(failmat);
}

void write_to_file(const cu::xtal::Lattice& lat, std::string filename)
{
    std::ofstream latstream;
    latstream.open(filename);
    latstream << lat.column_vector_matrix() << std::endl;
}

cu::xtal::Lattice fake_2d_as_lattice(const Eigen::Matrix2d& mat)
{
    Eigen::Matrix3d mat3 = Eigen::Matrix3d::Identity();
    mat3.block<2, 2>(0, 0) = mat;
    return mat3;
}

void write_all_lattices(const mush::MoireLattice& moire, int frame)
{
    const auto& L = moire.aligned_lattice;
    const auto& Lt = moire.rotated_lattice;

    const auto& K = moire.reciprocal_aligned_lattice;
    const auto& Kt = moire.reciprocal_rotated_lattice;

    cu::xtal::Lattice G(fake_2d_as_lattice(moire.full_reciprocal_difference));

    cu::xtal::Lattice Gz(fake_2d_as_lattice(moire.aligned_brillouin_zone_reciprocal_difference));
    cu::xtal::Lattice Gzt(fake_2d_as_lattice(moire.rotated_brillouin_zone_reciprocal_difference));

    const auto& M = moire.aligned_moire_lattice;
    const auto& Mt = moire.rotated_moire_lattice;

    /* std::cout << frame << " "; */
    /* std::cout << moire.input_degrees << " "; */
    /* std::cout << moire.is_within_brillouin_zone_overlap.at(&M).at(0) << " "; */
    /* std::cout << moire.is_within_brillouin_zone_overlap.at(&M).at(1) << " "; */
    /* std::cout << moire.is_within_brillouin_zone_overlap.at(&Mt).at(0) << " "; */
    /* std::cout << moire.is_within_brillouin_zone_overlap.at(&Mt).at(1) << "\n"; */

    /* xtal::Lattice Gd=make_moire_within_voronoi_overlap(G,K,Kt); */
    /* xtal::Lattice Md= Gd.reciprocal(); */

    std::string dir = "/home/mesto/programming/multishifter/tests/projects/frames/" + std::to_string(frame) + "/";

    write_to_file(L, dir + "L.txt");
    write_to_file(Lt, dir + "Lt.txt");
    write_to_file(K, dir + "K.txt");
    write_to_file(Kt, dir + "Kt.txt");
    write_to_file(G, dir + "G.txt");
    write_to_file(Gz, dir + "Gz.txt");
    write_to_file(Gzt, dir + "Gzt.txt");
    write_to_file(M, dir + "M.txt");
    write_to_file(Mt, dir + "Mt.txt");

    std::ofstream degstream(dir + "degrees.txt");
    degstream << moire.input_degrees << "\n";
    degstream << moire.is_within_brillouin_zone_overlap.at(&M).at(0) << "\n";
    degstream << moire.is_within_brillouin_zone_overlap.at(&M).at(1) << "\n";
    degstream << moire.is_within_brillouin_zone_overlap.at(&Mt).at(0) << "\n";
    degstream << moire.is_within_brillouin_zone_overlap.at(&Mt).at(1) << "\n";
    degstream.close();

    return;
}

} // namespace dirty

/* TEST_F(TwistTest, DirtyHacks) */
/* { */
/*     using namespace dirty; */

/*     std::ifstream failstream("./FAIL/target.txt"); */
/*     Eigen::Matrix3d failmat; */
/*     double val; */
/*     for(int i=0; i<3; ++i) */
/*     { */
/*         for(int j=0; j<3; ++j) */
/*         { */
/*             failstream>>val; */
/*             failmat(i,j)=val; */
/*         } */
/*     } */

/*     double degrees; */
/*     failstream>>degrees; */
/*     failstream.close(); */

/*     cu::xtal::Lattice faillat(failmat); */

/*     std::cout<<faillat.column_vector_matrix()<<"\n\n "<<degrees<<"\n"; */

/*     int frame = 5000; */
/*     mush::MoireLattice moire0(faillat, degrees); */
/*     write_all_lattices(moire0, frame); */
/* } */

/* TEST_F(TwistTest, DirtyHacks) */
/* { */
/*     using namespace dirty; */
/*     /1* auto L=make_hexagonal_lattice(); *1/ */
/*     /1* auto L=make_rectangular_lattice(); *1/ */
/*     /1* auto L = make_ugly_lattice(); *1/ */
/*     auto L = make_almost_ugly_lattice(); */
/*     /1* auto L = read_lattice(); *1/ */

/*     int frame = 1000; */
/*     for (int i = 0; i < 360; ++i) */
/*     { */

/*         double degrees = i + 0.25; */
/*         mush::MoireLattice moire0(L, degrees); */
/*         write_all_lattices(moire0, frame); */
/*         ++frame; */

/*         mush::MoireLattice moire1(L, degrees); */
/*         degrees = i + 0.75; */
/*         write_all_lattices(moire1, frame); */
/*         ++frame; */
/*     } */
/* } */
#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/xtal/frankenstein.hpp>
/* TEST_F(TwistTest, DirtyHacks) */
/* { */
/*     auto graph_path = autotools::input_filesdir / "graphene.vasp"; */
/*     auto graph = cu::xtal::Structure::from_poscar(graph_path); */

/*     std::vector<cu::xtal::Site> basis0; */
/*     std::vector<cu::xtal::Site> basis1; */

/*     for (const auto s : graph.basis_sites()) */
/*     { */
/*         cu::xtal::Coordinate c(s.cart()); */
/*         basis0.emplace_back(c, "C0"); */
/*         basis1.emplace_back(c, "C1"); */
/*     } */

/*     std::cout << "Degrees aligned_error rotated_error\n"; */
/*     int frame=1000; */
/*     for (int degrees = -360; degrees < 361; degrees += 1.0) */
/*     { */
/*         if(degrees%60==0) */
/*         { */
/*             continue; */
/*         } */

/*         MoireGenerator graph_moire(graph.lattice(), degrees); */

/*         dirty::write_all_lattices(graph_moire.moire, frame++); */
/*         using ZONE = MoireGenerator::ZONE; */
/*         using LATTICE = MoireGenerator::LATTICE; */

/*         auto [aligned_R, aligned_U] = cu::xtal::polar_decomposition(graph_moire.approximation_deformation(ZONE::ALIGNED, LATTICE::ALIGNED)); */
/*         auto [rot_R, rot_U] = cu::xtal::polar_decomposition(graph_moire.approximation_deformation(ZONE::ALIGNED, LATTICE::ROTATED)); */

/*         /1* std::cout<<rot_R<<"\n\n"; *1/ */

/*         auto rad_to_deg = [](double rad) { return 180 * rad / M_PI; }; */
/*         std::cout << degrees << "    "; */
/*         std::cout << rad_to_deg(std::asin(rot_R(1, 0))) << "    "; */
/*         std::cout << rad_to_deg(std::asin(aligned_R(1, 0))) << "    "; */
/*         std::cout<<std::endl; */
/*         /1* std::cout << rad_to_deg(rot_R(0, 0)) << "    "; *1/ */
/*         /1* std::cout << rad_to_deg(rot_R(0, 1)) << "    "; *1/ */
/*         /1* std::cout << rad_to_deg(rot_R(1, 0)) << "    "; *1/ */
/*         /1* std::cout << rad_to_deg(rot_R(1, 1)) << "\n"; *1/ */

/*         cu::xtal::Structure graph_aligned(graph.lattice(), basis0); */
/*         cu::xtal::Structure graph_rotated(graph.lattice(), basis1); */

/*         const auto& approx_aligned_lat = graph_moire.approximate_lattice(ZONE::ALIGNED, LATTICE::ALIGNED); */
/*         const auto& approx_rotated_lat = graph_moire.approximate_lattice(ZONE::ALIGNED, LATTICE::ROTATED); */

/*         const MoireApproximant::matrix_type& approx_aligned_transform = */
/*             graph_moire.approximate_moire_integer_transformation(ZONE::ALIGNED, LATTICE::ALIGNED); */
/*         const MoireApproximant::matrix_type& approx_rotated_transform = */
/*             graph_moire.approximate_moire_integer_transformation(ZONE::ALIGNED, LATTICE::ROTATED); */

/*         graph_aligned.set_lattice(approx_aligned_lat, cu::xtal::FRAC); */
/*         graph_rotated.set_lattice(approx_rotated_lat, cu::xtal::FRAC); */

/*         Lattice aligned_superlat = cu::xtal::make_superlattice(approx_aligned_lat, approx_aligned_transform.cast<int>()); */

/*         Lattice rotated_superlat = cu::xtal::make_superlattice(approx_rotated_lat, approx_rotated_transform.cast<int>()); */

/*         cu::xtal::write_poscar(graph_aligned, "./" + std::to_string(degrees) + "_aligned.vasp"); */
/*         cu::xtal::write_poscar(graph_rotated, "./" + std::to_string(degrees) + "_rotated.vasp"); */

/*         auto graph_bottom = cu::xtal::make_superstructure(graph_aligned, approx_aligned_transform.cast<int>()); */
/*         auto graph_top = cu::xtal::make_superstructure(graph_rotated, approx_rotated_transform.cast<int>()); */

/*         auto graph_twist = cu::xtal::frankenstein::stack({graph_bottom, graph_top}); */
/*         cu::xtal::write_poscar(graph_twist, "./" + std::to_string(degrees) + "_graphstack.vasp"); */
/*         cu::xtal::write_poscar(graph_bottom, "./" + std::to_string(degrees) + "_bottom.vasp"); */
/*         cu::xtal::write_poscar(graph_top, "./" + std::to_string(degrees) + "_top.vasp"); */
/*     } */
/* } */

TEST_F(TwistTest, DirtyHacks)
{
    auto graph_path = autotools::input_filesdir / "graphene.vasp";
    auto graph = cu::xtal::Structure::from_poscar(graph_path);

    for (int degrees = 0; degrees < 360; degrees += 1)
    {
        if(degrees%60==0)
        {
            continue;
        }

        using ZONE = MoireStructureGenerator::ZONE;
        using LATTICE = MoireStructureGenerator::LATTICE;

        MoireStructureGenerator gen(graph,degrees);

        auto graph_bottom=gen.layer(ZONE::ALIGNED,LATTICE::ALIGNED);
        auto graph_top=gen.layer(ZONE::ALIGNED,LATTICE::ROTATED);

        auto graph_twist = cu::xtal::frankenstein::stack({graph_bottom, graph_top});
        cu::xtal::write_poscar(graph_twist, "./" + std::to_string(degrees) + "_graphstack.vasp");
        cu::xtal::write_poscar(graph_bottom, "./" + std::to_string(degrees) + "_bottom.vasp");
        cu::xtal::write_poscar(graph_top, "./" + std::to_string(degrees) + "_top.vasp");
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
