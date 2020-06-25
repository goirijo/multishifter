#include "../../autotools.hh"
#include "casmutils/xtal/structure.hpp"
#include "casmutils/xtal/structure_tools.hpp"
#include "multishift/slab.hpp"
#include <cmath>
#include <fstream>
#include <gtest/gtest.h>

#include <multishift/definitions.hpp>
#include <multishift/twist.hpp>
#include <tuple>
#include <vector>

namespace casmutils::xtal
{
xtal::Lattice make_sliced_lattice(const xtal::Lattice& unit_lattice, const Eigen::Vector3i& miller_indexes);

namespace frankenstein
{
xtal::Structure stack(const std::vector<xtal::Structure>& sub_strucs);
}
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
            auto [moire_lat,aligned_lat,rot_lat]=make_aligned_moire_lattice(start_lat, angle);

            Eigen::Matrix3d aligned_to_moire_transform=aligned_lat.column_vector_matrix().inverse()*moire_lat.column_vector_matrix();
            Eigen::Matrix3d aligned_to_moire_transform_diff=aligned_to_moire_transform-aligned_to_moire_transform.unaryExpr([](double x){return std::round(x);});

            Eigen::Matrix3d rot_to_moire_transform=rot_lat.column_vector_matrix().inverse()*moire_lat.column_vector_matrix();
            Eigen::Matrix3d rot_to_moire_transform_diff=rot_to_moire_transform-rot_to_moire_transform.unaryExpr([](double x){return std::round(x);});

            // TODO: How to test the moire lattice? Plotting the two lattices seems
            // to give the correct results. I THINK this check makes sense
            EXPECT_TRUE(rot_to_moire_transform_diff.isApprox(rot_to_moire_transform_diff));

            //This block is for using a dumb python script for visualization
            dumpstream<<"****************************\n\n";
            dumpstream<<aligned_lat.column_vector_matrix()<<"\n\n";
            dumpstream<<rot_lat.column_vector_matrix()<<"\n\n";
            dumpstream<<moire_lat.column_vector_matrix()<<"\n\n";
        }
    }
    dumpstream.close();
}

TEST_F(TwistTest, ApproximantMoireLattice)
{
    std::ofstream dumpstream("./apprixmoirelatticeout.txt");
    for (const Lattice& start_lat : sliced_lattices)
    {
        for (double angle : {-3.00, -2.0, -1.0, -0.50, 0.50, 2.00, 5.00, 22.0})
        {
            auto [moire_lat,aligned_lat,rot_lat]=make_approximant_moire_lattice(start_lat, angle);

            Eigen::Matrix3d aligned_to_moire_transform=aligned_lat.column_vector_matrix().inverse()*moire_lat.column_vector_matrix();
            Eigen::Matrix3d aligned_to_moire_transform_diff=aligned_to_moire_transform-aligned_to_moire_transform.unaryExpr([](double x){return std::round(x);});

            Eigen::Matrix3d rot_to_moire_transform=rot_lat.column_vector_matrix().inverse()*moire_lat.column_vector_matrix();
            Eigen::Matrix3d rot_to_moire_transform_diff=rot_to_moire_transform-rot_to_moire_transform.unaryExpr([](double x){return std::round(x);});

            EXPECT_TRUE(almost_zero(aligned_to_moire_transform_diff.block<2,2>(0,0),1e-8));
            EXPECT_TRUE(almost_zero(rot_to_moire_transform_diff.block<2,2>(0,0),1e-8));

            //This block is for using a dumb python script for visualization
            dumpstream<<"****************************\n\n";
            dumpstream<<aligned_lat.column_vector_matrix()<<"\n\n";
            dumpstream<<rot_lat.column_vector_matrix()<<"\n\n";
            dumpstream<<moire_lat.column_vector_matrix()<<"\n\n";

            if(!almost_zero(aligned_to_moire_transform_diff.block<2,2>(0,0),1e-8) || !almost_zero(rot_to_moire_transform_diff.block<2,2>(0,0),1e-8))
            {
                std::cout<<"ALERT "<<angle<<"\n";

                std::cout<<"DEBUGGING: start_lat.column_vector_matrix().determinant() is "<<start_lat.column_vector_matrix().determinant()<<std::endl;
                std::cout<<"DEBUGGING: aligned_to_moire_transform_diff is "<<aligned_to_moire_transform_diff<<std::endl;
                std::cout<<"DEBUGGING: rot_to_moire_transform_diff is "<<rot_to_moire_transform_diff<<std::endl;
                
                

                auto tmp=mush::make_aligned_moire_lattice(start_lat,angle);
                std::cout<<std::get<1>(tmp).column_vector_matrix()<<"\n\n";

                std::cout<<mush::make_twist_rotation_matrix(start_lat,angle)<<"\n\n";
                std::cout<<mush::make_twist_rotation_matrix(start_lat,angle).inverse()<<"\n\n";
                
                std::cout<<start_lat.column_vector_matrix()<<"\n\n";
                std::cout<<aligned_lat.column_vector_matrix()<<"\n\n";
                std::cout<<aligned_to_moire_transform_diff<<"\n\n";
            }
        }
    }
    dumpstream.close();
}

TEST_F(TwistTest, MoireApproximantVectorMatch)
{
    for (const Lattice& start_lat : sliced_lattices)
    {
        for (double angle : {-3.00, -2.0, -1.0, -0.50, 0.50, 2.00, 5.00, 22.0})
        {
            MoireApproximant moire(start_lat,angle);
            Lattice aligned_super=cu::xtal::make_superlattice(moire.approximate_lattices[0],moire.approximate_moire_integer_transformations[0].cast<int>());
            Lattice rot_super=cu::xtal::make_superlattice(moire.approximate_lattices[1],moire.approximate_moire_integer_transformations[1].cast<int>());

            std::cout<<"DEBUGGING: aligned_super.column_vector_matrix()-rot_super.column_vector_matrix() is \n"<<aligned_super.column_vector_matrix()-rot_super.column_vector_matrix()<<std::endl<<std::endl;;
            

            almost_equal(aligned_super.a(),rot_super.a());
            almost_equal(aligned_super.b(),rot_super.b());
        }
    }
}

TEST_F(TwistTest, MoireApproximantAlignedRotatedMatch)
{
    for (const Lattice& start_lat : sliced_lattices)
    {
        for (double angle : {-3.00, -2.0, -1.0, -0.50, 0.50, 2.00, 5.00, 22.0})
        {
            //TODO: Make sure that the aligned and rotated lattices
            //(both perfect and approximated) are exclusively related
            //by a rotation
        }
    }
}

#include <casmutils/xtal/frankenstein.hpp>

TEST_F(TwistTest, MoireApproximantRotationEffect)
{
    for (const Lattice& start_lat : sliced_lattices)
    {
        for (double angle : {-3.00, -2.0, -1.0, -0.50, 0.50, 2.00, 5.00, 22.0})
        {
            MoireApproximant moire(start_lat,angle);
            auto [aligned_R,aligned_U]=cu::xtal::polar_decomposition(moire.approximation_deformations[0]);
            auto [rot_R,rot_U]=cu::xtal::polar_decomposition(moire.approximation_deformations[1]);

            std::cout<<moire.approximate_moire_integer_transformations[0]<<"\n\n";
            std::cout<<moire.approximate_moire_integer_transformations[1]<<"\n\n";

            for(const Eigen::Matrix3d& M : {aligned_R,aligned_U,rot_R,rot_U})
            {
                /* std::cout<<moire.approximation_deformations[0]<<std::endl<<std::endl; */
                /* std::cout<<moire.approximation_deformations[1]<<std::endl<<std::endl; */

                /* std::cout<<moire.aligned_lattice.column_vector_matrix()<<"\n\n"; */
                /* std::cout<<moire.approximate_lattices[0].column_vector_matrix()<<"\n\n"; */
                /* std::cout<<moire.approximation_deformations[0]<<"\n\n"; */

                /* std::cout<<moire.aligned_lattice.column_vector_matrix()-moire.approximate_lattices[0].column_vector_matrix(); */
                /* std::cout<<"\n\n"; */
                /* std::cout<<moire.rotated_lattice.column_vector_matrix()-moire.approximate_lattices[1].column_vector_matrix(); */
                /* std::cout<<"\n\n"; */

                /* std::cout<<M.determinant()<<std::endl; */
                /* std::cout<<M<<"\n\n"; */

                /* EXPECT_TRUE(almost_equal(M(0,1),0.0,1e-8)); */
                /* EXPECT_TRUE(almost_equal(M(0,2),0.0,1e-8)); */
                /* EXPECT_TRUE(almost_equal(M(2,0),0.0,1e-8)); */
                /* EXPECT_TRUE(almost_equal(M(2,1),0.0,1e-8)); */
                /* EXPECT_TRUE(almost_equal(M(2,1),1.0,1e-8)); */
            }
            std::cout<<"***************\n";
        }
    }

    //Just for testing
    auto hcp_path=autotools::input_filesdir/"hcp.vasp";
    auto hcp=cu::xtal::Structure::from_poscar(hcp_path);

    MoireApproximant hcp_moire(hcp.lattice(),5.0);
    auto hcp_aligned=hcp;
    auto hcp_rotated=hcp;

    hcp_aligned.set_lattice(hcp_moire.approximate_lattices[0],cu::xtal::FRAC);
    hcp_rotated.set_lattice(hcp_moire.approximate_lattices[1],cu::xtal::FRAC);


    
    Lattice aligned_superlat=cu::xtal::make_superlattice(hcp_moire.approximate_lattices[0],hcp_moire.approximate_moire_integer_transformations[0].cast<int>());
    Lattice rotated_superlat=cu::xtal::make_superlattice(hcp_moire.approximate_lattices[1],hcp_moire.approximate_moire_integer_transformations[1].cast<int>());

    std::cout<<"&&&&&&&&&&&&&&&&&&\n";
    std::cout<<aligned_superlat.column_vector_matrix()<<"\n\n";
    std::cout<<rotated_superlat.column_vector_matrix()<<"\n\n";
    std::cout<<"&&&&&&&&&&&&&&&&&&\n";






    cu::xtal::write_poscar(hcp_aligned,"./aligned.vasp");
    cu::xtal::write_poscar(hcp_rotated,"./rotated.vasp");

    auto hcp_bottom=cu::xtal::make_superstructure(hcp_aligned,hcp_moire.approximate_moire_integer_transformations[0].cast<int>());
    auto hcp_top=cu::xtal::make_superstructure(hcp_aligned,hcp_moire.approximate_moire_integer_transformations[1].cast<int>());

    auto hcp_twist=cu::xtal::frankenstein::stack({hcp_bottom,hcp_top});
    cu::xtal::write_poscar(hcp_twist,"./hcpstack.vasp");
    cu::xtal::write_poscar(hcp_bottom,"./bottom.vasp");
    cu::xtal::write_poscar(hcp_top,"./top.vasp");
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
