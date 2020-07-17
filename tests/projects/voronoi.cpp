#include <casm/crystallography/Coordinate.hh>
#include <casm/crystallography/Lattice.hh>
#include <casm/crystallography/BasicStructure.hh>
#include <casm/external/Eigen/Dense>
#include <fstream>
#include <iostream>
#include <string>

using namespace CASM;

double degrees_to_radians(double angle)
{
    return angle*M_PI/180.0;
}

double radians_to_degrees(double rad)
{
    return rad*108/M_PI;
}

Eigen::Matrix3d make_rotation_matrix(double angle)
{
    double rad=degrees_to_radians(angle);
    Eigen::Matrix3d rotation=Eigen::Matrix3d::Identity();

    double c=std::cos(rad);
    double s=std::sin(rad);

    rotation(0,0)=c;
    rotation(0,1)=-s;
    rotation(1,0)=s;
    rotation(1,1)=c;

    return rotation;
}

xtal::Lattice make_hexagonal_lattice()
{
    Eigen::Vector3d a(0,3,0);
    Eigen::Matrix3d R=make_rotation_matrix(60);
    Eigen::Vector3d b=R*a;
    Eigen::Vector3d c(0,0,1);

    return xtal::Lattice(a,b,c);
}

xtal::Lattice make_rectangular_lattice()
{
    Eigen::Vector3d a(4,0,0);
    Eigen::Vector3d b(0,3.749372649,0);
    Eigen::Vector3d c(0,0,1);
    return xtal::Lattice(a,b,c);
}

xtal::Lattice make_ugly_lattice()
{
    Eigen::Vector3d a(4,0,0);
    Eigen::Vector3d b(3.7,0,0);
    Eigen::Matrix3d R=make_rotation_matrix(40);
    b=R*b;
    Eigen::Vector3d c(0,0,1);
    return xtal::Lattice(a,b,c);
}

xtal::Lattice make_transformed_lattice(const xtal::Lattice& lat, const Eigen::Matrix3d& transformation)
{
    return xtal::Lattice(transformation*lat.lat_column_mat());
}

void write_to_file(const xtal::Lattice& lat, std::string filename)
{
    std::ofstream latstream;
    latstream.open(filename);
    latstream<<lat.lat_column_mat()<<std::endl;
}

xtal::Lattice make_moire_within_voronoi(const xtal::Lattice& moire, const xtal::Lattice& bz_lat)
{
    xtal::Coordinate Mz_a(moire.lat_column_mat().col(0),bz_lat,CART);
    xtal::Coordinate Mz_b(moire.lat_column_mat().col(1),bz_lat,CART);
    Mz_a.voronoi_within();
    Mz_b.voronoi_within();

    return xtal::Lattice(Mz_a.const_cart(),Mz_b.const_cart(),moire.lat_column_mat().col(2));
}

int main(int argc, char* argv[])
{
    auto L=make_hexagonal_lattice();
    /* auto L=make_rectangular_lattice(); */
    /* auto L=make_ugly_lattice(); */
    double angle=std::stod(argv[1]);
    Eigen::Matrix3d rot_m=make_rotation_matrix(angle);
    auto Lt=make_transformed_lattice(L,rot_m);

    auto K=L.reciprocal();
    auto Kt=Lt.reciprocal();

    Eigen::Matrix3d G_mat=Kt.lat_column_mat()-K.lat_column_mat();
    G_mat(2,2)=1;
    xtal::Lattice G(G_mat);

    xtal::Lattice Gz=make_moire_within_voronoi(G,K);
    xtal::Lattice Gzt=make_moire_within_voronoi(G,Kt);

    xtal::Lattice M=Gz.reciprocal();
    xtal::Lattice Mt=Gzt.reciprocal();

    write_to_file(L,"L.txt");
    write_to_file(Lt,"Lt.txt");
    write_to_file(K,"K.txt");
    write_to_file(Kt,"Kt.txt");
    write_to_file(G,"G.txt");
    write_to_file(Gz,"Gz.txt");
    write_to_file(Gzt,"Gzt.txt");
    write_to_file(M,"M.txt");
    write_to_file(Mt,"Mt.txt");

    return 0;
}
