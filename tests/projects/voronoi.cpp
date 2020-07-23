#include <casm/crystallography/BasicStructure.hh>
#include <casm/crystallography/Coordinate.hh>
#include <casm/crystallography/Lattice.hh>
#include <casm/external/Eigen/Dense>
#include <fstream>
#include <iostream>
#include <string>

using namespace CASM;

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

xtal::Lattice make_hexagonal_lattice()
{
    Eigen::Vector3d a(0, 3, 0);
    Eigen::Matrix3d R = make_rotation_matrix(60);
    Eigen::Vector3d b = R * a;
    Eigen::Vector3d c(0, 0, 1);

    return xtal::Lattice(a, b, c);
}

xtal::Lattice make_rectangular_lattice()
{
    Eigen::Vector3d a(5, 0, 0);
    Eigen::Vector3d b(0, 2, 0);
    Eigen::Vector3d c(0, 0, 1);
    return xtal::Lattice(a, b, c);
}

xtal::Lattice make_ugly_lattice()
{
    Eigen::Vector3d a(3, 0, 0);
    Eigen::Vector3d b(3.5, 0, 0);
    Eigen::Matrix3d R = make_rotation_matrix(50);
    b = R * b;
    Eigen::Vector3d c(0, 0, 1);
    return xtal::Lattice(a, b, c);
}

xtal::Lattice make_transformed_lattice(const xtal::Lattice& lat, const Eigen::Matrix3d& transformation)
{
    return xtal::Lattice(transformation * lat.lat_column_mat());
}

void write_to_file(const xtal::Lattice& lat, std::string filename)
{
    std::ofstream latstream;
    latstream.open(filename);
    latstream << lat.lat_column_mat() << std::endl;
}

xtal::Lattice make_moire_within_voronoi(const xtal::Lattice& moire, const xtal::Lattice& bz_lat)
{
    xtal::Coordinate Mz_a(moire.lat_column_mat().col(0), bz_lat, CART);
    xtal::Coordinate Mz_b(moire.lat_column_mat().col(1), bz_lat, CART);
    Mz_a.voronoi_within();
    Mz_b.voronoi_within();

    return xtal::Lattice(Mz_a.const_cart(), Mz_b.const_cart(), moire.lat_column_mat().col(2));
}

xtal::Lattice minimize_lattice_vectors(const xtal::Lattice lat1, const xtal::Lattice& lat2)
{
    Eigen::Matrix3d col_lat_mat;
    const Eigen::Matrix3d& mat1=lat1.lat_column_mat();
    const Eigen::Matrix3d& mat2=lat2.lat_column_mat();
    for(int i=0; i<3; ++i)
    {
        col_lat_mat.col(i)=mat1.col(i).norm()<mat2.col(i).norm()?mat1.col(i):mat2.col(i);
    }
    return xtal::Lattice(col_lat_mat);
}

Eigen::Vector3d bring_within_voronoi(const Eigen::Vector3d& v, const xtal::Lattice& lat)
{
    xtal::Coordinate c(v,lat,CART);
    c.voronoi_within();
    return c.const_cart();
}

bool is_within_voronoi(const Eigen::Vector3d& v, const xtal::Lattice& lat)
{
    Eigen::Vector3d v_within=bring_within_voronoi(v,lat);
    return almost_equal(v,v_within);
}

xtal::Lattice make_moire_within_voronoi_overlap(const xtal::Lattice& moire, const xtal::Lattice& bz_lat1, const xtal::Lattice& bz_lat2)
{
    Eigen::Vector3d Ma=moire.lat_column_mat().col(0);
    Eigen::Vector3d Mb=moire.lat_column_mat().col(1);

    Eigen::Vector3d Ma_within=Ma;
    Eigen::Vector3d Mb_within=Mb;

    if(!is_within_voronoi(Ma_within,bz_lat1))
    {
        Ma_within=bring_within_voronoi(Ma,bz_lat1);
    }

    if(!is_within_voronoi(Ma_within,bz_lat2))
    {
        Ma_within=bring_within_voronoi(Ma,bz_lat2);
    }

    assert(is_within_voronoi(Ma_within,bz_lat1)&&is_within_voronoi(Ma_within,bz_lat2));

    if(!is_within_voronoi(Mb_within,bz_lat1))
    {
        Mb_within=bring_within_voronoi(Mb,bz_lat1);
    }

    if(!is_within_voronoi(Mb_within,bz_lat2))
    {
        Mb_within=bring_within_voronoi(Mb,bz_lat2);
    }

    assert(is_within_voronoi(Mb_within,bz_lat1)&&is_within_voronoi(Mb_within,bz_lat2));

    return xtal::Lattice(Ma_within,Mb_within,moire.lat_column_mat().col(2));
}

xtal::Lattice make_moire_within_voronoi_double(const xtal::Lattice& moire, const xtal::Lattice& bz_lat1, const xtal::Lattice& bz_lat2)
{
    xtal::Coordinate Mz_a(moire.lat_column_mat().col(0), bz_lat1, CART);
    xtal::Coordinate Mz_b(moire.lat_column_mat().col(1), bz_lat1, CART);
    Mz_a.voronoi_within();
    Mz_b.voronoi_within();

    Mz_a.set_lattice(bz_lat2,CART);
    Mz_b.set_lattice(bz_lat2,CART);
    Mz_a.voronoi_within();
    Mz_b.voronoi_within();

    /* Mz_a.set_lattice(bz_lat1,CART); */
    /* Mz_b.set_lattice(bz_lat1,CART); */
    /* Mz_a.voronoi_within(); */
    /* Mz_b.voronoi_within(); */

    return xtal::Lattice(Mz_a.const_cart(), Mz_b.const_cart(), moire.lat_column_mat().col(2));
}


void write_all_lattices(const xtal::Lattice& L, double degrees, int frame)
{
    Eigen::Matrix3d rot_m = make_rotation_matrix(degrees);
    auto Lt = make_transformed_lattice(L, rot_m);

    auto K = L.reciprocal();
    auto Kt = Lt.reciprocal();

    Eigen::Matrix3d G_mat = Kt.lat_column_mat() - K.lat_column_mat();
    G_mat(2, 2) = 1;
    xtal::Lattice G(G_mat);

    xtal::Lattice Gz = make_moire_within_voronoi(G, K);
    xtal::Lattice Gzt = make_moire_within_voronoi(G, Kt);

    xtal::Lattice Gzx = make_moire_within_voronoi_double(G,K,Kt);
    xtal::Lattice Gzx2 = make_moire_within_voronoi_double(G,Kt,K);

    xtal::Lattice M = Gz.reciprocal();
    xtal::Lattice Mt = Gzt.reciprocal();

    xtal::Lattice Mx = Gzx.reciprocal();
    xtal::Lattice Mx2 = Gzx2.reciprocal();

    std::cout<<frame<<" ";
    std::cout<<degrees<<" ";
    std::cout<<is_within_voronoi(Gz[0],Kt)<<" ";
    std::cout<<is_within_voronoi(Gz[1],Kt)<<" ";
    std::cout<<is_within_voronoi(Gzt[0],K)<<" ";
    std::cout<<is_within_voronoi(Gzt[1],K)<<"\n";
    
    /* xtal::Lattice Gd=make_moire_within_voronoi_overlap(G,K,Kt); */
    /* xtal::Lattice Md= Gd.reciprocal(); */

    std::string dir = "frames/" + std::to_string(frame) + "/";

    write_to_file(L, dir + "L.txt");
    write_to_file(Lt, dir + "Lt.txt");
    write_to_file(K, dir + "K.txt");
    write_to_file(Kt, dir + "Kt.txt");
    write_to_file(G, dir + "G.txt");
    write_to_file(Gz, dir + "Gz.txt");
    write_to_file(Gzt, dir + "Gzt.txt");
    write_to_file(Gzx, dir + "Gzx.txt");
    write_to_file(Gzx2, dir + "Gzx2.txt");
    /* write_to_file(Gd, dir + "Gd.txt"); */
    write_to_file(M, dir + "M.txt");
    write_to_file(Mt, dir + "Mt.txt");
    write_to_file(Mx, dir + "Mx.txt");
    write_to_file(Mx2, dir + "Mx2.txt");
    /* write_to_file(Md, dir + "Md.txt"); */


    std::ofstream degstream(dir+"degrees.txt");
    degstream<<degrees<<"\n";
    degstream<<is_within_voronoi(Gz[0],Kt)<<"\n";
    degstream<<is_within_voronoi(Gz[1],Kt)<<"\n";
    degstream<<is_within_voronoi(Gzt[0],K)<<"\n";
    degstream<<is_within_voronoi(Gzt[1],K)<<"\n";
    degstream.close();

    return;
}

int main(int argc, char* argv[])
{
    /* auto L=make_hexagonal_lattice(); */
    auto L=make_rectangular_lattice();
    /* auto L = make_ugly_lattice(); */

    int frame = 1000;
    for (int i = 0; i < 360; ++i)
    {
        double degrees = i + 0.2;
        write_all_lattices(L, degrees, frame);
        ++frame;

        degrees = i + 0.7;
        write_all_lattices(L, degrees, frame);
        ++frame;
    }

    return 0;
}
