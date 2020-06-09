#include "./twist.hpp"
#include "casmutils/xtal/coordinate.hpp"
#include "casmutils/xtal/site.hpp"
#include "casmutils/xtal/structure_tools.hpp"
#include "multishift/slab.hpp"
#include <casmutils/xtal/lattice.hpp>
#include <cassert>
#include <cmath>
#include <multishift/definitions.hpp>
#include <tuple>

namespace
{
    using namespace mush;
    Eigen::Matrix2d make_2d_column_matrix(const Lattice lat3d) {
        Eigen::Vector2d a2d(lat3d.a()(0), lat3d.a()(1));
        Eigen::Vector2d b2d(lat3d.b()(0), lat3d.b()(1));

        Eigen::Matrix2d lat_mat2d;
        lat_mat2d.col(0) = a2d;
        lat_mat2d.col(1) = b2d;

        return lat_mat2d;
    };

    Eigen::Matrix3d make_3d_column_matrix(const Eigen::Matrix2d& lat2d)
    {
        Eigen::Matrix3d lat3d=Eigen::Matrix3d::Identity();
        lat3d.block<2,2>(0,0)=lat2d;
        return lat3d;
    }
}

namespace casmutils
{
namespace xtal
{
Lattice make_reciprocal(const Lattice& real_lattice)
{
    auto casm_reciprocal = real_lattice.__get().reciprocal();
    // TODO: Move constructor for lattice?
    return Lattice(casm_reciprocal);
}
} // namespace xtal
} // namespace casmutils

#include <casmutils/xtal/structure.hpp>
namespace mush
{
/// Returns orhtogonal unit vectors oriented such that the point along the
/// a vector, the ab plane normal, and whatever is perpendicular to that
Eigen::Matrix3d slab_unit_vectors(const Lattice& slab)
{
    Eigen::Matrix3d slab_span;
    slab_span.col(0) = slab.a().normalized();
    slab_span.col(2) = slab.a().cross(slab.b()).normalized();
    slab_span.col(1) = slab_span.col(2).cross(slab_span.col(0)).normalized();

    assert(slab_span.determinant() * slab.column_vector_matrix().determinant() > 0);
    assert(almost_equal(std::abs(slab_span.determinant()), 1.0, 1e-10));

    return slab_span;
}

/// Applies the given transformation the the *column* vector matrix representation
/// of the lattice
Lattice make_transformed_lattice(const Lattice& lat, const Eigen::Matrix3d& transform)
{
    return Lattice(transform * lat.column_vector_matrix());
}

Lattice make_aligned_lattice(const Lattice& lat)
{
    Eigen::Matrix3d lat_span_to_standard = slab_unit_vectors(lat).inverse();
    Lattice aligned_lat = make_transformed_lattice(lat, lat_span_to_standard);
    assert(almost_equal((aligned_lat.a().cross(aligned_lat.b())).normalized(), Eigen::Vector3d(0, 0, 1), 1e-10));
    return aligned_lat;
}

Eigen::Matrix3d make_twist_rotation_matrix(const Lattice& lat, double degrees)
{
    double rad = M_PI * degrees / 180.0;
    Eigen::Matrix3d z_rot_mat;
    z_rot_mat << std::cos(rad), -std::sin(rad), 0, std::sin(rad), std::cos(rad), 0, 0, 0, 1;

    Eigen::Matrix3d standard_to_lat_span = slab_unit_vectors(lat);
    Eigen::Matrix3d lat_span_to_standard = standard_to_lat_span.inverse();

    return standard_to_lat_span * z_rot_mat * lat_span_to_standard;
}

Lattice make_twisted_lattice(const Lattice& lat, double degrees)
{
    Eigen::Matrix3d twist_matrix = make_twist_rotation_matrix(lat, degrees);
    return make_transformed_lattice(lat, twist_matrix);
}

std::tuple<Lattice, Lattice, Lattice> make_aligned_moire_lattice(const Lattice& lat, double degrees)
{

    const Lattice aligned_lat = make_aligned_lattice(lat);
    assert(almost_equal(aligned_lat.a()(2), 0.0, 1e-10));
    assert(almost_equal(aligned_lat.a()(1), 0.0, 1e-10));
    assert(almost_equal(aligned_lat.b()(2), 0.0, 1e-10));

    const Lattice rotated_lat = make_twisted_lattice(aligned_lat, degrees);
    assert(almost_equal(rotated_lat.a()(2), 0.0, 1e-10));
    assert(almost_equal(rotated_lat.b()(2), 0.0, 1e-10));

    const Eigen::Matrix2d aligned_mat2d = ::make_2d_column_matrix(aligned_lat);
    const Eigen::Matrix2d rotated_mat2d = ::make_2d_column_matrix(rotated_lat);

    const Eigen::Matrix2d moire_recip_mat = rotated_mat2d.inverse() - aligned_mat2d.inverse();
    const Eigen::Matrix2d moire_mat = moire_recip_mat.inverse();

    Eigen::Matrix3d moire_mat3d(Eigen::Matrix3d::Identity());
    moire_mat3d.block<2, 2>(0, 0) = moire_mat;

    return std::make_tuple(moire_mat3d, aligned_lat, rotated_lat);
}

std::tuple<Lattice, Lattice, Lattice> make_approximant_moire_lattice(const Lattice& lat, double degrees)
{
    auto [moire_lat, aligned_lat, rot_lat] = make_aligned_moire_lattice(lat, degrees);

    Eigen::Matrix2d moire_lat_2d=::make_2d_column_matrix(moire_lat);
    Eigen::Matrix2d aligned_lat_2d=::make_2d_column_matrix(aligned_lat);
    Eigen::Matrix2d rot_lat_2d=::make_2d_column_matrix(rot_lat);

    Eigen::Matrix2d aligned_to_moire_transform = aligned_lat_2d.inverse() * moire_lat_2d;
    Eigen::Matrix2d aligned_to_moire_transform_round=aligned_to_moire_transform.unaryExpr([](double x){return std::round(x);});

    Eigen::Matrix2d rot_to_moire_transform = rot_lat_2d.inverse() * moire_lat_2d;
    Eigen::Matrix2d rot_to_moire_transform_round=rot_to_moire_transform.unaryExpr([](double x){return std::round(x);});

    Eigen::Matrix2d aligned_superlattice_2d=aligned_lat_2d*aligned_to_moire_transform_round;
    Eigen::Matrix2d rot_superlattice_2d=rot_lat_2d*rot_to_moire_transform_round;

    Eigen::Matrix2d approx_superlattice_2d=(aligned_superlattice_2d+rot_superlattice_2d)/2.0;
    Eigen::Matrix2d approx_aligned_2d=approx_superlattice_2d*aligned_to_moire_transform_round.inverse();
    Eigen::Matrix2d approx_rot_2d=approx_superlattice_2d*rot_to_moire_transform_round.inverse();

    Eigen::Matrix3d approx_superlattice_col_mat=::make_3d_column_matrix(approx_superlattice_2d);
    //The moire lattice has zero c vector, but we still want to keep the c vector around for the aligned and rotated lattices
    Eigen::Matrix3d approx_aligned_col_mat=::make_3d_column_matrix(approx_aligned_2d);
    approx_aligned_col_mat.col(2)=aligned_lat.c();
    Eigen::Matrix3d approx_rot_col_mat=::make_3d_column_matrix(approx_rot_2d);
    approx_rot_col_mat.col(2)=rot_lat.c();

    return std::make_tuple(approx_superlattice_col_mat,approx_aligned_col_mat,approx_rot_col_mat);
}
} // namespace mush
