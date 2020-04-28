#include "./twist.hpp"
#include "casmutils/xtal/coordinate.hpp"
#include "casmutils/xtal/site.hpp"
#include "casmutils/xtal/structure_tools.hpp"
#include <casmutils/xtal/lattice.hpp>
#include <cassert>
#include <cmath>
#include <multishift/definitions.hpp>

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
Eigen::Matrix3d slab_unit_vectors(const cu::xtal::Lattice& slab)
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
cu::xtal::Lattice make_transformed_lattice(const cu::xtal::Lattice& lat, const Eigen::Matrix3d& transform)
{
    return cu::xtal::Lattice(transform * lat.column_vector_matrix());
}

cu::xtal::Lattice make_aligned_lattice(const cu::xtal::Lattice& lat)
{
    Eigen::Matrix3d lat_span_to_standard = slab_unit_vectors(lat).inverse();
    cu::xtal::Lattice aligned_lat = make_transformed_lattice(lat, lat_span_to_standard);
    assert(almost_equal((aligned_lat.a().cross(aligned_lat.b())).normalized(), Eigen::Vector3d(0, 0, 1), 1e-10));
    return aligned_lat;
}

Eigen::Matrix3d make_twist_rotation_matrix(const cu::xtal::Lattice& lat, double degrees)
{
    double rad = M_PI * degrees / 180.0;
    Eigen::Matrix3d z_rot_mat;
    z_rot_mat << std::cos(rad), -std::sin(rad), 0, std::sin(rad), std::cos(rad), 0, 0, 0, 1;

    Eigen::Matrix3d standard_to_lat_span = slab_unit_vectors(lat);
    Eigen::Matrix3d lat_span_to_standard = standard_to_lat_span.inverse();

    return standard_to_lat_span * z_rot_mat * lat_span_to_standard;
}

cu::xtal::Lattice make_twisted_lattice(const cu::xtal::Lattice& lat, double degrees)
{
    Eigen::Matrix3d twist_matrix = make_twist_rotation_matrix(lat, degrees);
    return make_transformed_lattice(lat, twist_matrix);
}

cu::xtal::Lattice make_aligned_moire_lattice(const cu::xtal::Lattice& lat, double degrees)
{
    auto make_2d_column_matrix=[](const cu::xtal::Lattice lat3d){
    Eigen::Vector2d a2d(lat3d.a()(0),lat3d.a()(1));
    Eigen::Vector2d b2d(lat3d.b()(0),lat3d.b()(1));

    Eigen::Matrix2d lat_mat2d;
    lat_mat2d.col(0)=a2d;
    lat_mat2d.col(1)=b2d;

    return lat_mat2d;
    };

    cu::xtal::Lattice aligned_lat = make_aligned_lattice(lat);
    assert(almost_equal(aligned_lat.a()(2),0.0,1e-10));
    assert(almost_equal(aligned_lat.a()(1),0.0,1e-10));
    assert(almost_equal(aligned_lat.b()(2),0.0,1e-10));

    cu::xtal::Lattice rotated_lat=make_twisted_lattice(aligned_lat,degrees);
    assert(almost_equal(rotated_lat.a()(2),0.0,1e-10));
    assert(almost_equal(rotated_lat.b()(2),0.0,1e-10));

    Eigen::Matrix2d aligned_mat2d=make_2d_column_matrix(aligned_lat);
    Eigen::Matrix2d rotated_mat2d=make_2d_column_matrix(rotated_lat);

    Eigen::Matrix2d moire_recip_mat=rotated_mat2d.inverse()-aligned_mat2d.inverse();
    Eigen::Matrix2d moire_mat=moire_recip_mat.inverse();

    Eigen::Matrix3d moire_mat3d(Eigen::Matrix3d::Identity());
    moire_mat3d.block<2,2>(0,0)=moire_mat;

    return cu::xtal::Lattice(moire_mat3d);
}
} // namespace mush
