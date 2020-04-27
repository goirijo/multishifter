#include "./twist.hpp"
#include "casmutils/xtal/coordinate.hpp"
#include "casmutils/xtal/site.hpp"
#include "casmutils/xtal/structure_tools.hpp"
#include <casmutils/xtal/lattice.hpp>
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

    assert(slab_span.determinant()*slab.column_vector_matrix().determinant() > 0);
    assert(almost_equal(std::abs(slab_span.determinant()),1.0,1e-10));

    return slab_span;
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

    const Eigen::Matrix3d& lat_column_mat = lat.column_vector_matrix();
    Eigen::Matrix3d twisted_lat_column_mat = twist_matrix * lat_column_mat;

    return cu::xtal::Lattice(twisted_lat_column_mat);
}
} // namespace mush
