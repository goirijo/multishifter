#include "./twist.hpp"
#include "casmutils/xtal/coordinate.hpp"
#include "casmutils/xtal/site.hpp"
#include "casmutils/xtal/structure_tools.hpp"
#include "multishift/slab.hpp"
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/xtal/symmetry.hpp>
#include <cassert>
#include <cmath>
#include <multishift/definitions.hpp>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace
{
using namespace mush;
Eigen::Matrix2d make_2d_column_matrix(const Lattice lat3d)
{
    assert(almost_equal(lat3d.column_vector_matrix()(2, 0), 0.0, 1e-10));
    assert(almost_equal(lat3d.column_vector_matrix()(2, 1), 0.0, 1e-10));

    return lat3d.column_vector_matrix().block<2, 2>(0, 0);

    Eigen::Vector2d a2d(lat3d.a()(0), lat3d.a()(1));
    Eigen::Vector2d b2d(lat3d.b()(0), lat3d.b()(1));

    Eigen::Matrix2d lat_mat2d;
    lat_mat2d.col(0) = a2d;
    lat_mat2d.col(1) = b2d;

    return lat_mat2d;
};

Eigen::Matrix3d make_3d_column_matrix(const Eigen::Matrix2d& lat2d)
{
    Eigen::Matrix3d lat3d = Eigen::Matrix3d::Identity();
    lat3d.block<2, 2>(0, 0) = lat2d;
    return lat3d;
}

double coarse_round(double x)
{
    double trucated = std::round(100000.0 * x) / 100000.0;
    return std::round(trucated);
}
} // namespace

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

std::pair<Eigen::Matrix3d, Eigen::Matrix3d> polar_decomposition(Eigen::Matrix3d const& F)
{
    Eigen::Matrix3d U = Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(F.transpose() * F).operatorSqrt();
    Eigen::Matrix3d R = F * U.inverse();
    return std::make_pair(R, U);
}

namespace frankenstein
{
xtal::Structure stack(const std::vector<xtal::Structure>& sub_strucs)
{
    // TODO:
    // Assert that ab matches

    // Create a new lattice that has the same ab vectors. but summed up
    // the c vectors of every structure
    Eigen::Matrix3d stacked_lat_mat = sub_strucs[0].lattice().column_vector_matrix();
    for (int i = 1; i < sub_strucs.size(); i++)
    {
        Eigen::Matrix3d lat_mat = sub_strucs[i].lattice().column_vector_matrix();
        stacked_lat_mat.col(2) = stacked_lat_mat.col(2) + lat_mat.col(2);
    }

    // We now have a template lattice with the right shape, we'll put the
    // sites inside in a second. It already has the basis for the bottom of
    // the stack
    Lattice stacked_lat(stacked_lat_mat);
    std::vector<Site> stacked_basis = sub_strucs[0].basis_sites();

    // For each structure we stack, we'll take the basis, shift it up by the
    // approprate amount, and stick it into our template stacked structure
    Eigen::Vector3d c_shift = Eigen::Vector3d::Zero();
    for (int i = 1; i < sub_strucs.size(); i++)
    {
        // determine appropriate c-axis shift for position in stacking
        c_shift += sub_strucs[i].lattice().column_vector_matrix().col(2);

        // Shift each site of the basis by the appropriate c shift,
        // and adds them to the stacked structure
        for (const Site& s : sub_strucs[i].basis_sites())
        {
            Coordinate new_coord = Coordinate(s.cart() + c_shift);
            stacked_basis.emplace_back(new_coord, s.label());
        }
    }
    return Structure(stacked_lat, stacked_basis);
}
} // namespace frankenstein
} // namespace xtal
} // namespace casmutils

#include <casmutils/xtal/structure.hpp>
namespace mush
{
/// Returns orhtogonal unit vectors oriented such that the point along the
/// a vector, the ab plane normal, and whatever is perpendicular to that
Eigen::Matrix3d slab_unit_vectors(const Lattice& slab)
{
    if (slab.column_vector_matrix().determinant() < 0)
    {
        throw std::runtime_error("Encountered a left handed lattice when making slab unit vectors.");
    }

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

Eigen::Matrix3d make_alignment_matrix(const Lattice& lat)
{
    Eigen::Matrix3d lat_span_to_standard = slab_unit_vectors(lat).inverse();
    return lat_span_to_standard;
}

Lattice make_aligned_lattice(const Lattice& lat)
{
    Eigen::Matrix3d lat_span_to_standard = make_alignment_matrix(lat);
    Lattice aligned_lat = make_transformed_lattice(lat, lat_span_to_standard);

    assert(almost_equal(aligned_lat.a()(2), 0.0, 1e-10));
    assert(almost_equal(aligned_lat.a()(1), 0.0, 1e-10));
    assert(almost_equal(aligned_lat.b()(2), 0.0, 1e-10));

    assert(almost_equal((aligned_lat.a().cross(aligned_lat.b())).normalized(), Eigen::Vector3d(0, 0, 1), 1e-10));
    assert(almost_equal(lat.column_vector_matrix().determinant(), aligned_lat.column_vector_matrix().determinant(), 1e-10));
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
    const auto [moire_lat, aligned_lat, rot_lat] = make_aligned_moire_lattice(lat, degrees);

    Eigen::Matrix2d moire_lat_2d = ::make_2d_column_matrix(moire_lat);
    Eigen::Matrix2d aligned_lat_2d = ::make_2d_column_matrix(aligned_lat);
    Eigen::Matrix2d rot_lat_2d = ::make_2d_column_matrix(rot_lat);

    Eigen::Matrix2d aligned_to_moire_transform = aligned_lat_2d.inverse() * moire_lat_2d;
    Eigen::Matrix2d aligned_to_moire_transform_round = aligned_to_moire_transform.unaryExpr(&::coarse_round);
    // Assert not nan
    std::cout << "DEGREES: " << degrees << "\n";
    assert(aligned_to_moire_transform == aligned_to_moire_transform);

    Eigen::Matrix2d rot_to_moire_transform = rot_lat_2d.inverse() * moire_lat_2d;
    Eigen::Matrix2d rot_to_moire_transform_round = rot_to_moire_transform.unaryExpr(&::coarse_round);
    assert(rot_to_moire_transform == rot_to_moire_transform);

    Eigen::Matrix2d aligned_superlattice_2d = aligned_lat_2d * aligned_to_moire_transform_round;
    Eigen::Matrix2d rot_superlattice_2d = rot_lat_2d * rot_to_moire_transform_round;

    Eigen::Matrix2d approx_superlattice_2d = (aligned_superlattice_2d + rot_superlattice_2d) / 2.0;
    Eigen::Matrix2d approx_aligned_2d = approx_superlattice_2d * aligned_to_moire_transform_round.inverse();
    Eigen::Matrix2d approx_rot_2d = approx_superlattice_2d * rot_to_moire_transform_round.inverse();

    assert(approx_aligned_2d == approx_aligned_2d);
    assert(approx_rot_2d == approx_rot_2d);

    Eigen::Matrix3d approx_superlattice_col_mat = ::make_3d_column_matrix(approx_superlattice_2d);
    // The moire lattice has zero c vector, but we still want to keep the c vector around for the aligned and rotated lattices
    Eigen::Matrix3d approx_aligned_col_mat = ::make_3d_column_matrix(approx_aligned_2d);
    approx_aligned_col_mat.col(2) = aligned_lat.c();
    Eigen::Matrix3d approx_rot_col_mat = ::make_3d_column_matrix(approx_rot_2d);
    approx_rot_col_mat.col(2) = rot_lat.c();

    return std::make_tuple(approx_superlattice_col_mat, approx_aligned_col_mat, approx_rot_col_mat);
}

MoireApproximant::MoireApproximant(const Lattice& lat, double degrees)
    : input_lattice(lat),
      input_degrees(degrees),
      alignment_rotation(make_alignment_matrix(lat)),
      aligned_lattice(this->make_default_lattice()),
      rotated_lattice(this->make_default_lattice()),
      moire_lattice(this->make_default_lattice()),
      approximate_lattices({this->make_default_lattice(), this->make_default_lattice()})
{
    std::tie(moire_lattice, aligned_lattice, rotated_lattice) = make_aligned_moire_lattice(input_lattice, input_degrees);

    auto [approx_moire, approx_aligned, approx_rot] = make_approximant_moire_lattice(input_lattice, input_degrees);
    this->approximate_lattices[0] = approx_aligned;
    this->approximate_lattices[1] = approx_rot;

    Lattice approx_aligned_moire(approx_moire.a(), approx_moire.b(), aligned_lattice.c());
    Lattice approx_rot_moire(approx_moire.a(), approx_moire.b(), rotated_lattice.c());

    // TODO: Make casmutils compatible
    cu::xtal::Superlattice aligned_to_moire(approx_aligned.__get(), approx_aligned_moire.__get());
    cu::xtal::Superlattice rot_to_moire(approx_rot.__get(), approx_rot_moire.__get());

    this->approximate_moire_integer_transformations[0] = aligned_to_moire.transformation_matrix();
    this->approximate_moire_integer_transformations[1] = rot_to_moire.transformation_matrix();

    this->approximation_deformations[0] = approx_aligned.column_vector_matrix() * aligned_lattice.column_vector_matrix().inverse();
    this->approximation_deformations[1] = approx_rot.column_vector_matrix() * rotated_lattice.column_vector_matrix().inverse();
}

Lattice make_prismatic_lattice(const Lattice& lat)
{
    Eigen::Vector3d orthogonal_unit = lat.a().cross(lat.b()).normalized();
    Eigen::Vector3d new_c = lat.c().dot(orthogonal_unit) * orthogonal_unit;
    return Lattice(lat.a(), lat.b(), new_c);
}

MoirePrismaticApproximant::MoirePrismaticApproximant(const Lattice& lat, double degrees)
    : MoireApproximant(make_prismatic_lattice(lat), degrees)
{
}

ReducedAngleMoirePrismaticApproximant::ReducedAngleMoirePrismaticApproximant(const Lattice& lat, double degrees)
    :
        original_degrees(degrees),
        prismatic_aligned_lattice(make_prismatic_lattice(make_aligned_lattice(lat))),
        prismatic_rotated_lattice(make_twisted_lattice(prismatic_aligned_lattice,degrees)),
      reduced_angle_operation(
          this->find_reduced_angle_operation(prismatic_aligned_lattice, prismatic_rotated_lattice)),
      reduced_angle(
          calculate_rotation_angle(prismatic_aligned_lattice,
                                   make_transformed_lattice(prismatic_rotated_lattice, reduced_angle_operation.matrix))),
      reduced_angle_moire(lat, reduced_angle)
{
}

ReducedAngleMoirePrismaticApproximant::CartOp
ReducedAngleMoirePrismaticApproximant::find_reduced_angle_operation(const Lattice& aligned_lat, const Lattice& rotated_lat) const
{
    auto in_plane_ops = in_plane_rotation_point_operations(rotated_lat);

    const CartOp* best_op;
    double smallest_angle = 100000;

    for (const CartOp& op : in_plane_ops)
    {
        Lattice rerotated_lat = make_transformed_lattice(rotated_lat, op.matrix);
        double candidate_angle = calculate_rotation_angle(aligned_lat, rerotated_lat);
        if (std::abs(candidate_angle) < std::abs(smallest_angle))
        {
            smallest_angle = candidate_angle;
            best_op = &op;
        }
    }

    return *best_op;
}

double ReducedAngleMoirePrismaticApproximant::calculate_rotation_angle(const Lattice& aligned_lat, const Lattice& rotated_lat) const
{
    Eigen::MatrixXd rotation = rotated_lat.column_vector_matrix() * aligned_lat.column_vector_matrix().inverse();

    // This stuff better be aligned along all the right directions:
    assert(almost_equal(rotation(0, 2), 0.0, 1e-8));
    assert(almost_equal(rotation(1, 2), 0.0, 1e-8));
    assert(almost_equal(rotation(2, 0), 0.0, 1e-8));
    assert(almost_equal(rotation(2, 1), 0.0, 1e-8));
    assert(almost_equal(rotation(2, 2), 1.0, 1e-8));

    return std::atan2(rotation(1, 0), rotation(0, 0)) * 180 / M_PI;
}

std::vector<ReducedAngleMoirePrismaticApproximant::CartOp>
ReducedAngleMoirePrismaticApproximant::in_plane_rotation_point_operations(const Lattice& lat) const
{
    // TODO: Allow setting tol?
    auto point_group = cu::xtal::make_point_group(lat, 1e-5);
    std::vector<CartOp> point_group_subset;
    for (const CartOp& op : point_group)
    {
        // Discard anything improper
        if (op.matrix.determinant() < 0.5)
        {
            continue;
        }

        // Discard anything that alters the c vector (you're using prismatic lattices remember?)
        Eigen::Vector3d rotated_c = op.matrix * lat.c();
        if (!almost_equal(lat.c(), rotated_c))
        {
            continue;
        }

        point_group_subset.push_back(op);
    }

    return point_group_subset;
}
} // namespace mush
