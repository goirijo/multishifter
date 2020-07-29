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

bool is_within_voronoi(const Eigen::Vector3d& v, const cu::xtal::Lattice& lat)
{
    cu::xtal::Coordinate vw(v);
    vw.bring_within_wigner_seitz(lat);
    return almost_equal(v,vw.cart(),1e-13);
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
    Eigen::Matrix3d R = F*U.inverse();
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

MoireLattice::MoireLattice(const Lattice& lat, double degrees)
    : input_lattice(lat),
      input_degrees(degrees),
      aligned_lattice(make_prismatic_lattice(make_aligned_lattice(input_lattice))),
      reciprocal_aligned_lattice(cu::xtal::make_reciprocal(aligned_lattice)),
      rotated_lattice(make_twisted_lattice(aligned_lattice, input_degrees)),
      reciprocal_rotated_lattice(cu::xtal::make_reciprocal(rotated_lattice)),
      full_reciprocal_difference(this->calculate_reciprocal_difference()),
      aligned_brillouin_zone_reciprocal_difference(this->bring_vectors_into_voronoi(full_reciprocal_difference, reciprocal_aligned_lattice)),
      rotated_brillouin_zone_reciprocal_difference(this->bring_vectors_into_voronoi(full_reciprocal_difference, reciprocal_rotated_lattice)),
      aligned_moire_lattice(make_moire_lattice_from_reciprocal_difference(aligned_brillouin_zone_reciprocal_difference, aligned_lattice.c())),
      rotated_moire_lattice(make_moire_lattice_from_reciprocal_difference(rotated_brillouin_zone_reciprocal_difference, aligned_lattice.c())),
      moire_lattice(nullptr)
{
    //calculate overlap
    for(int i=0; i<2; ++i)
    {
        this->is_within_brillouin_zone_overlap[&aligned_moire_lattice][i]=this->is_within_voronoi(aligned_brillouin_zone_reciprocal_difference.col(i),reciprocal_rotated_lattice);
        this->is_within_brillouin_zone_overlap[&rotated_moire_lattice][i]=this->is_within_voronoi(rotated_brillouin_zone_reciprocal_difference.col(i),reciprocal_aligned_lattice);
    }
    
    //Point to appropriate moire lattice
    for(const Lattice* moire_ptr : {&aligned_moire_lattice, &rotated_moire_lattice})
    {
        if(is_within_brillouin_zone_overlap[moire_ptr][0] && is_within_brillouin_zone_overlap[moire_ptr][1])
        {
            moire_lattice=moire_ptr;
            break;
        }
    }
}

bool MoireLattice::is_within_voronoi(const Eigen::Vector2d& v, const cu::xtal::Lattice& lat) const
{
    Eigen::Vector3d v3(v(0),v(1),0);
    return ::is_within_voronoi(v3,lat);
}

Eigen::Matrix2d MoireLattice::calculate_reciprocal_difference() const
{
    Eigen::Matrix3d diff =
        this->reciprocal_rotated_lattice.column_vector_matrix() - this->reciprocal_aligned_lattice.column_vector_matrix();

    Eigen::Vector3d z3=Eigen::Vector3d::Zero();

    // There should be no difference in the c vector
    assert(almost_equal(diff.col(2), z3));
    // There should be no z components for anything
    assert(almost_equal(diff.row(2), z3.transpose()));

    return diff.block<2, 2>(0, 0);
}

Eigen::Matrix2d MoireLattice::bring_vectors_into_voronoi(const Eigen::Matrix2d& col_vectors, const Lattice& lat)
{
    // You must provide a prismatic aligned lattice for this to word
    if (!almost_equal(lat.a()(2), 0.0) && !almost_equal(lat.b()(2), 0.0) && !almost_equal(lat.c()(0), 0.0) &&
        !almost_equal(lat.c()(1), 0.0))
    {
        throw std::runtime_error("The provided lattice must be prismatic (perpendicular c vector) and aligned in the xy plane");
    }

    //TODO: Bring within voronoi. Function is missing in cu
    cu::xtal::Coordinate ka(col_vectors(0,0),col_vectors(1,0),0);
    cu::xtal::Coordinate kb(col_vectors(0,1),col_vectors(1,1),0);

    ka.bring_within_wigner_seitz(lat);
    kb.bring_within_wigner_seitz(lat);

    Eigen::Matrix2d col_vectors_within;
    col_vectors_within.col(0)=ka.cart().head(2);
    col_vectors_within.col(1)=kb.cart().head(2);

    return col_vectors_within;
}

Lattice MoireLattice::make_moire_lattice_from_reciprocal_difference(const Eigen::Matrix2d diff, const Eigen::Vector3d& real_c_vector)
{
    Eigen::Matrix3d phony_recip_mat=Eigen::Matrix3d::Identity();
    phony_recip_mat.block<2,2>(0,0)=diff;
    Lattice phony_recip_lat(phony_recip_mat.col(0),phony_recip_mat.col(1),phony_recip_mat.col(2));

    Lattice phony_real_moire_lattice=cu::xtal::make_reciprocal(phony_recip_lat);
    Eigen::Matrix3d final_real_moire_mat=phony_real_moire_lattice.column_vector_matrix();
    final_real_moire_mat.col(2)=real_c_vector;

    return Lattice(final_real_moire_mat.col(0),final_real_moire_mat.col(1),final_real_moire_mat.col(2));
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

std::tuple<Lattice, Lattice, Lattice> make_approximant_moire_lattice(const Lattice& moire_lat, const Lattice& aligned_lat, const Lattice& rot_lat)
{
    Eigen::Matrix2d moire_lat_2d = ::make_2d_column_matrix(moire_lat);
    Eigen::Matrix2d aligned_lat_2d = ::make_2d_column_matrix(aligned_lat);
    Eigen::Matrix2d rot_lat_2d = ::make_2d_column_matrix(rot_lat);

    Eigen::Matrix2d aligned_to_moire_transform = aligned_lat_2d.inverse() * moire_lat_2d;
    Eigen::Matrix2d aligned_to_moire_transform_round = aligned_to_moire_transform.unaryExpr(&::coarse_round);
    // Assert not nan
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

MoireApproximant::MoireApproximant(const Lattice& moire_lat, const Lattice& aligned_lat, const Lattice& rotated_lat):
    approximate_moire_lattice(this->default_lattice())
{
    auto [approx_moire, approx_aligned, approx_rot] = make_approximant_moire_lattice(moire_lat,aligned_lat,rotated_lat);
    approximate_moire_lattice=approx_moire;

    this->approximate_lattices.emplace(&aligned_lat, approx_aligned);
    this->approximate_lattices.emplace(&rotated_lat, approx_rot);

    Lattice approx_aligned_moire(approx_moire.a(), approx_moire.b(), aligned_lat.c());
    Lattice approx_rot_moire(approx_moire.a(), approx_moire.b(), rotated_lat.c());

    // TODO: Make casmutils compatible
    cu::xtal::Superlattice aligned_to_moire(approx_aligned.__get(), approx_aligned_moire.__get());
    cu::xtal::Superlattice rot_to_moire(approx_rot.__get(), approx_rot_moire.__get());

    this->approximate_moire_integer_transformations[&aligned_lat] = aligned_to_moire.transformation_matrix();
    this->approximate_moire_integer_transformations[&rotated_lat] = rot_to_moire.transformation_matrix();

    this->approximation_deformations[&aligned_lat] = approx_aligned.column_vector_matrix() * aligned_lat.column_vector_matrix().inverse();
    this->approximation_deformations[&rotated_lat] = approx_rot.column_vector_matrix() * rotated_lat.column_vector_matrix().inverse();
}

Lattice make_prismatic_lattice(const Lattice& lat)
{
    Eigen::Vector3d orthogonal_unit = lat.a().cross(lat.b()).normalized();
    Eigen::Vector3d new_c = lat.c().dot(orthogonal_unit) * orthogonal_unit;
    return Lattice(lat.a(), lat.b(), new_c);
}

MoireGenerator::MoireGenerator(const Lattice& input_lat, double degrees):
    moire(input_lat, degrees),
    aligned_moire_approximant(moire.aligned_moire_lattice,moire.aligned_lattice, moire.rotated_lattice),
    rotated_moire_approximant(moire.rotated_moire_lattice,moire.aligned_lattice,moire.rotated_lattice),
    aligned_key(&moire.aligned_lattice),
    rotated_key(&moire.rotated_lattice)
{
}

} // namespace mush
