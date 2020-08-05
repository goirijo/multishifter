#ifndef TWIST_HH
#define TWIST_HH

#include "casmutils/xtal/structure.hpp"
#include <array>
#include <casmutils/sym/cartesian.hpp>
#include <casmutils/xtal/lattice.hpp>
#include <multishift/definitions.hpp>
#include <tuple>
#include <unordered_map>

// TODO: Put in casm utilities
#include <casm/crystallography/Superlattice.hh>
#include <vector>

namespace mush
{
namespace cu = casmutils;

/// Create a rotation matrix that rotates the given lattice by the specified angle
/// within the plane spanned by the ab vectors of the lattice.
/// The rotation should be applied to column vectors.
Eigen::Matrix3d make_twist_rotation_matrix(const Lattice& lat, double degrees);

/// Create a new lattice that has been rotated about the ab normal by the specified angle
Lattice make_twisted_lattice(const Lattice& lat, double degrees);

/// Returns the same lattice, but rotated such that the a vector points along the
/// Cartesian x direction, and the b vector is parallel to the xy plane.
Lattice make_aligned_lattice(const Lattice& lat);

/// Returns the same lattice, but the c vector has been modified to be orthogonal to the
/// ab vectors. This will break periodicity, but not the thickness of the slab.
Lattice make_prismatic_lattice(const Lattice& lat);

/// Constructs Moire lattice for the given angle, and keeps information of the steps made along
/// the way. Before anything happens, the input lattice will be transformed to an aligned
/// prismatic one, breaking periodicity along the c axis.
/// The Moire lattice is likely not coincident with the smaller lattices. For constructing
/// periodic interference patterns, a small deformation must be introduced, accomplished
/// using ApproximantMoireLattice
struct MoireLattice
{
    enum class LATTICE
    {
        ALIGNED,
        ROTATED
    };

    using ZONE = LATTICE;

    MoireLattice(const Lattice& lat, double degrees);

    /// The lattice given at construction (not necessarily aligned)
    Lattice input_lattice;

    /// The degrees given at construction
    double input_degrees;

    /// The input lattice, after being transfromed to be prismatic and properly aligned along the xy plane
    Lattice aligned_lattice;
    Lattice reciprocal_aligned_lattice;

    /// The aligned lattice after being rotated by the input degrees
    Lattice rotated_lattice;
    Lattice reciprocal_rotated_lattice;

    /// Returns the real lattice of either ALIGNED or ROTATED
    const Lattice& real(LATTICE lat) const { return lat == LATTICE::ALIGNED ? aligned_lattice : rotated_lattice; } 

    /// Returns the reciprocal lattice of either ALIGNED or ROTATED
     const Lattice& reciprocal(LATTICE lat) const { return lat == LATTICE::ALIGNED ? reciprocal_aligned_lattice : reciprocal_rotated_lattice; }

    /// Result from subtracting the reciprocal aligned and rotated lattices.
    /// There is no c vector difference, it's zero by construction.
    Eigen::Matrix2d full_reciprocal_difference;

    /// The reciprocal difference, brought into the first Brillouin zone of the
    /// aligned lattice
    Eigen::Matrix2d aligned_brillouin_zone_reciprocal_difference;

    /// The reciprocal difference, brought into the first Brillouin zone of the
    /// rotated lattice
    Eigen::Matrix2d rotated_brillouin_zone_reciprocal_difference;

    /// Moire lattice constructed from the reciprocal difference brought within the
    /// aligned lattice brillouin zone (this is equivalent to using the rotated
    /// lattice brillouin zone, but results in different lattice vectors)
    Lattice aligned_moire_lattice;

    /// Moire lattice constructed from the reciprocal difference brought within the
    /// rotated lattice brillouin zone (this is equivalent to using the aligned
    /// lattice brillouin zone, but results in different lattice vectors)
    Lattice rotated_moire_lattice;
    
    /// Returns the moire lattice of either ALIGNED or ROTATED Brillouin zone
     const Lattice& moire(LATTICE lat) const { return lat == LATTICE::ALIGNED ? aligned_moire_lattice : rotated_moire_lattice; }

    /// Maps the address of a Moire lattice to an array that specifies if its reciprocal vectors
    /// fall within the Brillouin zone of the other (rotated/aligned) lattice. For example
    /// brillouin_zone_overlap[LATTICE::ALIGNED][0] returns true if the reciprocal <a> vector
    /// of the aligned moire lattice falls within the first Brillouin zone of the rotated
    /// moire lattice
    std::unordered_map<LATTICE, std::array<bool, 2>> is_within_brillouin_zone_overlap;

    /// Brings the given vectors (columns in matrix) into the first voronoi (Wigner Seitz) cell
    /// of the provided lattice. Note that the function assumes vectors are in the xy plane
    ///(no z component) and that the lattice used for the Brillouin zone is prismatic
    /// and aligned.
    static Eigen::Matrix2d bring_vectors_into_voronoi(const Eigen::Matrix2d& col_vectors, const Lattice& lat);

    /// Given the a and b reciprocal Moire lattice vectors, transform them into a real
    /// Moire lattice, assigning the c vector as the last provided parameter.
    /// Vectors must have no z component.
    static Lattice make_moire_lattice_from_reciprocal_difference(const Eigen::Matrix2d diff, const Eigen::Vector3d& real_c_vector);

private:
    Eigen::Matrix2d calculate_reciprocal_difference() const;
    bool is_within_voronoi(const Eigen::Vector2d& v, const cu::xtal::Lattice& lat) const;
};

/// Helper struct to convert an aligned and rotated lattice to a superlattice that's as close
/// as possible to the Moire lattice. The superlattices are calculated by finding the closest
/// integer transformation for each of the aligned and rotated lattices, then applying a deformation
/// to make them coincident. This class makes no checks whatsoever on the input parameters, and assumes
/// that you're giving something sensible that came from the MoireLattice class.
/// Members include unordered maps that use Lattice pointers as keys, the expected pointers are the
/// addresses of the aligned and rotated lattices given at construction.
struct MoireApproximant
{
    using LATTICE = MoireLattice::LATTICE;

    typedef Eigen::Matrix3l matrix_type;

    MoireApproximant(const Lattice& moire_lat, const Lattice& aligned_lat, const Lattice& rotated_lat);

    /// The moire lattice after straining it a bit to make the aligned and rotated lattices coincident
    Lattice approximate_moire_lattice;

    /// The aligned and rotated lattices with some strain introduced, such creating superlattices
    /// from them results in fully periodic Moire lattices
    std::unordered_map<LATTICE, Lattice> approximate_lattices;

    /// Integer transformation matrices that convert the approximate lattices into the Moire lattice.
    std::unordered_map<LATTICE, matrix_type> approximate_moire_integer_transformations;

    /// Deformation introduced by making the approximations to introduce complete periodicity
    std::unordered_map<LATTICE, Eigen::Matrix3d> approximation_deformations;

private:
    Lattice default_lattice() { return Lattice(Eigen::Matrix3d::Zero()); }
};

/// Interface class to access the Moire lattice, which can be relative do different Brillouin zones,
/// as well as relevant deformations on each lattice. Basically just a wrapper class for everythin
/// in MoireLattice and MoireApproximant
class MoireGenerator
{
public:
    using LATTICE = MoireLattice::LATTICE;
    using ZONE = MoireLattice::ZONE;

    MoireLattice moire;

    /// Approximations made using the Moire lattice generated using the fixed (aligned)
    /// Brillouin zone
    MoireApproximant aligned_moire_approximant;

    /// Approximations made using the Moire lattice generated using the rotated
    /// Brillouin zone
    MoireApproximant rotated_moire_approximant;

private:
    /* const Lattice* requested_key(LATTICE lat) { return lat == LATTICE::ALIGNED ? aligned_key : rotated_key; } */

    const MoireApproximant& requested_zone(ZONE brillouin) const
    {
        return brillouin == ZONE::ALIGNED ? aligned_moire_approximant : rotated_moire_approximant;
    }

public:
    MoireGenerator(const Lattice& input_lat, double degrees);

    const Lattice& approximate_lattice(ZONE brillouin, LATTICE lat) const { return requested_zone(brillouin).approximate_lattices.at(lat); }

     const MoireApproximant::matrix_type& approximate_moire_integer_transformation(ZONE bz, LATTICE lat) const
    {
        return requested_zone(bz).approximate_moire_integer_transformations.at(lat);
    }

    const Eigen::Matrix3d& approximation_deformation(ZONE bz, LATTICE lat) { return requested_zone(bz).approximation_deformations.at(lat); }

    double degrees() const {return moire.input_degrees;} 
};

/// Generates slab superstructures that can be stacked together to create bilayers with Moire
/// patterns
class MoireStructureGenerator : MoireGenerator
{
public:
    using ZONE = MoireGenerator::ZONE;
    using LATTICE = MoireGenerator::LATTICE;
    using Structure = cu::xtal::Structure;
    using MoireGenerator::degrees;

    MoireStructureGenerator(const Structure& slab_unit, double degrees);

    Structure layer(ZONE brillouin, LATTICE lat) const;

private:
    const Structure slab_unit;
};
} // namespace mush

#endif
