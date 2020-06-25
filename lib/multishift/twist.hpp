#ifndef TWIST_HH
#define TWIST_HH

#include <multishift/definitions.hpp>
#include <casmutils/xtal/lattice.hpp>
#include <tuple>
#include <array>

//TODO: Put in casm utilities
#include <casm/crystallography/Superlattice.hh>

namespace casmutils
{
    namespace xtal
    {
        using CASM::xtal::Superlattice;

        /// Retrun the reciprocal of the given lattice, where each reciprocal vector
        /// is perpendicular to the other two real ones
        Lattice make_reciprocal(const Lattice& real_lattice);

        /// For a deformation matrix F, return the polar decomposition of its rotational and
        /// strain components R (first) and U (second), where F=R*U
        std::pair<Eigen::Matrix3d,Eigen::Matrix3d> polar_decomposition(Eigen::Matrix3d const &F);
    }
}

namespace mush
{
    namespace cu = casmutils;
    
    ///Create a rotation matrix that rotates the given lattice by the specified angle
    ///within the plane spanned by the ab vectors of the lattice.
    ///The rotation should be applied to column vectors.
    Eigen::Matrix3d make_twist_rotation_matrix(const Lattice& lat, double degrees);

    ///Create a new lattice that has been rotated about the ab normal by the specified angle
    Lattice make_twisted_lattice(const Lattice& lat, double degrees);

    ///Creates a lattice whose lattice points follow the Moire pattern that results from
    ///twisting the given lattice by the specified angle. 
    ///Because the lattice will be rotated in the ab plane, the resulting Moire lattice
    ///will be two dimensional. The c vector of the returned lattice is therefore meaningless,
    ///and the a and b vectors will be aligned along the xy plane
    ///The second and thrid returned lattices are the original lattice and rotated
    ///lattice after alignment along the ab-plane.
    std::tuple<Lattice,Lattice,Lattice> make_aligned_moire_lattice(const Lattice& lat, double degrees);

    ///Returns the same lattice, but rotated such that the a vector points along the
    ///Cartesian x direction, and the b vector is parallel to the xy plane.
    Lattice make_aligned_lattice(const Lattice& lat);

    ///Returns Moire, plus two lattices that are the same under the rotation of the specified degrees.
    ///The lattices will differ from the input lattice by a small deformation that allows
    ///the construction of a commensurate Moire supercell that can fit both the original
    ///and rotated orientations.
    ///Returns approximate moire, approximate aligned, and apprpximate rotated
    std::tuple<Lattice,Lattice,Lattice> make_approximant_moire_lattice(const Lattice& lat, double degrees);
    
    ///representation of a Moire lattice, given an input lattice and rotation angle, as well
    ///as approximant versions of the same Moire lattice, where some strain has been introduced
    ///to enforce periodicity.
    struct MoireApproximant
    {
        typedef Eigen::Matrix3l matrix_type;

        MoireApproximant(const Lattice& lat, double degrees);

        ///Lattice given at construction
        Lattice input_lattice;

        ///Degrees given at construction
        double input_degrees;

        ///Rotation matrix to go from the input lattice to the alinged lattice
        Eigen::Matrix3d alignment_rotation;

        ///Same as the input lattice, but rotated such that the a and b vectors span the x-y plane (no z component)
        Lattice aligned_lattice;

        ///Same as the aligned lattice, but rotated by the specified degrees.
        Lattice rotated_lattice;

        ///The true Moire lattice resulting from the superposition of the input lattice and it's rotated version.
        ///The pattern resulting from this Moire lattice is only semi-periodic, and not useful for DFT calculations.
        Lattice moire_lattice;

        ///The aligned and rotated lattices with some strain introduced, such creating superlattices
        ///from them results in fully periodic Moire lattices
        std::array<Lattice,2> approximate_lattices;

        ///Integer transformation matrices that convert the approximate lattices into the Moire lattice.
        std::array<matrix_type,2> approximate_moire_integer_transformations;

        ///Deformation introduced by making the approximations to introduce complete periodicity
        std::array<Eigen::Matrix3d,2> approximation_deformations;

        private:

        Lattice make_default_lattice()
        {
            return cu::xtal::Lattice(Eigen::Matrix3d::Zero());
        }
    };
}

#endif
