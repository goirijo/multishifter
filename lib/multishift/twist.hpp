#ifndef TWIST_HH
#define TWIST_HH

#include <casmutils/xtal/lattice.hpp>

namespace casmutils
{
    namespace xtal
    {
        /// Retrun the reciprocal of the given lattice, where each reciprocal vector
        /// is perpendicular to the other two real ones
Lattice make_reciprocal(const Lattice& real_lattice);

    }
}

namespace mush
{
    namespace cu = casmutils;
    
    ///Create a rotation matrix that rotates the given lattice by the specified angle
    ///within the plane spanneed by the ab vectors of the lattice.
    ///The rotation should be applied to column vectors.
    Eigen::Matrix3d make_twist_rotation_matrix(const cu::xtal::Lattice& lat, double degrees);

    ///Create a new lattice that has been rotated about the ab normal by the specified angle
    cu::xtal::Lattice make_twisted_lattice(const cu::xtal::Lattice& lat, double degrees);

    ///Creates a lattice whose lattice points follow the Moire pattern that results from
    ///twisting the given lattice by the specified angle. 
    ///Because the lattice will be rotated in the ab plane, the resulting Moire lattice
    ///will be two dimensional. The c vector of the returned lattice is therefore meaningless,
    ///and the a and b vectors will be aligned along the xy plane
    cu::xtal::Lattice make_aligned_moire_lattice(const cu::xtal::Lattice& lat, double degrees);

    ///Returns the same lattice, but rotated such that the a vector points along the
    ///Cartesian x direction, and the b vector is parallel to the xy plane.
    cu::xtal::Lattice make_aligned_lattice(const cu::xtal::Lattice& lat);
    
}

#endif
