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
}

#endif
