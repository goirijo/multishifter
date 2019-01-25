#include "casmutils/structure.hpp"

/// This is some stuff you extracted out of an ancient casm version
namespace legacy
{
///return signed angle, in radians, between -pi and pi that describe separation in direction of two vectors
double signed_angle_between(const Eigen::Vector3d& u, const Eigen::Vector3d& v)
{
}

/**
 * Given a rotational matrix rotate the
 * ENTIRE structure in cartesian space. This yields the same
 * structure, just with different cartesian definitions of
 * the lattice vectors
 */

void reorient_structure(Rewrap::Structure* reorientstruc, const Eigen::Matrix3d reorientmat)
{

    auto reorientlatvecs = reorientmat * (reorientstruc->lattice().lat_column_mat());

    CASM::Lattice reorientlat(reorientlatvecs);
    reorientstruc->set_lattice(reorientlat, CASM::FRAC);

    return;
}

/**
 *	This function takes two structures and modifies one to
 *	have its lattice vector realigned such that they match
 *	the reference structure. The modified structure ends up
 *	having a and axb pointing in the same direction as the
 *	reference structure.
 */

void align_with(Rewrap::Structure* alignable_struc, const Rewrap::Structure& refstruc, bool override)
{

    Eigen::Vector3d rotaxis;
    double angle;
    Eigen::Matrix3d rotmat;

    // First align the a vectors
    rotaxis = alignable_struc->lattice.get(0).cross(refstruc.lattice.get(0));
    angle = alignable_struc->lattice.get(0).get_signed_angle(refstruc.lattice.get(0), rotaxis);

    rotmat = rotaxis.get_rotation_mat(angle);

    /* Rewrap::Structure aaligned = reorient(rotmat, override); */
    reorient_structure(alignable_struc);

    // Then turn axb along the a axis
    rotaxis = alignable_struc->lattice.get(0);
    angle = (alignable_struc->lattice.get(0).cross(alignable_struc->lattice.get(1)))
                .get_signed_angle(refstruc.lattice.get(0).cross(refstruc.lattice.get(1)), rotaxis);

    rotmat = rotaxis.get_rotation_mat(angle);
    return;
}

//***********************************************************
/**
 * Reorients structure to have a along (x,0,0), b along (x,y,0)
 * and c along (x,y,z).
 */
//***********************************************************

Rewrap::Structure Rewrap::Structure::align_standard(bool override) const
{
    if (!is_primitive() && !override)
    {
        std::cerr << "ERROR in Rewrap::Structure::align_with_standard" << std::endl;
        std::cerr << "This function is for primitive cells only. Reduce your structure and try again." << std::endl;
        exit(15);
    }

    Vector3<Vector3<double>> standardvecs;
    standardvecs.at(0) = Vector3<double>(1, 0, 0);
    standardvecs.at(1) = Vector3<double>(0, 1, 0);
    standardvecs.at(2) = Vector3<double>(0, 0, 1);

    Lattice standardlat(standardvecs);
    Rewrap::Structure standardstruc(standardlat);

    return align_with(standardstruc, override);
}
} // namespace legacy
