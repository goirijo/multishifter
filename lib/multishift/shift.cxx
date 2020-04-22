#include <casmutils/xtal/structure.hpp>
#include "./shift.hpp"
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/definitions.hpp>

namespace casmutils
{
    namespace xtal
    {
    }
}

namespace mush
{
    std::vector<cu::xtal::Structure> make_cleaved_structures(const cu::xtal::Structure& slab, const std::vector<double>& cleavage_values)
    {
        //TODO: This could go astray if you get a left handed lattice
        std::vector<cu::xtal::Structure> cleaved_structures;
        
        Eigen::Vector3d unit_normal=slab.lattice().a().cross(slab.lattice().b()).normalized();
        for(double cleave : cleavage_values)
        {
            Eigen::Vector3d new_c_vector=slab.lattice().c()+cleave*unit_normal;
            cu::xtal::Lattice cleaved_lattice(slab.lattice().a(),slab.lattice().b(),new_c_vector);

            cleaved_structures.emplace_back(slab.set_lattice(cleaved_lattice,cu::xtal::CART));
        }
        return cleaved_structures;
    }
}
