#include "./slice.hpp"
#include <filesystem>
#include <casmutils/xtal/structure_tools.hpp>

void write_slicer_structures(const mush::Slicer& slicer, const mush::fs::path& slices_path, std::ostream& log)
{
        log<<"Back up prim used...\n";
        mush::cu::xtal::write_poscar(slicer.prim,slices_path/"prim.vasp");

        log<<"Write sliced prim...\n";
        mush::cu::xtal::write_poscar(slicer.sliced_prim,slices_path/"sliced_prim.vasp");

        log<<"Write all possible basis translations...";
        for(int i=0; i<slicer.floored_sliced_prims.size(); ++i)
        {
            log<<"..";
            mush::fs::path target=slices_path/("sliced_prim_floor."+std::to_string(i)+".vasp");
            mush::cu::xtal::write_poscar(slicer.floored_sliced_prims[i],target);
        }
        log<<std::endl;

        return;
}

