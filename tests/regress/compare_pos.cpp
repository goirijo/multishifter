#include <casmutils/xtal/structure_tools.hpp>
#include <casmutils/xtal/lattice.hpp>
#include <filesystem>
#include <iostream>

//Super rudimentary way to check if two POSCAR files are the same,
//without having to worry about indentation, significant figures,
//or +0.0 vs -0.0
int main(int argc, char** argv)
{
    namespace fs=std::filesystem;
    namespace cu=casmutils;

    fs::path struc0_path=argv[1];
    fs::path struc1_path=argv[2];

    auto struc0=cu::xtal::Structure::from_poscar(struc0_path);
    auto struc1=cu::xtal::Structure::from_poscar(struc1_path);

    if(!struc0.lattice().column_vector_matrix().isApprox(struc1.lattice().column_vector_matrix()))
    {
        return 1;
    }

    if(struc0.basis_sites().size()!=struc1.basis_sites().size())
    {
        return 2;
    }

    int num_sites=struc0.basis_sites().size();

    for(int i=0; i<num_sites; ++i)
    {
        const auto& s0=struc0.basis_sites().at(i);
        const auto& s1=struc1.basis_sites().at(i);

        if(!s0.cart().isApprox(s1.cart()))
        {
            return 3;
        }

        if(s0.label()!=s1.label())
        {
            return 4;
        }
    }
    
    return 0;
}
