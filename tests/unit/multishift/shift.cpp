#include "../../autotools.hh"
#include "casmutils/xtal/coordinate.hpp"
#include "casmutils/xtal/structure.hpp"
#include "casmutils/xtal/structure_tools.hpp"

#include <algorithm>
#include <gtest/gtest.h>
#include <multishift/shift.hpp>
#include <casmutils/mapping/structure_mapping.hpp>
#include <vector>

namespace cu = casmutils;

class B2ModifiedC_Vector : public testing::Test
{
protected:
    using Structure = cu::xtal::Structure;
    std::unique_ptr<Structure> b2_ptr;

    virtual void SetUp() override { b2_ptr.reset(new Structure(Structure::from_poscar(mush::autotools::input_filesdir / "b2.vasp"))); }

    void compare_basis_to_b2(const Structure& mutated_b2)
    {
        EXPECT_TRUE(mutated_b2.basis_sites().size() == b2_ptr->basis_sites().size());

        for (int i = 0; i < mutated_b2.basis_sites().size(); ++i)
        {
            const cu::xtal::Site& cl_b = mutated_b2.basis_sites()[i];
            const cu::xtal::Site& b = b2_ptr->basis_sites()[i];
            EXPECT_EQ(cl_b.cart(), b.cart());
            EXPECT_EQ(cl_b.label(), b.label());
        }
    }
};

class CleavingTest : public B2ModifiedC_Vector
{
protected:
    void SetUp() override
    {
        B2ModifiedC_Vector::SetUp();
        cleaved_structures = mush::make_cleaved_structures(*b2_ptr, cleavage_values);
    }

    std::vector<double> cleavage_values{-0.1, 0, 5};
    std::vector<Structure> cleaved_structures;
};

TEST_F(CleavingTest, Consistent_AB_Vectors)
{
    for (const Structure& cleaved_struc : cleaved_structures)
    {
        EXPECT_EQ(b2_ptr->lattice().a(), cleaved_struc.lattice().a());
        EXPECT_EQ(b2_ptr->lattice().b(), cleaved_struc.lattice().b());
    }
}

TEST_F(CleavingTest, C_VectorLengths)
{
    for (int i = 0; i < cleavage_values.size(); ++i)
    {
        if (cleavage_values[i] < 0)
        {
            EXPECT_TRUE(b2_ptr->lattice().c().norm() > cleaved_structures[i].lattice().c().norm());
        }
        else if (cleavage_values[i] == 0)
        {
            EXPECT_EQ(b2_ptr->lattice().c().norm(), cleaved_structures[i].lattice().c().norm());
        }
        else
        {
            EXPECT_TRUE(b2_ptr->lattice().c().norm() < cleaved_structures[i].lattice().c().norm());
        }
    }
}

TEST_F(CleavingTest, SitesHaventMoved)
{
    for (const Structure& cleaved_struc : cleaved_structures)
    {
        compare_basis_to_b2(cleaved_struc);
    }
}

class ShiftingTest : public B2ModifiedC_Vector
{
protected:
    void SetUp() override
    {
        B2ModifiedC_Vector::SetUp();

        std::tie(shift_values, shift_records) = mush::make_uniform_in_plane_shift_vectors(b2_ptr->lattice(), as, bs);
        std::tie(wigner_seitz_shift_values, wigner_seitz_shift_records) =
            mush::make_uniform_in_plane_wigner_seitz_shift_vectors(b2_ptr->lattice(), as, bs);

        shifted_structures = mush::make_shifted_structures(*b2_ptr, shift_values);
        wigner_seitz_shifted_structures = mush::make_shifted_structures(*b2_ptr, wigner_seitz_shift_values);
    }

    int as = 8;
    int bs = 5;
    std::vector<Eigen::Vector3d> shift_values;
    std::vector<mush::ShiftRecord> shift_records;

    std::vector<Eigen::Vector3d> wigner_seitz_shift_values;
    std::vector<mush::ShiftRecord> wigner_seitz_shift_records;

    std::vector<Structure> shifted_structures;
    std::vector<Structure> wigner_seitz_shifted_structures;
};

TEST_F(ShiftingTest, ShiftRecords) { EXPECT_EQ(shift_records, wigner_seitz_shift_records); }

TEST_F(ShiftingTest, ShiftValues)
{
    EXPECT_EQ(shift_values.size(), as * bs);

    // Shift at half along a has record a=5, b=0
    Eigen::Vector3d expected_shift = b2_ptr->lattice().a() * 0.5;
    auto is_half_a_record = [this](const mush::ShiftRecord& record) { return record.a == as/2 && record.b == 0; };
    mush::ShiftRecord half_a_record = *(std::find_if(shift_records.begin(), shift_records.end(), is_half_a_record));

    EXPECT_EQ(expected_shift, shift_values[half_a_record.index]);
    EXPECT_EQ(Eigen::Vector3d::Zero(), shift_values[0]);
}

TEST_F(ShiftingTest, WignerSeitzShiftMatches)
{
    EXPECT_EQ(shift_values.size(),wigner_seitz_shift_values.size());
    int didnt_match=0;
    for(int i=0; i<shift_values.size(); ++i)
    {
        cu::xtal::Coordinate ws_shift(wigner_seitz_shift_values[i]);
        if(!shift_values[i].isApprox(ws_shift.cart()))
        {
            ++didnt_match;
        }

        ws_shift.bring_within(b2_ptr->lattice());
        EXPECT_TRUE(shift_values[i].isApprox(ws_shift.cart()));
    }
    
    EXPECT_TRUE(didnt_match>=as*bs/4);
}

TEST_F(ShiftingTest, WignerSeitzStructureMatches)
{
    EXPECT_EQ(shifted_structures.size(),wigner_seitz_shifted_structures.size());
    cu::mapping::MappingInput map_strategy;
    map_strategy.k_best_maps=0;
    map_strategy.use_crystal_symmetry=true;
    //This is a casm bug, you should be able to pass a much smaller min_cost
    map_strategy.min_cost=1e-5;

    for(int i=0; i<shifted_structures.size(); ++i)
    {
        cu::mapping::StructureMapper_f map_to_shifted(shifted_structures[i],map_strategy);
        auto results=map_to_shifted(wigner_seitz_shifted_structures[i]);

        EXPECT_TRUE(results.size()>0);
        EXPECT_TRUE(results[0].cost<1e-10);
    }
}

TEST_F(ShiftingTest, SitesHaventMoved)
{
    for(const auto& structure_set : {shifted_structures,wigner_seitz_shifted_structures})
    {
        for(const Structure& struc : structure_set)
        {
            compare_basis_to_b2(struc);
        }
    }
}

TEST_F(ShiftingTest, CategorizeShifts)
{
    auto index_map=mush::categorize_equivalently_shifted_structures(shifted_structures);

    //The origin only maps to itself
    EXPECT_TRUE(index_map[0].size()==1);
    EXPECT_TRUE(index_map[0][0]==0);

    //Everyone mapped to at least themselves
    for(const auto& indexes : index_map)
    {
        EXPECT_TRUE(indexes.size()>0);
    }

    //TODO: Is this arrangement useful elsewhere?
    std::vector<std::vector<std::size_t>> arranged_structure_indexes(as,std::vector<std::size_t>(bs));
    for(const auto& record : shift_records)
    {
        int a_ix=record.a;
        int b_ix=record.b;
        int ix=record.index;

        arranged_structure_indexes[a_ix][b_ix]=ix;
    }

    //Expect 2 fold symmetry for the b2 slice along both axis
    //start counting at 1 because the equivalent is the periodic image, which
    //doesn't get generated
    for(int a=1; a<arranged_structure_indexes.size()/2; ++a)
    {
        for(int b=1; b<arranged_structure_indexes[a].size()/2; ++b)
        {
            int pa_pb_ix=arranged_structure_indexes[a][b];          //bottom left
            int ma_pb_ix=arranged_structure_indexes[as-a][b];       //bottom righ
            int pa_mb_ix=arranged_structure_indexes[a][bs-b];       //top left
            int ma_mb_ix=arranged_structure_indexes[as-a][bs-b];    //top right

            EXPECT_EQ(index_map[pa_pb_ix],index_map[ma_pb_ix]);
            EXPECT_EQ(index_map[ma_pb_ix],index_map[pa_mb_ix]);
            EXPECT_EQ(index_map[pa_mb_ix],index_map[ma_mb_ix]);
            EXPECT_EQ(index_map[ma_mb_ix],index_map[pa_pb_ix]); //redundant
        }
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
