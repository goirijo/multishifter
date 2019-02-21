#include "casm/crystallography/Coordinate.hh"
#include "./fourier.hpp"
#include "./exceptions.hpp"

namespace mush
{
    Interpolator::InterGrid Interpolator::_k_grid(const InterGrid& init_values, const Lattice& reciprocal_lattice)
    {
        InterGrid k_values;

        /* std::cout<<"REAL"<<std::endl; */
        /* std::cout<<m_real_lat.lat_column_mat()<<std::endl; */
        /* std::cout<<"RECIP"<<std::endl; */
        /* std::cout<<m_recip_lat.lat_column_mat()<<std::endl; */
        auto num_as=init_values.size();
        auto num_bs=init_values[0].size();

        //This is the value you need to translate the directly constructed grid to have
        //the k-points centered around the gamma point
        //TODO: Decide if using this makes any sense at all
        /* double ka_centrize=(1.0-1.0/num_as)/2.0; */
        /* double kb_centrize=(1.0-1.0/num_bs)/2.0; */
        /* CASM::Coordinate center_shift(ka_centrize,kb_centrize,0.0,reciprocal_lattice,CASM::FRAC); */

        for(int a=0; a<init_values.size(); ++a)
        {
            InterGrid::value_type kb_row;
            double ka_frac=static_cast<double>(a)/num_as;

            const auto& a_row=init_values[a];
            if(a_row.size()!=num_bs)
            {
                throw except::DimensionalMismatch(a_row.size(),num_bs,"The grid of values to interpolate is not uniform.");
            }

            for(int b=0; b<a_row.size(); ++b)
            {
                double kb_frac=static_cast<double>(b)/num_bs;
                CASM::Coordinate k_coord(ka_frac,kb_frac,0.0,reciprocal_lattice,CASM::FRAC);
                /* k_coord-=center_shift; */

                /* std::cout<<k_coord.const_cart().transpose()<<"    "; */
                k_coord.voronoi_within();
                /* std::cout<<k_coord.const_cart().transpose(); */
                /* auto val=k_coord.voronoi_number(); */
                /* std::cout<<"    "<<val; */
                /* std::cout<<std::endl; */

                kb_row.emplace_back(k_coord.const_frac()(0),k_coord.const_frac()(1),0.0);
            }
            k_values.emplace_back(std::move(kb_row));
        }
        return k_values;
    }

    std::vector<InterPoint*> Interpolator::_unrolled_values(InterGrid* grid_values_ptr)
    {
        std::vector<InterPoint*> unrolled_values;
        auto& grid_values=*grid_values_ptr;
        for(auto& row : grid_values)
        {
            for(auto& val : row)
            {
                unrolled_values.push_back(&val);
            }
        }
        return unrolled_values;
    }

    std::vector<const InterPoint*> Interpolator::_unrolled_values(const InterGrid& grid_values)
    {
        //TODO: There's some clever way to avoid the code duplication with some template or const_cast magic, but
        //I can't be bothered right now.
        std::vector<const InterPoint*> unrolled_values;
        for(auto& row : grid_values)
        {
            for(auto& val : row)
            {
                unrolled_values.push_back(&val);
            }
        }
        return unrolled_values;
    }

    void Interpolator::_take_fourier_transform() 
    {
        std::complex<double> im(0,1);
        std::complex<double> normalization=1.0/this->size();

        auto real_vals=_unrolled_values(m_real_ipoints);
        auto k_vals=_unrolled_values(&m_k_values);

        assert(real_vals.size()==k_vals.size()==this->size());

        for(auto k_val : k_vals)
        {
            k_val->value=0.0;
            auto k_vec=k_val->cart(m_recip_lat);
            for(auto r_val : real_vals)
            {
                auto r_vec=r_val->cart(m_real_lat);
                k_val->value+=normalization*r_val->value*std::exp(-im*k_vec.dot(r_vec));
            }
        }

        return;
    }

    Interpolator::Interpolator(const Lattice& init_lat, const InterGrid& init_values):
        m_real_lat(phony_lat_mat_reorient(init_lat.lat_column_mat())),
        m_recip_lat(this->m_real_lat.get_reciprocal()),
        m_real_ipoints(init_values),
        m_k_values(this->_k_grid(init_values, this->m_recip_lat))
    {
        this->_take_fourier_transform();
    }

    int Interpolator::size() const
    {
        for(const auto& k_row : m_k_values)
        {
            assert(k_row.size()==m_k_values[0].size());
        }

        assert(m_k_values.size()==m_real_ipoints.size());

        return m_k_values.size()*m_k_values[0].size();
    }
}
