#ifndef MULTISHIFTFOURIER_HH
#define MULTISHIFTFOURIER_HH

#include "./define.hpp"
#include "./surf.hpp"
#include "casm/crystallography/Lattice.hh"
#include <complex>
#include <vector>

namespace mush
{
    /**
     * Given a lattice (only ab-vectors matter) and a list of values
     * and grid points, creates a reciprocal space from which to sample
     * plane waves uniformly, and calculates coefficients to
     * reproduce the given data.
     *
     * The reciprocal grid to sample in is determined by the shape of
     * the 2d arrangement of the interpolated values.
     */

    class Interpolator
    {
        public:

            typedef std::vector<std::vector<InterPoint>> InterGrid;

            Interpolator(const Lattice& init_lat, const InterGrid& init_values);

            ///Number of k-points/data points
            int size() const;

            /// Use the Fourier basis to reconstruct the signal at an arbitrary resolution
            std::pair<Lattice, InterGrid> interpolate(int a_dim, int b_dim) const;

        private:

            /// The real lattice where the gamma surface happened on the ab-plane
            /// This is not the true lattice of the structure, just one with the
            /// correct ab-vectors oriented along the xy-plane
            Lattice m_real_lat;

            /// The reciprocal lattice of the reduced real lattice
            Lattice m_recip_lat;

            /// Holds initial values at the real interpolation points
            InterGrid m_real_ipoints;

            /// Holds coefficients at each of the reciprocal k points
            InterGrid m_k_values;

            /// Constructs grid of k points, setting all coefficients to zero
            static InterGrid _k_grid(const InterGrid& init_values, const Lattice& reciprocal_lattice);

            /// Sets the coefficients for each of the k points
            void _take_fourier_transform();

            /// Helper to avoid so many loops. Unrolls 2d vector into a 1d vector.
            static std::vector<InterPoint*> _unrolled_values(InterGrid* grid_values);

            /// Helper to avoid so many loops. Unrolls 2d vector into a 1d vector, put pointers are const.
            static std::vector<const InterPoint*> _unrolled_values(const InterGrid& grid_values);
    };
}

#endif
