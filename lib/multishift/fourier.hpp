#ifndef MULTISHIFTFOURIER_HH
#define MULTISHIFTFOURIER_HH

#include <complex>
#include <vector>

namespace mush
{
    /**
     * One point on the surface that should be used for the interpolation.
     * Only fractional coordinates needed.
     */

    class InterPoint
    {
        public:

            InterPoint(double init_a_frac, double init_b_frac, std::complex<double> init_value):
                a_frac(init_a_frac), b_frac(init_b_frac), value(init_value){}

            double a_frac;
            double b_frac;
            std::complex<double> value;

        private:
    };

    /**
     * Given a lattice (only ab-vectors matter) and a list of values
     * and grid points, creates a reciprocal space from which to sample
     * plane waves uniformly, and calculates coefficients to
     * reproduce the given data.
     */

    class Interpolator
    {
        public:

        private:

            /// Holds initial values at the real interpolation points
            std::vector<InterPoint> m_real_ipoints;

            /// Holds coefficients at each of the reciprocal k points
            std::vector<InterPoint> m_k_values;
    };
}

#endif
