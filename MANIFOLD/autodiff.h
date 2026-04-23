
#ifndef AUTODIFF_H
#define AUTODIFF_H

#include <boost/math/differentiation/autodiff.hpp>

using namespace boost::math::differentiation;

#define MAX_DIFF_ORDER 1
using real1 = detail::fvar<double, MAX_DIFF_ORDER>;
using real2 = detail::fvar<real1, 0>;
using real3 = detail::fvar<real2, 0>;
using real = promote<real1, real2, real3>;

template <typename F, typename... P>
double evaluate(F func, const double pos[3], P... params)
{
    const auto variables =
        make_ftuple<double, MAX_DIFF_ORDER, MAX_DIFF_ORDER, MAX_DIFF_ORDER>(
            pos[0], pos[1], pos[2]);
    const auto &x = std::get<0>(variables);
    const auto &y = std::get<1>(variables);
    const auto &z = std::get<2>(variables);

    const real u = func(x, y, z, params...);
    return u.derivative(0, 0, 0);
}

template <typename F, typename... P>
void gradient(F func, double grad[3], const double pos[3], P... params)
{

    const auto variables =
        make_ftuple<double, MAX_DIFF_ORDER, MAX_DIFF_ORDER, MAX_DIFF_ORDER>(
            pos[0], pos[1], pos[2]);
    const auto &x = std::get<0>(variables);
    const auto &y = std::get<1>(variables);
    const auto &z = std::get<2>(variables);

    const real u = func(x, y, z, params...);

    grad[0] = u.derivative(1, 0, 0);
    grad[1] = u.derivative(0, 1, 0);
    grad[2] = u.derivative(0, 0, 1);
}

#endif
