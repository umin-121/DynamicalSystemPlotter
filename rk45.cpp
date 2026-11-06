#include "rk45.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <string>
rk45::rk45(rk45opts options) : opts(options) {}

// https://math.okstate.edu/people/yqwang/teaching/math4513_fall11/Notes/rungekutta.pdf
 traj rk45::integrate(const ContDynamicalSystem& sys, const State& x0, double t0, double t1) const
{
    const size_t n = x0.size();

    traj result;                    
    result.t.push_back(t0);
    result.x.push_back(x0);

    double t = t0;
    State w = x0;
    double h = opts.h0;
    State tmp(n);

    
    auto stage = [&](double tt, const State& y) {
        State k = sys.rhs(y, tt);   
        for (double& v : k) v *= h;
        return k;
        };

    int steps = 0;
    while (t < t1)
    {
        if (++steps > opts.maxSteps)
            throw  "rk45: maxSteps exceeded";

        bool last = false;
        if (t + h >= t1) {
            h = t1 - t;
            last = true;
        }

        State k1 = stage(t, w);

        for (size_t i = 0; i < n; ++i)
            tmp[i] = w[i] + k1[i] / 4;
        State k2 = stage(t + h / 4, tmp);

        for (size_t i = 0; i < n; ++i)
            tmp[i] = w[i] + 3 * k1[i] / 32 + 9 * k2[i] / 32;
        State k3 = stage(t + 3 * h / 8, tmp);

        for (size_t i = 0; i < n; ++i)
            tmp[i] = w[i] + 1932 * k1[i] / 2197 - 7200 * k2[i] / 2197 + 7296 * k3[i] / 2197;
        State k4 = stage(t + 12 * h / 13, tmp);

        for (size_t i = 0; i < n; ++i)
            tmp[i] = w[i] + 439 * k1[i] / 216 - 8 * k2[i] + 3680 * k3[i] / 513 - 845 * k4[i] / 4104;
        State k5 = stage(t + h, tmp);

        for (size_t i = 0; i < n; ++i)
            tmp[i] = w[i] - 8 * k1[i] / 27 + 2 * k2[i] - 3544 * k3[i] / 2565
            + 1859 * k4[i] / 4104 - 11 * k5[i] / 40;
        State k6 = stage(t + h / 2, tmp);

        State w1(n), w2(n);
        for (size_t i = 0; i < n; ++i)
        {
            w1[i] = w[i] + 25 * k1[i] / 216 + 1408 * k3[i] / 2565
                + 2197 * k4[i] / 4104 - k5[i] / 5;
            w2[i] = w[i] + 16 * k1[i] / 135 + 6656 * k3[i] / 12825 + 28561 * k4[i] / 56430
                - 9 * k5[i] / 50 + 2 * k6[i] / 55;
        }

        double R = 0.0;
        for (size_t i = 0; i < n; ++i)
        {
            double err = k1[i] / 360 - 128 * k3[i] / 4275 - 2197 * k4[i] / 75240
                + k5[i] / 50 + 2 * k6[i] / 55;
            R = std::max(R, std::abs(err) / h);
        }

        double delta;
        if (R == 0.0)
            delta = 4.0;
        else
            delta = std::clamp(0.84 * std::pow(opts.epsilon / R, 0.25), 0.1, 4.0);

        if (R <= opts.epsilon)
        {
            
            if (last) {
                t = t1;
            }
            else
            {
                t = t + h;
            }
            w = std::move(w1);
            result.t.push_back(t);
            result.x.push_back(w);
            if (last)
                break;
        }

        h = std::min(delta * h, opts.hMax);
        if (h < opts.hMin)
            throw std::string("rk45: step underflow at t=" + std::to_string(t)
                + "  t1-t=" + std::to_string(t1 - t)
                + "  R=" + std::to_string(R));
    }
    return result;
}