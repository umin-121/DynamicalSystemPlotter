#include "DiscreteDynamicalSystem.h"
#include <utility>



DiscreteDynamicalSystem::DiscreteDynamicalSystem(Map function)
    : f(std::move(function))
{
}

State DiscreteDynamicalSystem::evaluate(const State& x) const
{
    return f(x);
}

std::vector<State> DiscreteDynamicalSystem::orbit(const State& x0, int steps) const
{
    std::vector<State> result;
    result.reserve(steps + 1);
    result.push_back(x0);
    for (int i = 0; i < steps; ++i)
        result.push_back(f(result.back()));
    return result;
}
//https://arxiv.org/pdf/0811.0882
std::vector<double> DiscreteDynamicalSystem::lyap(const DiscreteDynamicalSystem& system, const State& x0, const DJacobian& D, int steps, int transient=1000, double d0=1e-8) {
    const size_t n = x0.size(); 
    State x = x0; 
    for (int i = 0; i < transient; ++i)
        x = system.evaluate(x);
    State w(n, 1.0 / std::sqrt(double(n)));
    std::vector<double> X1; 
    X1.reserve(steps);
    double sum = 0.0;
    for (int k = 1; k <= steps; ++k)
    {
        matrix A = D(x);
        State wn(n, 0.0);
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j) {
                wn[i] = A[i][j] * w[j];
            };
        x = system.evaluate(x);
        double alpha = 0.0;
        for (double v : wn) {
            alpha += v * v;
        }
        alpha = std::sqrt(alpha);
        sum += std::log(alpha);
        X1.push_back(sum / k);
        for (size_t i = 0; i < n; ++i) {
            w[i] = wn[i] / alpha;
        }
        return X1;
    };

}
