#pragma once
#include <functional>
#include <vector>

    using State = std::vector<double>;
    using matrix = std::vector<std::vector<double>>;
    using DJacobian = std::function<matrix(const State&)>;
    using Map = std::function<State(const State&)>;
    class DiscreteDynamicalSystem
    {
    public:
     

        explicit DiscreteDynamicalSystem(Map function);

        State evaluate(const State& x) const;

        
        std::vector<State> orbit(const State& x0, int steps) const;

        std::vector<double> lyap(const DiscreteDynamicalSystem& system, const State& x0, const DJacobian& D, int steps, int transient, double d0);

        

    private:
        Map f;
    };


