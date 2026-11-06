#pragma once
#include <functional>
#include <vector>

using State = std::vector<double>;
using matrix = std::vector<std::vector<double>>;
using CJacobian = std::function<matrix(const State&, double)>;  
using Field = std::function<State(const State&, double)>;   

struct traj {
    std::vector<double> t;
    std::vector<State>  x;
};

class rk45;   
class ContDynamicalSystem
{
public:
    explicit ContDynamicalSystem(Field function);

    State rhs(const State& x, double t) const;

    traj trajectory(const rk45& sol, const State& x0, double t0, double t1) const;

    std::vector<double> Clyap(const rk45& sol, const State& x0, const CJacobian& D,
        double tau, int segments, int transientseg = 1000) const;

private:
    Field f;
};