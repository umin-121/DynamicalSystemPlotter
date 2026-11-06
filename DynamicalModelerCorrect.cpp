#define WITHOUT_NUMPY
#define _USE_MATH_DEFINES
#include <iostream>
#include "DiscreteDynamicalSystem.h"
#include <numbers>
#include <cmath>
#include <algorithm>
#include "matplotlibcpp.h"
#include "ContDynamicalSystem.h"
#include "rk45.h"
namespace plt = matplotlibcpp;

double lyapapprox(std::vector<double>& X1)
{
    return X1.back();
};


std::vector<double> plotcomp(const std::vector<State>& orbit, int i) {
    std::vector<double> out;
    out.reserve(orbit.size());
    for (const State& s : orbit)
        out.push_back(s[i]);
    return out;
}
int main()
{
    
    double r = std::numbers::pi;
    DiscreteDynamicalSystem exp([r](const State& x) {
        return State{r*std::pow(std::numbers::e,x[0])};
        }
    );
   
    double muq = 1;
    DiscreteDynamicalSystem quad([muq](const State & x) {
        return State{ muq * x[0] * (1 - x[0]) };
    });

    DJacobian quadJ = [muq](const State& x) {
        return matrix{ { muq * (1 - 2 * x[0]) }
        };
        };

    std::cout << "Here new stuff";
    const double sigma = 10.0, rho = 28.0, beta = 8.0 / 3.0;

    ContDynamicalSystem lorenz([=](const State& s, double) {
        return State{ sigma * (s[1] -s[0]),
                      s[0] * (rho - s[2]) - s[1],
                      s[0] * s[1] - beta * s[2] };
        });

    CJacobian J = [=](const State& s, double) {
        return matrix{ { -sigma,       sigma,  0.0 },
                       { rho - s[2],   -1.0,   -s[0] },
                       { s[1],         s[0],   -beta } };
        };

    rk45 solver(rk45opts{ 1e-9 });          
    
    auto X1 = lorenz.Clyap(solver, { 1.0, 1.0, 1.0 }, J, 0.5, 1000, 100);

    std::cout << X1.back() << '\n';
    std::vector<double> k(X1.size());
    for (size_t i = 0; i < k.size(); ++i) k[i] = 0.5 * (i + 1);   

    plt::plot(k, X1, "b-");
    plt::plot(std::vector<double>{k.front(), k.back()},
        std::vector<double>{0.9056, 0.9056}, "r--");        
    plt::xlabel("time");
    plt::ylabel("X1");
    plt::savefig("C:\\Users\\gabea\\Desktop\\lorenz_lyap.png");

    traj tra = lorenz.trajectory(solver, { 1.0,1.0,1.0 }, 0.0, 50.0);
    std::vector<double> X = plotcomp(tra.x, 0);
    std::vector<double> Y = plotcomp(tra.x, 1);
    std::vector<double> Z = plotcomp(tra.x, 2);
    
    

    plt::plot(X, Y, "b-");  plt::xlabel("x"); plt::ylabel("y");
    plt::savefig("C:\\Users\\gabea\\Desktop\\lorenz_xy.png");  plt::clf();

    plt::plot(X, Z, "b-");  plt::xlabel("x"); plt::ylabel("z");
    plt::savefig("C:\\Users\\gabea\\Desktop\\lorenz_xz.png");  plt::clf();

    plt::plot(Y, Z, "b-");  plt::xlabel("y"); plt::ylabel("z");
    plt::savefig("C:\\Users\\gabea\\Desktop\\lorenz_yz.png"); plt::clf();
   
    matrix A = { { 0.0, 1.0 },
             { -1.0, 0.0 } };

    ContDynamicalSystem osc([A](const State& v, double) {
        return State{ A[0][0] * v[0] + A[0][1] * v[1],
                      A[1][0] * v[0] + A[1][1] * v[1] };
        });

    CJacobian C = [A](const State&, double) { 
        return A; 
        };  

    rk45 Csolver;
    auto trc = osc.trajectory(Csolver, { 1.0, 0.0 }, 0.0, 20.0);
    auto XC = plotcomp(trc.x, 0);
    auto YC = plotcomp(trc.x, 1);
    plt::plot(XC, YC, "b-");
    plt::axis("equal");         
    plt::xlabel("x"); plt::ylabel("y");
    plt::savefig("C:\\Users\\gabea\\Desktop\\oscillator.png"); plt::clf();


    const double a = 1.0, b = 1.0, c = 1.0, d = 1.0;

    ContDynamicalSystem lv([=](const State& s, double) {
        double x = s[0];     
        double y = s[1];     
        return State{ a * x - b * x * y,
                      d * x * y - c * y };
        });

    rk45 LVsol;
    auto trLV = lv.trajectory(LVsol, { 2.0, 1.0 }, 0.0, 30.0);
    auto XLV = plotcomp(trLV.x, 0);
    auto YLV = plotcomp(trLV.x, 1);
    plt::plot(XLV, YLV, "b-");
    plt::xlabel("x"); plt::ylabel("y");
    plt::save("C:\\Users\\gabea\\Desktop\\lotkavolterra.png");
    

    std::cout << "Here quad ";
    std::vector<double> X2 = quad.lyap(quad, { 0.3 }, quadJ, 10000000, 1000, 1e-8);
    std::cout << lyapapprox(X2)- std::numeric_limits<double>::epsilon() << "\n";









    /*

    std::vector<double> nt(X2.size()), ref(X2.size(), std::log(muq));
    for (size_t i = 0; i < nt.size(); i++) 
        nt[i] = static_cast<double>(i + 1);
    plt::plot(nt, X2);
    plt::xlabel("n");
    plt::ylabel("X1(n)");
    plt::legend();
    plt::savefig("C:\\Users\\gabea\\Desktop\\X1_tent.png");

    double mut = 1.9;
    DiscreteDynamicalSystem tent([mut](const State& x) {
        return State{mut*std::min(x[0], 1.0 - x[0]) };
        }
    );
    
    Jacobian tentJ = [mut](const State& x) {
        double slope;
        if (x[0] < 0.5)
            slope = mut;
        else
            slope = -mut; 
        return matrix{ { slope } };

           
     };

    std::cout << "Here";
    std::vector<double> X1 = tent.lyap(tent, { 0.3 }, tentJ, 10000, 1000, 1e-8);
    std::cout << lyapapprox(X1)<< "\n";




    std::vector<State> orbit = quad.orbit({ .12 }, 50);
    for (int i = 0; i < 10; ++i)
        std::cout << orbit[i][0] << "\n";
    std::vector<double> n(orbit.size());
    for (size_t i = 0; i < n.size(); ++i)
        n[i] = static_cast<double>(i);
    plt::plot(n, plotcomp(orbit, 0), "o-");
    plt::xlabel("n");
    plt::ylabel("Tent Map");
    plt::savefig("C:\\Users\\gabea\\Desktop\\my_plot1.png"); // save the figure
    */
}