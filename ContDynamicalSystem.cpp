#include "ContDynamicalSystem.h"
#include "RK45.h"
#include <cmath>
#include <utility>

ContDynamicalSystem::ContDynamicalSystem(Field field): f(std::move(field)) {
}

State ContDynamicalSystem::rhs(const State& x, double t) const
{
    return f(x,t);
}

 traj ContDynamicalSystem::trajectory(const rk45& sol, const State& x0,
    double t0, double t1) const
{
    return sol.integrate(*this, x0, t0, t1);
}

 std::vector<double> ContDynamicalSystem::Clyap(const rk45& sol, const State& x0, const CJacobian& D, double tau, int segments,
     int transientseg ) const
 {
     const size_t n = x0.size();

     ContDynamicalSystem aug([this, &D, n](const State& s, double t) {
         State x(s.begin(), s.begin() + n);
         State w(s.begin() + n, s.end());

         State ds = rhs(x, t);
         matrix A = D(x, t);
         for (size_t i = 0; i < n; ++i) {
             double sum = 0.0;
             for (size_t j = 0; j < n; ++j)
                 sum += A[i][j] * w[j];
             ds.push_back(sum);
         }
         return ds;
         });



     State x = x0;
     State w(n, 1.0 / std::sqrt(double(n)));
     double t = 0.0;
     double sum = 0.0;
     std::vector<double> X1;
     X1.reserve(segments);

     for (int k = 1; k <= transientseg + segments; ++k)
     {
         
         State s = x;
         s.insert(s.end(), w.begin(), w.end());
         s = sol.integrate(aug, s, t, t + tau).x.back();
         t += tau;

         x.assign(s.begin(), s.begin() + n);
         w.assign(s.begin() + n, s.end());

         double alpha = 0.0;
         for (double v : w) alpha += v * v;
         alpha = std::sqrt(alpha);

         if (k > transientseg)
         {
             sum += std::log(alpha);
             int m = k - transientseg;
             X1.push_back(sum / (m * tau));      
         }

         for (double& v : w) v /= alpha;        
     }
     return X1;
 };


