#include "lagrangeBasis.hpp"

namespace numericTools
{
double (**LagrangeBasis::N[3])(double) = {LagrangeBasis::N0, LagrangeBasis::N1, LagrangeBasis::N2};
double (*LagrangeBasis::N0[1])(double) = {LagrangeBasis::N00};
double (*LagrangeBasis::N1[2])(double) = {LagrangeBasis::N10, LagrangeBasis::N11};
double (*LagrangeBasis::N2[3])(double) = {LagrangeBasis::N20, LagrangeBasis::N21, LagrangeBasis::N22};

double LagrangeBasis::N00(double x) { return 0; }
// Lagrange interpolation basis first order eqn (7.59)
double LagrangeBasis::N10(double x) { return 1. - x; }
double LagrangeBasis::N11(double x) { return x; }
// Lagrange interpolation basis second order eqn (7.60)
double LagrangeBasis::N20(double x) { return 1. - 3. * x + 2. * x * x; }
double LagrangeBasis::N21(double x) { return 4. * x - 4. * x * x; }
double LagrangeBasis::N22(double x) { return 2. * x * x - x; }
} // namespace numericTools