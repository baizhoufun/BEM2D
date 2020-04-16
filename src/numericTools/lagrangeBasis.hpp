#ifndef LAGRANGEBASIS_HPP
#define LAGRANGEBASIS_HPP

#define MAX_NUMBER_LAGRANGE_BASIS 3

namespace numericTools
{

class LagrangeBasis
{
public:
    static double (**N[MAX_NUMBER_LAGRANGE_BASIS])(double); // Lagrange interpolation basis functions
private:
    static double N00(double x);
    static double N10(double x);
    static double N11(double x);
    static double N20(double x);
    static double N21(double x);
    static double N22(double x);
    static double (*N0[1])(double);
    static double (*N1[2])(double);
    static double (*N2[3])(double);
};
} // namespace numericTools

#endif