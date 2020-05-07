#include <cmath>
#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include "numericTools/quadratureRules.hpp"

using namespace numericTools;

Eigen::IOFormat fmt(Eigen::FullPrecision, 0, "\t", "\n", "", "", "", "");
std::ofstream file;

int main()
{
    const auto f = [](double x) -> double { return cos(10 * x) * exp(-x); };
    int orderMax = 21;
    Eigen::MatrixX4d quadratureIntegral(orderMax, 4);
    quadratureIntegral.setZero();

    double tstar = 0.2; // tstar

    for (int m = 1; m < orderMax; m++)
    {
        quadratureIntegral(m, 0) = m;
        const double *ab;
        const double *wt;

        for (int k = 0; k < m; k++)
        {
            // standard Gauss-Legendre : integrate f(t) over [0,1]
            ab = QuadratureRules::abascissa(m, QuadratureType::GAUSS_LEGENDRE);
            wt = QuadratureRules::weight(m, QuadratureType::GAUSS_LEGENDRE);
            quadratureIntegral(m, 1) += wt[k] * f(ab[k]);
            // integrate f(t)log(ts-t) over [0,1], ts = 0.5
            ab = QuadratureRules::abascissa(m, QuadratureType::GAUSS_LEGENDRE);
            wt = QuadratureRules::weight(m, QuadratureType::GAUSS_LEGENDRE);
            quadratureIntegral(m, 2) += tstar * log(tstar) * wt[k] * f(tstar * ab[k]);
            quadratureIntegral(m, 2) += (1 - tstar) * log(1 - tstar) * wt[k] * f(tstar + (1 - tstar) * ab[k]); // Gauss-Legendre part
            ab = QuadratureRules::abascissa(m, QuadratureType::GAUSS_LEGENDRE_LOG);
            wt = QuadratureRules::weight(m, QuadratureType::GAUSS_LEGENDRE_LOG);
            quadratureIntegral(m, 2) -= tstar * wt[k] * f(tstar * (1 - ab[k]));
            quadratureIntegral(m, 2) -= (1 - tstar) * wt[k] * f(tstar + (1 - tstar) * ab[k]); // Log-weighted part
            ab = QuadratureRules::abascissa(m, QuadratureType::GAUSS_LEGENDRE);
            wt = QuadratureRules::weight(m, QuadratureType::GAUSS_LEGENDRE);
            // integrate f(t)log(ts-t) over [0,1], only use Gauss-Legendre
            quadratureIntegral(m, 3) += tstar * wt[k] * f(tstar * ab[k]) * log(tstar - tstar * ab[k]);
            quadratureIntegral(m, 3) += (1 - tstar) * wt[k] * f(tstar + (1 - tstar) * ab[k]) * log((1 - tstar) * ab[k]);
        }
    }

    file.open("../resources/testNumericTools/quadratureIntegral.txt");
    file << quadratureIntegral.format(fmt);
    file.close();

    return 0;
}