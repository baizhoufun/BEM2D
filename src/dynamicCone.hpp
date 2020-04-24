#ifndef DYNAMICCONE_HPP
#define DYNAMICCONE_HPP
#include <vector>
#include <cmath>
#include <eigen3/Eigen/Core>
#include "bem2D/boundaryElement.hpp"

class DynamicCone
{
public:
    struct Parameter
    {
        double a[5], b[5], c[5];
        double rCutoff = 50;
        int n[2];
    } parameter;
    bem2D::BoundaryElement bem[3];

public:
    DynamicCone();
    Eigen::VectorXd setVacuumBC(const bem2D::BoundaryElement &interface, const bem2D::BoundaryElement &vacuum) const;
    Eigen::VectorXd setLiquidBC(const bem2D::BoundaryElement &interface, const bem2D::BoundaryElement &liquid) const;

private:
    static void computeABC(const double c1, const double b0, double *a, double *b, double *c);
    static double farFieldVelocityPotential(double r, double z, const double (&a)[5]);
    static double farFieldElectricFlux(double nr, double nz, double r, double z, const double (&b)[5]);
};

#endif