#ifndef DYNAMICCONE_HPP
#define DYNAMICCONE_HPP
#include <vector>
#include <cmath>
#include <eigen3/Eigen/Core>
#include "bem2D/boundaryElement.hpp"

namespace dynamicCone
{

enum class BoundaryType
{
    Interface,
    Liquid,
    Vacuum
};

class PatchedCone
{
public:
    struct Parameter
    {
        double a[5], b[5], c[5];
        double rCutoff = 50;
        int splineNodeNumber[3];
        int n[3];
        double rc = 10;
    };
    bem2D::BoundaryElement _bem[3];

private:
    Parameter _parameter;
    //const Parameter &parameter() const; // 1st coordinate

public:
    PatchedCone(double c1, double b0, double rCut, int n1, int n2, int n3);
    const Parameter &parameter() const; // 1st coordinate

    Eigen::VectorXd setVacuumBC(const bem2D::BoundaryElement &interface, const bem2D::BoundaryElement &vacuum) const;
    Eigen::VectorXd setLiquidBC(const bem2D::BoundaryElement &interface, const bem2D::BoundaryElement &liquid) const;
    const Eigen::MatrixXd farFieldShapePadding(double rEnd, double rIncrement) const; // STILL WORKING ON
    void printABC(const std::string &name) const;
    void setBoundaryElement(BoundaryType type, const Eigen::MatrixX2d &xy, int shift);
    void solveBVP(
        const bem2D::BoundaryElement &bem0, const Eigen::VectorXd &given, BoundaryType type,
        const Eigen::MatrixXd &S, const Eigen::MatrixXd &D, Eigen::VectorXd &p, Eigen::VectorXd &q);

private:
    static void computeABC(const double c1, const double b0, double *a, double *b, double *c);
    static double farFieldVelocityPotential(double r, double z, const double (&a)[5]);
    static double farFieldElectricFlux(double nr, double nz, double r, double z, const double (&b)[5]);
    static double farFieldShape(double r, const double (&c)[5]);
    static double farFieldShapeDr(double r, const double (&c)[5]);
    static double farFieldShapeArc(double r0, double r1, const double (&c)[5]);
};
} // namespace dynamicCone

#endif