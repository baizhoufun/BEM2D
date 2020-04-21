#ifndef GEOMETRY2D_HPP
#define GEOMETRY2D_HPP

#include <eigen3/Eigen/Dense>
#include "legendrePoly.hpp"

namespace numericTools
{
enum class HarmonicType
{
    Inner,
    Outer
};

enum class CircleType
{
    UpWrap,
    DownWrap
};

class Geometry2D
{

public:
    static Eigen::MatrixX2d c3Cone(double r_c, double r_star, const double (&c)[5], int n, double (*f)(double t) = nullptr);
    static Eigen::MatrixX2d circle(double x0, double y0, int n, CircleType type, double (*f)(double t) = nullptr);
    static Eigen::MatrixX2d line(double x0, double y0, double x1, double y1, int n, double (*f)(double t) = nullptr);
    static Eigen::MatrixX2d benchmarkShape(double angle0, double angle1, int n, double (*f)(double t) = nullptr);
    static bool isOnZAxis(double r);
    static void nrnz(double dr, double dz, double &nr, double &nz);
    static double R(double r, double z);
    static double cosTheta(double r, double z);
    static void zonalHarmonics(double r, double z, int l, LegendrePolyType lType, HarmonicType hType, double &phi, double &dPhiDr, double &dPhiDz);
    static double c3Cone(double r, double rc, const double (&c)[5]);
    static double finiteDifference(int order, double h[3], double y[4], int location);

    static void estimateDerivativeCircleBegin(const Eigen::MatrixX2d &xy, double &dx, double &ddx, double &dy, double &ddy);
    static void estimateDerivativeConeEnd(const Eigen::MatrixX2d &xy, const double (&c)[5], double &dx, double &ddx, double &dy, double &ddy);
};

} // namespace numericTools
#endif