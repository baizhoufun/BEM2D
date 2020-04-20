#include <cmath>
#include "geometry2D.hpp"

namespace numericTools
{

Eigen::MatrixX2d Geometry2D::circle(double x0, double y0, int n, bool up)
{
    Eigen::MatrixX2d xy(n, 2);
    xy.setZero();
    double radius = sqrt(x0 * x0 + y0 * y0);
    double angle0 = acos(y0 / radius);
    double angle1 = M_PI;

    xy(0, 0) = x0;
    xy(0, 1) = y0;
    if (up)
        angle1 = 0;

    for (int i = 1; i < n; i++)
    {
        double t = ((double)i) / (n - 1.);
        double theta = angle0 + (angle1 - angle0) * t;
        xy(i, 0) = radius * sin(theta);
        xy(i, 1) = radius * cos(theta);
    }
    // snap end point to axis
    xy(n - 1, 0) = 0;
    if (up)
        xy(n - 1, 1) = radius;
    else
        xy(n - 1, 1) = -radius;
    return xy;
};

Eigen::MatrixX2d Geometry2D::line(double x0, double y0, double x1, double y1, int n)
{
    Eigen::MatrixX2d xy(n, 2);
    xy.setZero();
    for (int i = 0; i < n; i++)
    {
        double t = ((double)i) / (n - 1.);
        t = pow(t, 1.0);
        double xt = x0 + (x1 - x0) * t;
        double yt = y0 + (y1 - y0) * t;
        xy(i, 0) = xt;
        xy(i, 1) = yt;
    }
    return xy;
};

//return n-knots
Eigen::MatrixX2d Geometry2D::benchmarkShape(double angle0, double angle1, int n)
{
    Eigen::MatrixX2d xy(n, 2);
    xy.setZero();
    for (int i = 0; i < n; i++)
    {

        double t = (static_cast<double>(i)) / (n - 1.);
        double theta = angle0 + (angle1 - angle0) * t;
        xy(i, 0) = 2. * sin(theta) * (1 + 0.25 * cos(12. * theta - M_PI));
        xy(i, 1) = 2. * cos(theta) * (1 + 0.25 * cos(12. * theta - M_PI));
    }
    if (abs(angle0) < 1e-12)
        xy(0, 0) = 0.0;
    if (abs(angle1 - M_PI) < 1e-12)
        xy(n - 1, 0) = 0.0;

    return xy;
}

void Geometry2D::nrnz(double dr, double dz, double &nr, double &nz)
{
    double J = sqrt(dz * dz + dr * dr);
    nr = -dz / J;
    nz = dr / J;
};

double Geometry2D::R(double r, double z)
{
    return std::sqrt(r * r + z * z);
};

double Geometry2D::cosTheta(double r, double z)
{
    return z / Geometry2D::R(r, z);
};

void Geometry2D::zonalHarmonics(double r, double z, int L, LegendrePolyType lType, HarmonicType hType, double &phi, double &dPhiDr, double &dPhiDz)
{
    double R = Geometry2D::R(r, z);
    double cosTheta = Geometry2D::cosTheta(r, z);

    double Pl = LegendrePoly::legendreP(L, cosTheta, lType);
    double Pl1 = LegendrePoly::legendreP(L + 1, cosTheta, lType);

    double l = (double)L;
    if (lType == LegendrePolyType::HALF_INTEGER)
        l = l + 0.5;

    if (hType == HarmonicType::Inner)
    {
        phi = pow(R, l) * Pl;
        if (std::abs(r) < 1e-14)
            dPhiDr = 0;
        else
            dPhiDr = ((l * r * r - (1 + l) * z * z) * Pl + R * (1 + l) * z * Pl1) * pow(R, -2 + l) / r;
        dPhiDz = ((z + 2.0 * l * z) * Pl - R * (1 + l) * Pl1) * pow(R, -2 + l);
    }
    else
    {
        phi = pow(R, -l - 1.0) * Pl;
    }
};

} // namespace numericTools
