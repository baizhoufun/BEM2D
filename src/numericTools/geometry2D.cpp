#include <cmath>
#include "geometry2D.hpp"

namespace numericTools
{

Eigen::MatrixX2d Geometry2D::circle(double x0, double y0, int n, CircleType type, double (*f)(double t))
{
    Eigen::MatrixX2d xy(n, 2);
    xy.setZero();
    double radius = sqrt(x0 * x0 + y0 * y0);
    double angle0 = acos(y0 / radius);
    double angle1 = M_PI;

    xy(0, 0) = x0;
    xy(0, 1) = y0;

    if (type == CircleType::UpWrap)
        angle1 = 0;
    else if (type == CircleType::DownWrap)
        angle1 = M_PI;

    for (int i = 1; i < n; i++)
    {
        double t = ((double)i) / (n - 1.);
        if (f != nullptr)
            t = (*f)(t);
        double theta = angle0 + (angle1 - angle0) * t;
        xy(i, 0) = radius * sin(theta);
        xy(i, 1) = radius * cos(theta);
    }
    // snap end point to axis
    xy(n - 1, 0) = 0;
    if (type == CircleType::UpWrap)
        xy(n - 1, 1) = radius;
    else if (type == CircleType::DownWrap)
        xy(n - 1, 1) = -radius;
    return xy;
};

Eigen::MatrixX2d Geometry2D::line(double x0, double y0, double x1, double y1, int n, double (*f)(double t))
{
    Eigen::MatrixX2d xy(n, 2);
    xy.setZero();
    for (int i = 0; i < n; i++)
    {
        double t = ((double)i) / (n - 1.);
        //        t = pow(t, 1.0);
        if (f != nullptr)
            t = (*f)(t);
        double xt = x0 + (x1 - x0) * t;
        double yt = y0 + (y1 - y0) * t;
        xy(i, 0) = xt;
        xy(i, 1) = yt;
    }
    return xy;
};

//return n-knots
Eigen::MatrixX2d Geometry2D::benchmarkShape(double angle0, double angle1, int n, double (*f)(double t))
{
    Eigen::MatrixX2d xy(n, 2);
    xy.setZero();
    for (int i = 0; i < n; i++)
    {
        double t = (static_cast<double>(i)) / (n - 1.);
        if (f != nullptr)
            t = (*f)(t);
        double theta = angle0 + (angle1 - angle0) * t;
        xy(i, 0) = 2. * sin(theta) * (1 + 0.25 * cos(8. * theta - M_PI));
        xy(i, 1) = 2. * cos(theta) * (1 + 0.25 * cos(8. * theta - M_PI));
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

    if (Geometry2D::R(r, z) < 1e-13)
        return 0;
    else
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
        if (isOnZAxis(r))
            dPhiDr = 0;
        else
            dPhiDr = ((l * r * r - (1 + l) * z * z) * Pl + R * (1 + l) * z * Pl1) * pow(R, -2 + l) / r;
        dPhiDz = ((z + 2.0 * l * z) * Pl - R * (1 + l) * Pl1) * pow(R, -2 + l);
    }
    else
    {
        phi = pow(R, -l - 1.0) * Pl;
        if (isOnZAxis(r))
            dPhiDr = 0;
        else
            dPhiDr = -pow(R, -3.0 - l) * (1 + l) * (R * R * Pl - R * z * Pl1) / r;
        dPhiDz = -pow(R, -2.0 - l) * (1 + l) * Pl1;
    }
};

bool Geometry2D::isOnZAxis(double r)
{
    if (std::abs(r) < 1e-13)
        return true;
    else
        return false;
};

Eigen::MatrixX2d Geometry2D::c3Cone(double r_c, double r_star, const double (&c)[5], int n, double (*f)(double t))
{
    Eigen::MatrixX2d xy(n, 2);
    xy.setZero();
    for (int i = 0; i < n; i++)
    {
        double t = ((double)i) / (n - 1.);
        if (f != nullptr)
            t = (*f)(t);

        double r = r_star * t;
        xy(i, 0) = r;
        xy(i, 1) = c3Cone(r, r_c, c);
    }
    return xy;
}

double Geometry2D::c3Cone(double r, double rc, const double (&c)[5])
{
    //double f0, f2, f4, f6;
    double f0 = (231 * c[4]) / (16. * pow(rc, 5)) + (1045 * c[3]) / (128. * pow(rc, 3.5)) + (195 * c[1]) / (128. * sqrt(rc)) + (5 * c[0] * rc) / 16.;
    double f2 = (-495 * c[4]) / (16. * pow(rc, 7)) - (1995 * c[3]) / (128. * pow(rc, 5.5)) - (117 * c[1]) / (128. * pow(rc, 2.5)) + (15 * c[0]) / (16. * rc);
    double f4 = (385 * c[4]) / (16. * pow(rc, 9)) + (1463 * c[3]) / (128. * pow(rc, 7.5)) + (65 * c[1]) / (128. * pow(rc, 4.5)) - (5 * c[0]) / (16. * pow(rc, 3));
    double f6 = (-105 * c[4]) / (16. * pow(rc, 11)) - (385 * c[3]) / (128. * pow(rc, 9.5)) - (15 * c[1]) / (128. * pow(rc, 6.5)) + c[0] / (16. * pow(rc, 5));
    if (r < rc)
    {
        return f0 + f2 * pow(r, 2.0) + f4 * pow(r, 4.0) + f6 * pow(r, 6.0);
    }
    else
    {
        return c[0] * r + c[1] / sqrt(r) + c[3] / pow(r, 3.5) + c[4] / pow(r, 5.0);
    }
}

double Geometry2D::curvatureC3Cone(double r, double rc, const double (&c)[5])
{
    //double f0 = (231 * c[4]) / (16. * pow(rc, 5)) + (1045 * c[3]) / (128. * pow(rc, 3.5)) + (195 * c[1]) / (128. * sqrt(rc)) + (5 * c[0] * rc) / 16.;
    double f2 = (-495 * c[4]) / (16. * pow(rc, 7)) - (1995 * c[3]) / (128. * pow(rc, 5.5)) - (117 * c[1]) / (128. * pow(rc, 2.5)) + (15 * c[0]) / (16. * rc);
    double f4 = (385 * c[4]) / (16. * pow(rc, 9)) + (1463 * c[3]) / (128. * pow(rc, 7.5)) + (65 * c[1]) / (128. * pow(rc, 4.5)) - (5 * c[0]) / (16. * pow(rc, 3));
    double f6 = (-105 * c[4]) / (16. * pow(rc, 11)) - (385 * c[3]) / (128. * pow(rc, 9.5)) - (15 * c[1]) / (128. * pow(rc, 6.5)) + c[0] / (16. * pow(rc, 5));
    if (r < rc)
    {
        return (2 * (f2 + 6 * f4 * pow(r, 2) + 15 * f6 * pow(r, 4) + (f2 + 2 * f4 * pow(r, 2) + 3 * f6 * pow(r, 4)) * (1 + pow(2 * f2 * r + 4 * f4 * pow(r, 3) + 6 * f6 * pow(r, 5), 2)))) / pow(1 + pow(2 * f2 * r + 4 * f4 * pow(r, 3) + 6 * f6 * pow(r, 5), 2), 1.5);
    }
    else
    {
        return ((3 * c[1]) / (4. * pow(r, 2.5)) + (15.75 * c[3]) / pow(r, 5.5) + (30 * c[4]) / pow(r, 7) + ((c[0] - c[1] / (2. * pow(r, 1.5)) - (3.5 * c[3]) / pow(r, 4.5) - (5 * c[4]) / pow(r, 6)) * (1 + pow(c[0] - c[1] / (2. * pow(r, 1.5)) - (3.5 * c[3]) / pow(r, 4.5) - (5 * c[4]) / pow(r, 6), 2))) / r) / pow(1 + pow(c[0] - c[1] / (2. * pow(r, 1.5)) - (3.5 * c[3]) / pow(r, 4.5) - (5 * c[4]) / pow(r, 6), 2), 1.5);
    }
}

void Geometry2D::estimateDerivativeConeEnd(const Eigen::MatrixX2d &xy, const double (&c)[5], double &dx, double &ddx, double &dy, double &ddy)
{
    const int n = xy.rows();
    double x0 = xy(n - 4, 0), y0 = xy(n - 4, 1), x1 = xy(n - 3, 0), y1 = xy(n - 3, 1), x2 = xy(n - 2, 0), y2 = xy(n - 2, 1), x3 = xy(n - 1, 0), y3 = xy(n - 1, 1);
    double h0 = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
    double h1 = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
    double h2 = sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2));
    double h[3] = {h0, h1, h2};
    double x[4] = {x0, x1, x2, x3};

    double r = x3;
    //double curvature = (0.75 * pow(r, 4.5) * c[1] + 15.75 * pow(r, 1.5) * c[3] + pow(r, 6) * (1. + pow(c[0] + (-0.5 * pow(r, 4.5) * c[1] - 3.5 * pow(r, 1.5) * c[3] - 5. * c[4]) / pow(r, 6), 2)) * (c[0] + (-0.5 * pow(r, 4.5) * c[1] - 3.5 * pow(r, 1.5) * c[3] - 5. * c[4]) / pow(r, 6)) + 30. * c[4]) / (pow(r, 7) * pow(1. + pow(c[0] + (-0.5 * pow(r, 4.5) * c[1] - 3.5 * pow(r, 1.5) * c[3] - 5. * c[4]) / pow(r, 6), 2), 1.5));
    double dydx = c[0] - c[1] / (2. * pow(r, 1.5)) - (7. * c[3]) / (2. * pow(r, 4.5)) - (5. * c[4]) / pow(r, 6);
    double ddyddx = (3. * c[1]) / (4. * pow(r, 2.5)) + (63. * c[3]) / (4. * pow(r, 5.5)) + (30. * c[4]) / (pow(r, 7.0));

    // finite difference formula non uniform grid at x3
    dx = finiteDifference(1, h, x, 3);
    ddx = finiteDifference(2, h, x, 3);
    dy = dx * dydx; // eqn 7.92
    //ddy = (ddx * dy + curvature * pow(dx * dx + dy * dy, 1.5)) / dx - dy * (dx * dx + dy * dy) / dx / x3; // eqn 7.92
    ddy = dydx * ddx + dx * dx * ddyddx;
}

double Geometry2D::finiteDifference(int order, double h[3], double y[4], int location)
{
    double h0 = h[0], h1 = h[1], h2 = h[2];
    double y0 = y[0], y1 = y[1], y2 = y[2], y3 = y[3];
    switch (order)
    {
    case 1:
    {
        switch (location)
        {
        case 0:
        {
            return ((-1. * (2. * h0 + h1)) / (h0 * (h0 + h1)) - 1. / (h0 + h1 + h2)) * y0 + ((h0 + h1) / (h0 * h1) + (h0 + h1) / (h1 * (h1 + h2))) * y1 + ((-1. * h0) / (h1 * (h0 + h1)) - (1. * h0) / (h1 * h2)) * y2 + (h0 / (h1 * h2) + (-1. * h0 - 1. * h1) / (h1 * (h1 + h2)) + 1 / (h0 + h1 + h2)) * y3;
            break;
        }
        case 3:
        {
            return (1 / (h0 + h1) - (1. * h2) / (h0 * (h0 + h1)) - 1. / (h0 + h1 + h2)) * y0 + (1 / h1 + h2 / (h0 * h1) - 1. / (h1 + h2)) * y1 + ((-1. * (h0 + 2. * h1)) / (h1 * (h0 + h1)) - 1. / h2 - (1. * h2) / (h1 * (h0 + h1))) * y2 + (1 / h2 + 1 / (h1 + h2) + 1 / (h0 + h1 + h2)) * y3;
            break;
        }
        }
        break;
    }
    case 2:
    {
        switch (location)
        {
        case 0:
        {
            return (2. / (h0 * (h0 + h1)) + (2. * (2. * h0 + h1)) / (h0 * (h0 + h1) * (h0 + h1 + h2))) * y0 + (-2. / (h0 * h1) - (2. * (2. * h0 + h1)) / (h0 * h1 * (h1 + h2))) * y1 + (2. / (h1 * (h0 + h1)) + (2. * (2. * h0 + h1)) / (h1 * (h0 + h1) * h2)) * y2 + ((-2. * (2. * h0 + h1)) / (h1 * (h0 + h1) * h2) + (2. * (2. * h0 + h1)) / (h0 * h1 * (h1 + h2)) - (2. * (2. * h0 + h1)) / (h0 * (h0 + h1) * (h0 + h1 + h2))) * y3;
            break;
        }
        case 3:
        {
            return (-4. / (h0 * (h0 + h1)) + (2. * (2. * h0 + h1)) / (h0 * (h0 + h1) * (h0 + h1 + h2))) * y0 + (4. / (h0 * h1) + (2. * (h0 - 1. * h1)) / (h0 * h1 * (h1 + h2))) * y1 + (-4. / (h1 * (h0 + h1)) - (2. * (h0 + 2. * h1)) / (h1 * (h0 + h1) * h2)) * y2 + ((2. * (h0 + 2. * h1)) / (h1 * (h0 + h1) * h2) - (2. * (h0 - 1. * h1)) / (h0 * h1 * (h1 + h2)) - (2. * (2. * h0 + h1)) / (h0 * (h0 + h1) * (h0 + h1 + h2))) * y3;
            break;
        }
        }
        return 0;
        break;
    }
    default:
        break;
    }
    return 0;
};

void Geometry2D::estimateDerivativeCircleBegin(const Eigen::MatrixX2d &xy, double &dx, double &ddx, double &dy, double &ddy)
{
    double x0 = xy(0, 0), y0 = xy(0, 1), x1 = xy(1, 0), y1 = xy(1, 1), x2 = xy(2, 0), y2 = xy(2, 1), x3 = xy(3, 0), y3 = xy(3, 1);
    double h0 = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
    double h1 = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
    double h2 = sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2));
    double h[3] = {h0, h1, h2};
    double x[4] = {x0, x1, x2, x3};

    //double curvature = -2.0 / sqrt(x0 * x0 + y0 * y0);
    double curvature = -1.0 / sqrt(x0 * x0 + y0 * y0);
    double slope = -x0 / y0;

    // finite difference formula non uniform grid at x0
    dx = finiteDifference(1, h, x, 0);
    ddx = finiteDifference(2, h, x, 0);
    dy = dx * slope; // eqn 7.92

    //ddy = (ddx * dy + curvature * pow(dx * dx + dy * dy, 1.5)) / dx - dy * (dx * dx + dy * dy) / dx / x0; // eqn 7.92
    ddy = (ddx * dy + curvature * pow(dx * dx + dy * dy, 1.5)) / dx; // eqn 7.92
};

double Geometry2D::curvatureAxisymmetric(double r, double z, double dr, double dz, double ddr, double ddz)
{
    if (isOnZAxis(r))
        return ddz / dr / dr * 2.;
    else
        return (dr * ddz - dz * ddr) / pow(dr * dr + dz * dz, 1.5) + dz / r / sqrt(dr * dr + dz * dz);
}

} // namespace numericTools
