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

} // namespace numericTools
