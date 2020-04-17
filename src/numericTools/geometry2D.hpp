#ifndef GEOMETRY2D_HPP
#define GEOMETRY2D_HPP

#include <eigen3/Eigen/Dense>

namespace numericTools
{

class Geometry2D
{
public:
    static Eigen::MatrixX2d circle(double x0, double y0, int n, bool up);
    static Eigen::MatrixX2d line(double x0, double y0, double x1, double y1, int n);
    static Eigen::MatrixX2d benchmarkShape(double angle0, double angle1, int n);
};
} // namespace numericTools
#endif