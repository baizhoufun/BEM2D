#include <cmath>
#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include "boundaryElement.hpp"
#include "numericTools/geometry2D.hpp"

Eigen::IOFormat fmt(Eigen::FullPrecision, 0, "\t", "\n", "", "", "", "");
std::ofstream file;

int main()
{
    int N = 128;

    Eigen::MatrixXd knots = numericTools::Geometry2D::benchmarkShape(0, M_PI, N + 1);

    bem2D::BoundaryElement bem;

    //spline::Quintic sp;
    bem.sp().init(knots);
    bem.sp().setNode();
    bem.sp().setBC(0, spline::BCType::Odd, spline::BCType::Odd);
    bem.sp().setBC(1, spline::BCType::Even, spline::BCType::Even);
    bem.sp().update();
    bem.sp().write("../resources/testBEM/spline.txt");
    bem.quadratureOrder(4);
    bem.elementOrder(2);
    bem.indexShift(0);
    bem.initialize();

    std::cin.get();
    return 0;
}