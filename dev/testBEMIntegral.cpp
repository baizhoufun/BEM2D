#include <cmath>
#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include "bem2D/boundaryElement.hpp"
#include "numericTools/geometry2D.hpp"

Eigen::IOFormat fmt(Eigen::FullPrecision, 0, "\t", "\n", "", "", "", "");

int main()
{

    int N = 32;
    Eigen::MatrixXd knots = numericTools::Geometry2D::benchmarkShape(0, M_PI, N + 1);
    bem2D::BoundaryElement bem;

    bem.sp().init(knots);
    bem.sp().setNode();
    bem.sp().setBC(0, spline::BCType::Odd, spline::BCType::Odd);
    bem.sp().setBC(1, spline::BCType::Even, spline::BCType::Even);
    bem.sp().update();

    for (int g = 0; g < 10; g++)
    {
        bem.sp().setNode(bem.sp().arcIncrement());
        bem.sp().update();
    }

    bem.sp().write("../resources/testBEMIntegral/spline.txt");
    bem.quadratureOrder(16);
    bem.elementOrder(2);
    bem.indexShift(0);
    bem.initialize();

    double rp = 3.0, zp = 4.0;
    int iElement = 8;
    Eigen::MatrixXd output(6, 4);
    output.setZero();

    output.col(0) = bem.axis(zp, iElement);
    output.col(1) = bem.regular(rp, zp, iElement);
    output.col(2) = bem.singular(bem.element()[iElement].t()[0], iElement);
    output.col(3) = bem.singular(bem.element()[iElement].t()[1], iElement);

    std::cout << output.format(fmt) << "\n";
    printf("tau = %2.16f\n", bem.element()[iElement].t()[1]);
    std::cin.get();
    return 0;
}
