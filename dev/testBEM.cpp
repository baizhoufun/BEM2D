#include <cmath>
#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include "boundaryElement.hpp"

Eigen::IOFormat fmt(Eigen::FullPrecision, 0, "\t", "\n", "", "", "", "");
std::ofstream file;

int main()
{
    Eigen::size_t N = 16;
    int dim = 2;
    Eigen::MatrixXd knots;
    knots.resize(N + 1, dim);
    for (int i = 0; i < N + 1; i++)
    {
        double theta = (double)i / N * M_PI;
        knots(i, 0) = sin(theta) * (1 + 0.25 * cos(8 * theta - M_PI)); // x (or r) coord
        knots(i, 1) = cos(theta) * (1 + 0.25 * cos(8 * theta - M_PI)); // y (or z) coord
    }

    bem2D::BoundaryElement bem;

    //spline::Quintic sp;
    bem.sp().init(knots);
    bem.sp().setNode();

    bem.sp().setBC(0, spline::BCType::Odd, spline::BCType::Odd);
    bem.sp().setBC(1, spline::BCType::Even, spline::BCType::Even);
    bem.sp().update();
    bem.sp().write("../resources/testBEM/spline.txt");
    bem.initialize();

    std::cin.get();
    return 0;
}