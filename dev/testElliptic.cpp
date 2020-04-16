#include <cmath>
#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include "ellipticIntegral.hpp"

using namespace numericTools;

Eigen::IOFormat fmt(Eigen::FullPrecision, 0, "\t", "\n", "", "", "", "");
std::ofstream file;

int main()
{
    int mNum = 100;
    Eigen::MatrixX3d ellipticIntegral(mNum, 3);
    for (int i = 0; i < mNum; i++)
    {
        double m = (double)i / (mNum);
        double Km, Em;
        EllipticIntegral::ellipticKE(m, Km, Em);
        ellipticIntegral(i, 0) = m;
        ellipticIntegral(i, 1) = Km;
        ellipticIntegral(i, 2) = Em;
    }
    file.open("../resources/testNumericTools/ellipticIntegral.txt");
    file << ellipticIntegral.format(fmt);
    file.close();

    return 0;
}