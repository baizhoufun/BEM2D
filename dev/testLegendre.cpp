#include <cmath>
#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include "numericTools/legendrePoly.hpp"

using namespace numericTools;

Eigen::IOFormat fmt(Eigen::FullPrecision, 0, "\t", "\n", "", "", "", "");
std::ofstream file;

int main()
{
    int mNum = 100;
    Eigen::MatrixXd legendrePoly;
    legendrePoly.resize(mNum, 12);

    for (int i = 0; i < mNum; i++)
    {
        double x = 2 * ((double)(i + 1) / (mNum)-0.5);
        legendrePoly(i, 0) = x;
        legendrePoly(i, 0 + 1) = LegendrePoly::legendreP(0, x, LegendrePolyType::INTEGER);
        legendrePoly(i, 2 + 1) = LegendrePoly::legendreP(1, x, LegendrePolyType::INTEGER);
        legendrePoly(i, 4 + 1) = LegendrePoly::legendreP(2, x, LegendrePolyType::INTEGER);
        legendrePoly(i, 6 + 1) = LegendrePoly::legendreP(3, x, LegendrePolyType::INTEGER);
        legendrePoly(i, 8 + 1) = LegendrePoly::legendreP(4, x, LegendrePolyType::INTEGER);
        legendrePoly(i, 10 + 1) = LegendrePoly::legendreP(5, x, LegendrePolyType::INTEGER);
        legendrePoly(i, 1 + 1) = LegendrePoly::legendreP(0, x, LegendrePolyType::HALF_INTEGER);
        legendrePoly(i, 3 + 1) = LegendrePoly::legendreP(1, x, LegendrePolyType::HALF_INTEGER);
        legendrePoly(i, 5 + 1) = LegendrePoly::legendreP(2, x, LegendrePolyType::HALF_INTEGER);
        legendrePoly(i, 7 + 1) = LegendrePoly::legendreP(3, x, LegendrePolyType::HALF_INTEGER);
        legendrePoly(i, 9 + 1) = LegendrePoly::legendreP(4, x, LegendrePolyType::HALF_INTEGER);
    }
    file.open("../resources/testNumericTools/legendrePoly.txt");
    file << legendrePoly.format(fmt);
    file.close();

    return 0;
}