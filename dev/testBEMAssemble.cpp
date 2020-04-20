#include <cmath>
#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include "boundaryElement.hpp"
#include "numericTools/geometry2D.hpp"

Eigen::IOFormat fmt(Eigen::FullPrecision, 0, "\t", "\n", "", "", "", "");

auto phi = [](double r, double z) {
    double phi_, tmp1, tmp2;
    numericTools::Geometry2D::zonalHarmonics(r, z, 2,
                                             numericTools::LegendrePolyType::INTEGER, numericTools::HarmonicType::Inner,
                                             phi_, tmp1, tmp2);
    return phi_;
}; // eqn 7.88: phi

auto dPhiDn = [](double nr, double nz, double r, double z) {
    double tmp0, dPhiDr_, dPhiDz_;
    numericTools::Geometry2D::zonalHarmonics(r, z, 2,
                                             numericTools::LegendrePolyType::INTEGER, numericTools::HarmonicType::Inner,
                                             tmp0, dPhiDr_, dPhiDz_);
    return dPhiDr_ * nr + dPhiDz_ * nz; // nr * (-r) + nz * (2. * z);
};                                      //eqn 7.88: dphi/dn

int main()
{
    int N = 128;
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

    bem.sp().write("../resources/testBEMAssemble/spline.txt");
    bem.quadratureOrder(6);
    bem.elementOrder(2);
    bem.indexShift(0);
    bem.initialize();

    int n0 = bem.node().size();
    Eigen::MatrixXd S, D, B;
    S.setZero(n0, n0);
    B.setZero(n0, n0);
    D.setZero(n0, n0);

    B.setIdentity();
    B *= 0.5;

    Eigen::VectorXd rhs, lhs, errorCurvature;
    rhs.setZero(n0);
    lhs.setZero(n0);

    bem2D::BoundaryElement::assembleMatrix(bem, bem, S, D);
    D = D + B; // equation (7.80)

    for (int k = 0; k < n0; k++)
    {
        double r, dr, z, dz, nr, nz;

        r = bem.node()[k].x;
        z = bem.node()[k].y;
        dr = bem.node()[k].dx;
        dz = bem.node()[k].dy;
        numericTools::Geometry2D::nrnz(dr, dz, nr, nz);
        rhs(k) = dPhiDn(nr, nz, r, z);
        lhs(k) = phi(r, z);
        //    errorCurvature(k) = Bem::Geometry::curvBenchmark(acos(z / sqrt(r * r + z * z))) - Bem::Geometry::totalCurv(r, z, dr, dz, ddr, ddz);
    }
    std::cout << "solving ..." << Eigen::nbThreads() << "\n";
    Eigen::VectorXd answer = S.partialPivLu().solve(D * lhs);
    std::cout << (answer - rhs).cwiseAbs().maxCoeff() << "\n";
    std::cout << "=====================\n";
    double drp = 0.1, dzp = 0.1;
    for (int ix = 1; ix < 10; ix++)
    {
        for (int iy = 1; iy < 10; iy++)
        {
            double rp = drp * ix, zp = dzp * iy;
            std::cout << bem2D::BoundaryElement::interiorField(rp, zp, bem, answer, lhs) - phi(rp, zp) << "\n";
        }
    }

    //std::cout <<  << "\n";

    std::cin.get();
    return 0;
}
