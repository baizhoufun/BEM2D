#include <cmath>
#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include "boundaryElement.hpp"
#include "numericTools/geometry2D.hpp"

Eigen::IOFormat fmt(Eigen::FullPrecision, 0, "\t", "\n", "", "", "", "");

auto phi = [](double r, double z) {
    double phi_, dPhiDr_, dPhiDz_;
    numericTools::Geometry2D::zonalHarmonics(r, -z, 1,
                                             numericTools::LegendrePolyType::HALF_INTEGER, numericTools::HarmonicType::Inner,
                                             phi_, dPhiDr_, dPhiDz_);
    //phi_ = -0.5 * r * r + z * z;
    return phi_;
}; // eqn 7.88: phi

auto dPhiDn = [](double nr, double nz, double r, double z) {
    double phi_, dPhiDr_, dPhiDz_, dPhiDn_;
    numericTools::Geometry2D::zonalHarmonics(r, -z, 1,
                                             numericTools::LegendrePolyType::HALF_INTEGER, numericTools::HarmonicType::Inner,
                                             phi_, dPhiDr_, dPhiDz_);
    dPhiDn_ = dPhiDr_ * nr - dPhiDz_ * nz;
    //dPhiDn_ = nr * (-r) + nz * (2. * z);
    return dPhiDn_;
}; //eqn 7.88: dphi/dn

int main()
{
    int globalMin = 5, globalMax = 11; //17;

    for (int global = globalMin; global < globalMax; global++)
    {

        int N = (int)pow(sqrt(2.0), global) + 1;
        double c[5] = {-0.86, -2, -0.0, -3.0, -2.0};
        double dx, dy, ddx, ddy;

        Eigen::MatrixXd knot[2];
        knot[0] = numericTools::Geometry2D::c3Cone(35, 50, c, N);
        knot[1] = numericTools::Geometry2D::circle(knot[0](knot[0].rows() - 1, 0), knot[0](knot[0].rows() - 1, 1), N / 1.5 + 1,
                                                   numericTools::CircleType::DownWrap);

        bem2D::BoundaryElement bem[2];

        bem[0].sp().init(knot[0]);
        numericTools::Geometry2D::estimateDerivativeConeEnd(knot[0], c, dx, ddx, dy, ddy);
        bem[0].sp().setNode();
        bem[0].sp().setBC(0, spline::BCType::Odd, spline::BCType::Mix, 0, 0, dx, ddx);
        bem[0].sp().setBC(1, spline::BCType::Even, spline::BCType::Mix, 0, 0, dy, ddy);
        // bem[0].sp()
        //     .setBC(0, spline::BCType::Odd, spline::BCType::Mix, 0, 0,
        //            bem[0].sp().estimateDerivative(0, spline::BCLocation::End, 1),
        //            bem[0].sp().estimateDerivative(0, spline::BCLocation::End, 2));

        // bem[0].sp().setBC(1, spline::BCType::Even, spline::BCType::Mix, 0, 0,
        //                   c[0] * bem[0].sp().estimateDerivative(0, spline::BCLocation::End, 1),
        //                   0);
        bem[0].sp().update();

        bem[1].sp().init(knot[1]);
        numericTools::Geometry2D::estimateDerivativeCircleBegin(knot[1], dx, ddx, dy, ddy);
        bem[1].sp().setNode();
        bem[1].sp().setBC(0, spline::BCType::Mix, spline::BCType::Odd, dx, ddx, 0, 0);
        bem[1].sp().setBC(1, spline::BCType::Mix, spline::BCType::Even, dy, ddy, 0, 0);
        bem[1].sp().update();

        for (int g = 0; g < -1; g++)
        {
            bem[0].sp().setNode(bem[0].sp().arcIncrement());
            bem[0].sp().update();
            bem[1].sp().setNode(bem[1].sp().arcIncrement());
            bem[1].sp().update();
        }

        bem[0].sp().write("../resources/testBEMCone/spline0.txt");
        bem[1].sp().write("../resources/testBEMCone/spline1.txt");

        bem[0].quadratureOrder(20);
        bem[0].elementOrder(2);
        bem[0].initialize();

        bem[1].quadratureOrder(20);
        bem[1].elementOrder(2);
        bem[1].initialize();

        int n[2];
        n[0] = bem[0].node().size();
        n[1] = bem[1].node().size();
        bem[0].indexShift(0);
        bem[1].indexShift(n[0]);
        int nTotal = n[0] + n[1];

        Eigen::MatrixXd S, D, B, L, R;

        S.setZero(nTotal, nTotal);
        B.setZero(nTotal, nTotal);
        D.setZero(nTotal, nTotal);
        R.setZero(nTotal, nTotal);
        L.setZero(nTotal, nTotal);

        Eigen::VectorXd rhs, lhs, p, q;
        rhs.setZero(nTotal);
        lhs.setZero(nTotal);
        p.setZero(nTotal);
        q.setZero(nTotal);

        bem2D::BoundaryElement::assembleMatrix(bem[0], bem[0], S, D);
        bem2D::BoundaryElement::assembleMatrix(bem[0], bem[1], S, D);
        bem2D::BoundaryElement::assembleMatrix(bem[1], bem[0], S, D);
        bem2D::BoundaryElement::assembleMatrix(bem[1], bem[1], S, D);

        B.setIdentity();
        B *= 0.5;
        B(n[0] - 1, n[0] - 1) = -D.row(n[0] - 1).sum();
        B(n[0], n[0]) = -D.row(n[0]).sum();
        D = D + B; // equation (7.80)

        for (int k = 0; k < nTotal; k++)
        {

            double r, dr, z, dz, nr, nz;
            if (k < n[0])
            {
                r = bem[0].node()[k].x;
                z = bem[0].node()[k].y;
                dr = bem[0].node()[k].dx;
                dz = bem[0].node()[k].dy;
                numericTools::Geometry2D::nrnz(dr, dz, nr, nz);
                rhs(k) = dPhiDn(nr, nz, r, z);
                lhs(k) = phi(r, z);
                //errorCurvature(k) = Bem::Geometry::curvBenchmark(acos(z / sqrt(r * r + z * z))) - Bem::Geometry::totalCurv(r, z, dr, dz, ddr, ddz);
            }
            else
            {
                int kk = k - n[0];
                r = bem[1].node()[kk].x;
                z = bem[1].node()[kk].y;
                dr = bem[1].node()[kk].dx;
                dz = bem[1].node()[kk].dy;
                numericTools::Geometry2D::nrnz(dr, dz, nr, nz);
                lhs(k) = dPhiDn(nr, nz, r, z);
                rhs(k) = phi(r, z);
            }
        }

        D.row(n[0] - 1) *= 0.;
        D(n[0] - 1, n[0] - 1) = 1.0;
        D(n[0] - 1, n[0]) = -1.0;
        S.row(n[0] - 1) *= 0.;

        //block matrix swapping in equation(7.84)
        bem2D::BoundaryElement::swapSDLR(S, D, n[0], L, R);

        Eigen::VectorXd answer = R.partialPivLu().solve(L * lhs);
        Eigen::VectorXd errorDirichlet(n[1]), errorNeumann(n[0]);

        for (int k = 0; k < nTotal; k++)
        {
            if (k < n[0])
            {
                p(k) = lhs(k);
                q(k) = answer(k);
                errorNeumann(k) = answer(k) - rhs(k);
            }
            else
            {
                p(k) = answer(k);
                q(k) = lhs(k);
                errorDirichlet(k - n[0]) = answer(k) - rhs(k);
            }
        }

        std::cout << nTotal << "\t";
        std::cout << errorDirichlet.cwiseAbs().maxCoeff() << "\t";
        std::cout << errorNeumann.cwiseAbs().maxCoeff() << std::endl;

        // std::cout << (answer - rhs).cwiseAbs().maxCoeff() << "\n";
        // std::cout << "=====================\n";
        // double drp = 0.1, dzp = 0.1;
        // for (int ix = 1; ix < 10; ix++)
        // {
        //     for (int iy = 1; iy < 10; iy++)
        //     {
        //         double rp = drp * ix, zp = dzp * iy;
        //         std::cout << bem2D::BoundaryElement::interiorField(rp, zp, bem, answer, lhs) - phi(rp, zp) << "\n";
        //     }
        // }
    }
    std::cin.get();
    return 0;
}
