#include <cmath>
#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include "boundaryElement.hpp"
#include "numericTools/geometry2D.hpp"
#include "io/ioEigen.hpp"
#include "io/utilities.hpp"

Eigen::IOFormat fmt(Eigen::FullPrecision, 0, "\t", "\n", "", "", "", "");

auto phi = [](double r, double z) {
    double phi_ = 0, phi_tmp = 0, dPhiDr_, dPhiDz_;
    numericTools::Geometry2D::zonalHarmonics(r, -z, 0,
                                             numericTools::LegendrePolyType::HALF_INTEGER, numericTools::HarmonicType::Inner,
                                             phi_tmp, dPhiDr_, dPhiDz_);

    phi_ += phi_tmp;
    numericTools::Geometry2D::zonalHarmonics(r, -z, 1,
                                             numericTools::LegendrePolyType::INTEGER, numericTools::HarmonicType::Outer,
                                             phi_tmp, dPhiDr_, dPhiDz_);

    phi_ += 500.0 * phi_tmp;

    return phi_; //phi_ = -0.5 * r * r + z * z;
};               // eqn 7.88: phi

auto dPhiDn = [](double nr, double nz, double r, double z) {
    double phi_, dPhiDr_ = 0, dPhiDz_ = 0, dPhiDn_ = 0, dPhiDz_tmp = 0, dPhiDr_tmp = 0;
    numericTools::Geometry2D::zonalHarmonics(r, -z, 0,
                                             numericTools::LegendrePolyType::HALF_INTEGER, numericTools::HarmonicType::Inner,
                                             phi_, dPhiDr_tmp, dPhiDz_tmp);
    dPhiDr_ += dPhiDr_tmp;
    dPhiDz_ += dPhiDz_tmp;
    numericTools::Geometry2D::zonalHarmonics(r, -z, 1,
                                             numericTools::LegendrePolyType::INTEGER, numericTools::HarmonicType::Outer,
                                             phi_, dPhiDr_tmp, dPhiDz_tmp);
    dPhiDr_ += 500.0 * dPhiDr_tmp;
    dPhiDz_ += 500.0 * dPhiDz_tmp;
    dPhiDn_ = dPhiDr_ * nr - dPhiDz_ * nz;

    return dPhiDn_; //dPhiDn_ = nr * (-r) + nz * (2. * z);
};                  //eqn 7.88: dphi/dn

double c[5] = {-0.86, -2, -0.0, -3.0, -2.0};
double rc = 10;

void test(int N, int order, const std::string tag)
{
    double dx, dy, ddx, ddy;

    Eigen::MatrixXd knot[2];
    knot[0] = numericTools::Geometry2D::c3Cone(rc, 50, c, N);
    knot[1] = numericTools::Geometry2D::circle(knot[0](knot[0].rows() - 1, 0), knot[0](knot[0].rows() - 1, 1), N * 0.86 + 1,
                                               numericTools::CircleType::DownWrap);

    bem2D::BoundaryElement bem[2];

    bem[0].sp().init(knot[0]);
    numericTools::Geometry2D::estimateDerivativeConeEnd(knot[0], c, dx, ddx, dy, ddy);
    bem[0].sp().setNode();
    bem[0].sp().setBC(0, spline::BCType::Odd, spline::BCType::Mix, 0, 0, dx, ddx);
    bem[0].sp().setBC(1, spline::BCType::Even, spline::BCType::Mix, 0, 0, dy, ddy);
    bem[0].sp().update();

    bem[1].sp().init(knot[1]);
    numericTools::Geometry2D::estimateDerivativeCircleBegin(knot[1], dx, ddx, dy, ddy);
    bem[1].sp().setNode();
    bem[1].sp().setBC(0, spline::BCType::Mix, spline::BCType::Odd, dx, ddx, 0, 0);
    bem[1].sp().setBC(1, spline::BCType::Mix, spline::BCType::Even, dy, ddy, 0, 0);
    bem[1].sp().update();

    for (int g = 0; g < 20; g++)
    {
        bem[0].sp().setNode(bem[0].sp().arcIncrement());
        bem[0].sp().update();
        bem[1].sp().setNode(bem[1].sp().arcIncrement());
        bem[1].sp().update();
    }

    bem[0].sp().write("../resources/testBEMCone/spline0_" + tag + ".txt");
    bem[1].sp().write("../resources/testBEMCone/spline1_" + tag + ".txt");

    bem[0].quadratureOrder(16);
    bem[0].elementOrder(order);
    bem[0].initialize();

    bem[1].quadratureOrder(16);
    bem[1].elementOrder(order);
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

    Eigen::VectorXd errorCurvature(n[0]);

    for (int k = 0; k < nTotal; k++)
    {

        double r, dr, z, dz, nr, nz, ddr, ddz;
        if (k < n[0])
        {
            r = bem[0].node()[k].x;
            z = bem[0].node()[k].y;
            dr = bem[0].node()[k].dx;
            dz = bem[0].node()[k].dy;
            ddr = bem[0].node()[k].ddx;
            ddz = bem[0].node()[k].ddy;
            numericTools::Geometry2D::nrnz(dr, dz, nr, nz);
            rhs(k) = dPhiDn(nr, nz, r, z);
            lhs(k) = phi(r, z);
            errorCurvature(k) = numericTools::Geometry2D::curvatureAxisymmetric(r, z, dr, dz, ddr, ddz) - numericTools::Geometry2D::curvatureC3Cone(r, rc, c);
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

    printf("%d\t%2.16f\t%2.16f\t%2.16f\n", nTotal, errorDirichlet.cwiseAbs().maxCoeff(), errorNeumann.cwiseAbs().maxCoeff(), errorCurvature.cwiseAbs().maxCoeff());

    bem2D::BoundaryElement::scan(bem[0], bem[1], q, p, "../resources/testBEMCone/fieldPointCoord.txt", "../resources/testBEMCone/scan_" + tag + ".txt");
};

int main()
{
    int globalMin = 12, globalMax = 34; //17;

    for (int global = globalMin; global < globalMax; global++)
    {
        int N = (int)pow(pow(2.0, 0.25), global) + 1;
        const auto tag = bem2D::io::Utilities::padZero(global - globalMin, 10);
        //std::cout << bem2D::io::Utilities::padZero(N, 100) << ":::::::\n";
        test(N, 2, tag);
    }

    std::cin.get();
    std::cout << "\n";
    return 0;
}