#include <cmath>
#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include "dynamicCone.hpp"
#include "numericTools/geometry2D.hpp"
#include "numericTools/legendrePoly.hpp"
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
};

int main()
{
    double c1 = -0.5;
    double b0 = 1;
    double rCut = 100.0;
    dynamicCone::PatchedCone patchedCone(c1, b0, rCut, 500, 400, 240);

    int nTotal = patchedCone.parameter().n[0] + patchedCone.parameter().n[1];

    Eigen::MatrixXd S, D, B, L, R;
    Eigen::VectorXd diagB;

    S.setZero(nTotal, nTotal);
    D.setZero(nTotal, nTotal);
    R.setZero(nTotal, nTotal);
    L.setZero(nTotal, nTotal);
    diagB.setOnes(nTotal);

    bem2D::BoundaryElement::assembleMatrix(patchedCone._bem[0], patchedCone._bem[0], S, D);
    bem2D::BoundaryElement::assembleMatrix(patchedCone._bem[0], patchedCone._bem[1], S, D);
    bem2D::BoundaryElement::assembleMatrix(patchedCone._bem[1], patchedCone._bem[0], S, D);
    bem2D::BoundaryElement::assembleMatrix(patchedCone._bem[1], patchedCone._bem[1], S, D);

    diagB = diagB * 0.5;
    diagB(patchedCone.parameter().n[0] - 1) = -D.row(patchedCone.parameter().n[0] - 1).sum();
    diagB(patchedCone.parameter().n[0]) = -D.row(patchedCone.parameter().n[0]).sum();
    for (int k = 0; k < D.rows(); k++)
        D(k, k) += diagB(k); // equation (7.80)

    std::cout << D.rowwise().sum().cwiseAbs().maxCoeff() << std::endl;

    Eigen::VectorXd rhs, lhs, p, q;
    rhs.setZero(nTotal);
    lhs.setZero(nTotal);
    p.setZero(nTotal);
    q.setZero(nTotal);
    Eigen::VectorXd errorCurvature(patchedCone.parameter().n[0]);

    for (int k = 0; k < nTotal; k++)
    {
        double r, dr, z, dz, nr, nz, ddr, ddz;
        if (k < patchedCone.parameter().n[0])
        {
            r = patchedCone._bem[0].node()[k].x;
            z = patchedCone._bem[0].node()[k].y;
            dr = patchedCone._bem[0].node()[k].dx;
            dz = patchedCone._bem[0].node()[k].dy;
            ddr = patchedCone._bem[0].node()[k].ddx;
            ddz = patchedCone._bem[0].node()[k].ddy;
            numericTools::Geometry2D::nrnz(dr, dz, nr, nz);
            rhs(k) = dPhiDn(nr, nz, r, z);
            lhs(k) = phi(r, z);
            errorCurvature(k) = numericTools::Geometry2D::curvatureAxisymmetric(r, z, dr, dz, ddr, ddz) -
                                numericTools::Geometry2D::curvatureC3Cone(r, patchedCone.parameter().rc, patchedCone.parameter().c);
        }
        else
        {
            int kk = k - patchedCone.parameter().n[0];
            r = patchedCone._bem[1].node()[kk].x;
            z = patchedCone._bem[1].node()[kk].y;
            dr = patchedCone._bem[1].node()[kk].dx;
            dz = patchedCone._bem[1].node()[kk].dy;
            numericTools::Geometry2D::nrnz(dr, dz, nr, nz);
            lhs(k) = dPhiDn(nr, nz, r, z);
            rhs(k) = phi(r, z);
        }
    }

    D.row(patchedCone.parameter().n[0] - 1) *= 0.;
    D(patchedCone.parameter().n[0] - 1, patchedCone.parameter().n[0] - 1) = 1.0;
    D(patchedCone.parameter().n[0] - 1, patchedCone.parameter().n[0]) = -1.0;
    S.row(patchedCone.parameter().n[0] - 1) *= 0.;

    //block matrix swapping in equation(7.84)
    bem2D::BoundaryElement::swapSDLR(S, D, patchedCone.parameter().n[0], L, R);

    Eigen::VectorXd answer = R.partialPivLu().solve(L * lhs);
    Eigen::VectorXd errorDirichlet(patchedCone.parameter().n[1]), errorNeumann(patchedCone.parameter().n[0]);

    for (int k = 0; k < nTotal; k++)
    {
        if (k < patchedCone.parameter().n[0])
        {
            p(k) = lhs(k);
            q(k) = answer(k);
            errorNeumann(k) = answer(k) - rhs(k);
        }
        else
        {
            p(k) = answer(k);
            q(k) = lhs(k);
            errorDirichlet(k - patchedCone.parameter().n[0]) = answer(k) - rhs(k);
        }
    }

    printf("%d\t%2.16f\t%2.16f\t%2.16f\n", nTotal, errorDirichlet.cwiseAbs().maxCoeff(), errorNeumann.cwiseAbs().maxCoeff(), errorCurvature.cwiseAbs().maxCoeff());

    // Eigen::MatrixXd knot[3];
    // knot[0] = numericTools::Geometry2D::c3Cone(
    //     patchedCone.parameter().rc, patchedCone.parameter().rCutoff,
    //     patchedCone.parameter().c,
    //     patchedCone.parameter().splineNodeNumber[0]);
    // knot[1] = numericTools::Geometry2D::circle(
    //     knot[0](knot[0].rows() - 1, 0),
    //     knot[0](knot[0].rows() - 1, 1),
    //     patchedCone.parameter().splineNodeNumber[1],
    //     numericTools::CircleType::DownWrap);
    // knot[2] = numericTools::Geometry2D::circle(
    //     knot[0](knot[0].rows() - 1, 0),
    //     knot[0](knot[0].rows() - 1, 1),
    //     patchedCone.parameter().splineNodeNumber[2],
    //     numericTools::CircleType::UpWrap);

    // patchedCone.setBoundaryElement(dynamicCone::BoundaryType::Interface, knot[0], 0);
    // patchedCone.setBoundaryElement(dynamicCone::BoundaryType::Liquid, knot[1], 0);
    // patchedCone.setBoundaryElement(dynamicCone::BoundaryType::Vacuum, knot[2], 0);

    patchedCone.printABC("ABC");

    return 0;
}