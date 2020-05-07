#include "dynamicCone.hpp"
#include "numericTools/geometry2D.hpp"
#include "numericTools/quadratureRules.hpp"

namespace dynamicCone
{

void PatchedCone::computeABC(const double c1, const double b0, double *a, double *b, double *c)
{
    a[0] = 2.7142175397111330 * c1;
    b[0] = b0;
    c[0] = -0.8604366861256783;

    a[1] = 0.1442586135181731 * a[0] * a[0] - 0.4749808176761397 * b0 * b0 - c[0];
    b[1] = -0.848581976487259 * b0 * c1; // old
    c[1] = c1;

    a[2] = -2.336057766096800 * c1 - 1.155514902883830 * c1 * c1 * c1 + 1.584211046805990 * b0 * b0 * c1;
    a[3] = -0.433421293527112 + 1.563930669354330 * c1 * c1 + 1.356140305325190 * pow(c1, 4.0) + 0.478517022152372 * b0 * b0 - 0.132076194633969 * pow(b0, 4.0) + 0.448670228959027 * b0 * b0 * c1 * c1;
    a[4] = -1.723725053118940 * c1 + 4.882389215855330 * pow(c1, 3.0) - 1.096390164067280 * pow(c1, 5.0) + 1.980228762128670 * b0 * b0 * c1 - 0.409457597247434 * pow(b0, 4.0) * c1 - 17.405016393743000 * b0 * b0 * pow(c1, 3.0);

    c[2] = 0.0;
    c[3] = -0.275783438603136 * c1 + 0.210659453420660 * b0 * b0 * c1 - 0.069483517708871 * c1 * c1 * c1;
    c[4] = -0.045843694202325 + 0.050613544746824 * b0 * b0 - 0.013969919726216 * pow(b0, 4.0) - 0.587515210204774 * c1 * c1 + 0.579439247828955 * b0 * b0 * c1 * c1 - 0.139013862957991 * pow(c1, 4.0);

    b[2] = -1.290655029188520 * b0 * c1 * c1;
    b[3] = 1.887300354735200 * b0 * c1 - 1.441629936818880 * b0 * b0 * b0 * c1 - 5.508686000286550 * b0 * c1 * c1 * c1;
    b[4] = -0.957571779297224 * b0 + 1.057203241210380 * pow(b0, 3.0) - 0.291800238214520 * pow(b0, 5.0) + 16.866940785206700 * b0 * c1 * c1 - 10.154738911572600 * pow(b0, 3.0) * c1 * c1 - 64.833275833138900 * b0 * pow(c1, 4.0);
};

double PatchedCone::farFieldVelocityPotential(double r, double z, const double (&a)[5])
{
    double phi_ = 0, phi_tmp = 0, dPhiDr_ = 0, dPhiDz_ = 0;
    numericTools::Geometry2D::zonalHarmonics(r, z, 0,
                                             numericTools::LegendrePolyType::HALF_INTEGER, numericTools::HarmonicType::Inner,
                                             phi_tmp, dPhiDr_, dPhiDz_);
    phi_ += a[0] * phi_tmp;
    numericTools::Geometry2D::zonalHarmonics(r, z, 0,
                                             numericTools::LegendrePolyType::INTEGER, numericTools::HarmonicType::Outer,
                                             phi_tmp, dPhiDr_, dPhiDz_);
    phi_ += a[1] * phi_tmp;
    // numericTools::Geometry2D::zonalHarmonics(r, z, 1,
    //                                          numericTools::LegendrePolyType::HALF_INTEGER, numericTools::HarmonicType::Outer,
    //                                          phi_tmp, dPhiDr_, dPhiDz_);
    // phi_ += a[2] * phi_tmp;
    // numericTools::Geometry2D::zonalHarmonics(r, z, 3,
    //                                          numericTools::LegendrePolyType::INTEGER, numericTools::HarmonicType::Outer,
    //                                          phi_tmp, dPhiDr_, dPhiDz_);
    // phi_ += a[3] * phi_tmp;
    // numericTools::Geometry2D::zonalHarmonics(r, z, 4,
    //                                          numericTools::LegendrePolyType::HALF_INTEGER, numericTools::HarmonicType::Outer,
    //                                          phi_tmp, dPhiDr_, dPhiDz_);
    phi_ += a[4] * phi_tmp;
    //printf("%f\n", phi_);
    return phi_;
};

double PatchedCone::farFieldElectricFlux(double nr, double nz, double r, double z, const double (&b)[5])
{
    double phi_tmp = 0, dPhiDn_ = 0, dPhiDz_tmp = 0, dPhiDr_tmp = 0;
    numericTools::Geometry2D::zonalHarmonics(r, z, 0,
                                             numericTools::LegendrePolyType::HALF_INTEGER, numericTools::HarmonicType::Inner,
                                             phi_tmp, dPhiDr_tmp, dPhiDz_tmp);
    //printf("%16.16f\t%16.16f\t%16.16f\n", phi_tmp, dPhiDr_tmp, dPhiDz_tmp);
    dPhiDn_ += b[0] * (nr * dPhiDr_tmp + nz * dPhiDz_tmp);
    numericTools::Geometry2D::zonalHarmonics(r, z, 0,
                                             numericTools::LegendrePolyType::INTEGER, numericTools::HarmonicType::Outer,
                                             phi_tmp, dPhiDr_tmp, dPhiDz_tmp);
    dPhiDn_ += b[1] * (nr * dPhiDr_tmp + nz * dPhiDz_tmp);
    numericTools::Geometry2D::zonalHarmonics(r, z, 1,
                                             numericTools::LegendrePolyType::HALF_INTEGER, numericTools::HarmonicType::Outer,
                                             phi_tmp, dPhiDr_tmp, dPhiDz_tmp);
    dPhiDn_ += b[2] * (nr * dPhiDr_tmp + nz * dPhiDz_tmp);
    numericTools::Geometry2D::zonalHarmonics(r, z, 3,
                                             numericTools::LegendrePolyType::INTEGER, numericTools::HarmonicType::Outer,
                                             phi_tmp, dPhiDr_tmp, dPhiDz_tmp);
    dPhiDn_ += b[3] * (nr * dPhiDr_tmp + nz * dPhiDz_tmp);
    numericTools::Geometry2D::zonalHarmonics(r, z, 4,
                                             numericTools::LegendrePolyType::HALF_INTEGER, numericTools::HarmonicType::Outer,
                                             phi_tmp, dPhiDr_tmp, dPhiDz_tmp);
    dPhiDn_ += b[4] * (nr * dPhiDr_tmp + nz * dPhiDz_tmp);
    return dPhiDn_;
};

Eigen::VectorXd PatchedCone::setVacuumBC(const bem2D::BoundaryElement &interface, const bem2D::BoundaryElement &vacuum) const
{
    int nCone = interface.node().size();
    int nPatch = vacuum.node().size();
    int nTotal = nCone + nPatch;
    Eigen::VectorXd vacuumBC;
    vacuumBC.setZero(nTotal);
    for (int k = 0; k < nTotal; k++)
    { // vacuum domain boundary  value problem (6.59)
        double r, z, dr, dz, nr, nz;

        if (k < nCone)
        {
            r = interface.node()[k].x;
            z = interface.node()[k].y;
            //dr = interface.node()[k].dx;
            //dz = interface.node()[k].dy;
            //vacuumBC(k) = 0.0;

            vacuumBC(k) = farFieldVelocityPotential(r, z, parameter().b);
            printf("%12.12f\n", vacuumBC(k));
        }

        else
        {
            r = vacuum.node()[k - nCone].x;
            z = vacuum.node()[k - nCone].y;
            dr = vacuum.node()[k - nCone].dx;
            dz = vacuum.node()[k - nCone].dy;
            numericTools::Geometry2D::nrnz(dr, dz, nr, nz);
            vacuumBC(k) = farFieldElectricFlux(nr, nz, r, z, parameter().b);
        }
    }
    return vacuumBC;
};

Eigen::VectorXd PatchedCone::setLiquidBC(const bem2D::BoundaryElement &interface, const bem2D::BoundaryElement &liquid) const
{

    int nCone = interface.node().size();
    int nPatch = liquid.node().size();
    int nTotal = nCone + nPatch;
    Eigen::VectorXd liquidBC;

    liquidBC.setZero(nTotal);
    for (int k = 0; k < nTotal; k++)
    { // liquid domain boundary value problem (6.58)
        double r, z, dr, dz, nr, nz;
        if (k < nCone)
        {
            r = liquid.node()[k - nCone].x;
            z = liquid.node()[k - nCone].y;
            dr = liquid.node()[k - nCone].dx;
            dz = liquid.node()[k - nCone].dy;
            numericTools::Geometry2D::nrnz(dr, dz, nr, nz);
            liquidBC(k) = -2. / 3. * (nr * r + nz * z);
        }
        else
        {
            r = liquid.node()[k - nCone].x;
            z = liquid.node()[k - nCone].y;
            liquidBC(k) = farFieldVelocityPotential(r, z, parameter().a);
        }
    }
    return liquidBC;
};

void PatchedCone::solveBVP(
    const bem2D::BoundaryElement &bem0, const Eigen::VectorXd &given, BoundaryType type,
    const Eigen::MatrixXd &SRaw, const Eigen::MatrixXd &DRaw, Eigen::VectorXd &p, Eigen::VectorXd &q)
{
    int nTotal = 0;
    const int &n0 = parameter().n[0];
    if (type == BoundaryType::Liquid)
        nTotal = n0 + parameter().n[1];
    else if (type == BoundaryType::Vacuum)
        nTotal = n0 + parameter().n[2];

    Eigen::MatrixXd L, R, D, S;
    Eigen::VectorXd diagB;
    R.setZero(nTotal, nTotal);
    L.setZero(nTotal, nTotal);
    D.setZero(nTotal, nTotal);
    S.setZero(nTotal, nTotal);
    diagB.setZero(nTotal);
    diagB.setOnes();

    if (type == BoundaryType::Liquid)
        D = DRaw;
    else if (type == BoundaryType::Vacuum)
        D = -DRaw;
    S = SRaw;

    diagB = diagB * 0.5;
    diagB(n0 - 1) = -D.row(n0 - 1).sum();
    diagB(n0) = -D.row(n0).sum();
    for (int k = 0; k < D.rows(); k++)
        D(k, k) += diagB(k); // equation (7.80)

    D.row(n0 - 1) *= 0.;
    D(n0 - 1, n0 - 1) = 1.0;
    D(n0 - 1, n0) = -1.0;
    S.row(n0 - 1) *= 0.;

    bem2D::BoundaryElement::swapSDLR(S, D, n0, L, R);

    Eigen::VectorXd answer;
    answer.setZero(nTotal);
    p.setZero(nTotal);
    q.setZero(nTotal);

    if (type == BoundaryType::Liquid)
    {
        answer = L.partialPivLu().solve(R * given);
        rearrangeNodalVectors(answer, given, n0, p, q);
    }
    else if (type == BoundaryType::Vacuum)
    {
        answer = R.partialPivLu().solve(L * given);
        rearrangeNodalVectors(given, answer, n0, p, q);
    }
}

void PatchedCone::rearrangeNodalVectors(const Eigen::VectorXd &u0, const Eigen::VectorXd &u1, int n0, Eigen::VectorXd &v0, Eigen::VectorXd &v1)
{
    int nTotal = v0.size();
    for (int k = 0; k < nTotal; k++)
    {
        if (k < n0)
        {
            v0(k) = u0(k);
            v1(k) = u1(k);
        }
        else
        {
            v0(k) = u1(k);
            v1(k) = u0(k);
        }
    }
};

double PatchedCone::farFieldShape(double r, const double (&c)[5])
{
    return c[0] * r + c[1] / sqrt(r) + c[3] / pow(r, 3.5) + c[4] / pow(r, 5.0);
};

double PatchedCone::farFieldShapeDr(double r, const double (&c)[5])
{
    return c[0] - (5.0 * c[4]) / pow(r, 6.0) - (7 * c[3]) / (2. * pow(r, 4.5)) - c[1] / (2. * pow(r, 1.5));
};

double PatchedCone::farFieldShapeArc(double r0, double r1, const double (&c)[5])
{
    const int nqd = 16;
    const double *ab;
    const double *w;
    ab = numericTools::QuadratureRules::abascissa(nqd, numericTools::QuadratureType::GAUSS_LEGENDRE);
    w = numericTools::QuadratureRules::weight(nqd, numericTools::QuadratureType::GAUSS_LEGENDRE);

    double arc = 0;

    for (int k = 0; k < nqd; k++)
    {
        double r = r0 + (r1 - r0) * ab[k];
        double dShapeDr = farFieldShapeDr(r, c);
        arc += w[k] * (r1 - r0) * sqrt(1 + dShapeDr * dShapeDr);
    }

    return arc;
};

const Eigen::MatrixXd PatchedCone::farFieldShapePadding(double rEnd, double rIncrement) const // STILL WORKING ON
{
    Eigen::MatrixXd pad;
    pad.resize(5, 2);
    pad.setZero();
    pad(0, 0) = 0;

    pad(0, 1) = farFieldVelocityPotential(rEnd, farFieldShape(rEnd, parameter().c), parameter().a);

    for (int k = 1; k < 5; k++)
    {
        double rNow = rEnd + k * rIncrement;
        double zNow = farFieldShape(rNow, parameter().c);
        double rLast = rEnd + (k - 1) * rIncrement;

        pad(k, 0) = pad(k - 1, 0) + farFieldShapeArc(rLast, rNow, parameter().c);
        pad(k, 1) = farFieldVelocityPotential(rNow, zNow, parameter().a);
    }
    return pad;
}

PatchedCone::PatchedCone(double c1, double b0, double rCut, int n0, int n1, int n2)
{
    computeABC(c1, b0, _parameter.a, _parameter.b, _parameter.c);
    _parameter.rCutoff = rCut;
    _parameter.splineNodeNumber[0] = n0;
    _parameter.splineNodeNumber[1] = n1;
    _parameter.splineNodeNumber[2] = n2;

    Eigen::MatrixXd knot[3];
    knot[0] = numericTools::Geometry2D::c3Cone(
        parameter().rc, parameter().rCutoff, parameter().c, parameter().splineNodeNumber[0]);
    knot[1] = numericTools::Geometry2D::circle(
        knot[0](knot[0].rows() - 1, 0),
        knot[0](knot[0].rows() - 1, 1),
        parameter().splineNodeNumber[1], numericTools::CircleType::DownWrap);
    knot[2] = numericTools::Geometry2D::circle(
        knot[0](knot[0].rows() - 1, 0),
        knot[0](knot[0].rows() - 1, 1),
        parameter().splineNodeNumber[2],
        numericTools::CircleType::UpWrap);

    setBoundaryElement(BoundaryType::Interface, knot[0], 0);
    setBoundaryElement(BoundaryType::Liquid, knot[1], 0);
    setBoundaryElement(BoundaryType::Vacuum, knot[2], 0);

    for (int k = 0; k <= 2; k++)
        _parameter.n[k] = _bem[k].node().size();
    _bem[0].indexShift(0);
    _bem[1].indexShift(parameter().n[0]);
    _bem[2].indexShift(parameter().n[0]);
};

void PatchedCone::printABC(const std::string &name) const
{
    const double *a = parameter().a;
    const double *b = parameter().b;
    const double *c = parameter().c;
    printf("%+18.15f\t%+18.15f\t%+18.15f\t%+18.15f\t%+18.15f\t\n", a[0], a[1], a[2], a[3], a[4]);
    printf("%+18.15f\t%+18.15f\t%+18.15f\t%+18.15f\t%+18.15f\t\n", b[0], b[1], b[2], b[3], b[4]);
    printf("%+18.15f\t%+18.15f\t%+18.15f\t%+18.15f\t%+18.15f\t\n", c[0], c[1], c[2], c[3], c[4]);

    //Eigen::MatrixXd outputMat(3, 5);
    //outputMat.row(0) << a[0], a[1], a[2], a[3], a[4];
    //outputMat.row(1) << b[0], b[1], b[2], b[3], b[4];
    //outputMat.row(2) << c[0], c[1], c[2], c[3], c[4];

    //Bem::writeMatrix(name, outputMat);
    //    printf("wtf %d\t%d\t%d\n", _bem[0].sp().node().rows(), _bem[1].sp().node().rows(), _bem[2].sp().node().rows());
    //printf("wtf %d\t%d\t%d\n", parameter().n[0], parameter().n[1], parameter().n[2]);
}

const PatchedCone::Parameter &PatchedCone::parameter() const
{
    return _parameter;
}; // 1st coordinate

void PatchedCone::setBoundaryElement(BoundaryType type, const Eigen::MatrixX2d &xy, int shift)
{
    switch (type)
    {
    case BoundaryType::Interface:
    {
        double dx, dy, ddx, ddy;
        // auto knot = numericTools::Geometry2D::c3Cone(parameter().rc, parameter().rCutoff,
        //                                              parameter().c, parameter().n[0]);
        _bem[0].sp().init(xy);
        numericTools::Geometry2D::estimateDerivativeConeEnd(xy, parameter().c, dx, ddx, dy, ddy);
        _bem[0].sp().setNode();
        _bem[0].sp().setBC(0, spline::BCType::Odd, spline::BCType::Mix, 0, 0, dx, ddx);
        _bem[0].sp().setBC(1, spline::BCType::Even, spline::BCType::Mix, 0, 0, dy, ddy);
        _bem[0].sp().update();

        _bem[0].quadratureOrder(16);
        _bem[0].elementOrder(2);
        _bem[0].initialize();

        break;
    }
    case BoundaryType::Vacuum:
    {
        double dx, dy, ddx, ddy;
        _bem[2].sp().init(xy);
        numericTools::Geometry2D::estimateDerivativeCircleBegin(xy, dx, ddx, dy, ddy);
        _bem[2].sp().setNode();
        _bem[2].sp().setBC(0, spline::BCType::Mix, spline::BCType::Odd, dx, ddx, 0, 0);
        _bem[2].sp().setBC(1, spline::BCType::Mix, spline::BCType::Even, dy, ddy, 0, 0);
        _bem[2].sp().update();

        _bem[2].quadratureOrder(16);
        _bem[2].elementOrder(2);
        _bem[2].initialize();
        break;
    }
    case BoundaryType::Liquid:
    {
        double dx, dy, ddx, ddy;
        _bem[1].sp().init(xy);
        numericTools::Geometry2D::estimateDerivativeCircleBegin(xy, dx, ddx, dy, ddy);
        _bem[1].sp().setNode();
        _bem[1].sp().setBC(0, spline::BCType::Mix, spline::BCType::Odd, dx, ddx, 0, 0);
        _bem[1].sp().setBC(1, spline::BCType::Mix, spline::BCType::Even, dy, ddy, 0, 0);
        _bem[1].sp().update();

        _bem[1].quadratureOrder(16);
        _bem[1].elementOrder(2);
        _bem[1].initialize();
        break;
    }
    default:
        break;
    }
};
} // namespace dynamicCone