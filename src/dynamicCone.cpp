#include "dynamicCone.hpp"
#include "numericTools/geometry2D.hpp"
#include "numericTools/quadratureRules.hpp"

void DynamicCone::computeABC(const double c1, const double b0, double *a, double *b, double *c)
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

double DynamicCone::farFieldVelocityPotential(double r, double z, const double (&a)[5])
{
    double phi_ = 0, phi_tmp = 0, dPhiDr_, dPhiDz_;
    numericTools::Geometry2D::zonalHarmonics(r, -z, 0,
                                             numericTools::LegendrePolyType::HALF_INTEGER, numericTools::HarmonicType::Inner,
                                             phi_tmp, dPhiDr_, dPhiDz_);
    phi_ += a[0] * phi_tmp;
    numericTools::Geometry2D::zonalHarmonics(r, -z, 0,
                                             numericTools::LegendrePolyType::INTEGER, numericTools::HarmonicType::Outer,
                                             phi_tmp, dPhiDr_, dPhiDz_);
    phi_ += a[1] * phi_tmp;
    numericTools::Geometry2D::zonalHarmonics(r, -z, 1,
                                             numericTools::LegendrePolyType::HALF_INTEGER, numericTools::HarmonicType::Outer,
                                             phi_tmp, dPhiDr_, dPhiDz_);
    phi_ += a[2] * phi_tmp;
    numericTools::Geometry2D::zonalHarmonics(r, -z, 3,
                                             numericTools::LegendrePolyType::INTEGER, numericTools::HarmonicType::Outer,
                                             phi_tmp, dPhiDr_, dPhiDz_);
    phi_ += a[3] * phi_tmp;
    numericTools::Geometry2D::zonalHarmonics(r, -z, 4,
                                             numericTools::LegendrePolyType::HALF_INTEGER, numericTools::HarmonicType::Outer,
                                             phi_tmp, dPhiDr_, dPhiDz_);
    phi_ += a[4] * phi_tmp;
    return phi_;
};

double DynamicCone::farFieldElectricFlux(double nr, double nz, double r, double z, const double (&b)[5])
{
    double phi_tmp = 0, dPhiDn_ = 0, dPhiDz_tmp = 0, dPhiDr_tmp = 0;
    numericTools::Geometry2D::zonalHarmonics(r, z, 0,
                                             numericTools::LegendrePolyType::HALF_INTEGER, numericTools::HarmonicType::Inner,
                                             phi_tmp, dPhiDr_tmp, dPhiDz_tmp);
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

Eigen::VectorXd DynamicCone::setVacuumBC(const bem2D::BoundaryElement &interface, const bem2D::BoundaryElement &vacuum) const
{
    int nCone = interface.node().size();
    int nPatch = vacuum.node().size();
    int nTotal = nCone + nPatch;
    Eigen::VectorXd vacuumBC;
    vacuumBC.setZero(nTotal);
    for (int k = 0; k < nTotal; k++)
    { // vacuum domain boundary  value problem (6.59)
        double r, z, dr, dz, ddr, ddz, nr, nz;

        if (k < nCone)
            vacuumBC(k) = 0.0;
        else
        {
            r = vacuum.node()[k - nCone].x;
            z = vacuum.node()[k - nCone].y;
            dr = vacuum.node()[k - nCone].dx;
            dz = vacuum.node()[k - nCone].dy;
            numericTools::Geometry2D::nrnz(dr, dz, nr, nz);
            vacuumBC(k) = farFieldElectricFlux(nr, nz, r, z, parameter.b);
        }
    }
    return vacuumBC;
};

Eigen::VectorXd DynamicCone::setLiquidBC(const bem2D::BoundaryElement &interface, const bem2D::BoundaryElement &liquid) const
{

    int nCone = interface.node().size();
    int nPatch = liquid.node().size();
    int nTotal = nCone + nPatch;
    Eigen::VectorXd liquidBC;

    liquidBC.setZero(nTotal);
    for (int k = 0; k < nTotal; k++)
    { // liquid domain boundary value problem (6.58)
        double r, z, dr, dz, ddr, ddz, nr, nz;
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
            liquidBC(k) = farFieldVelocityPotential(r, z, parameter.a);
        }
    }
    return liquidBC;
};

double DynamicCone::farFieldShape(double r, const double (&c)[5])
{
    return c[0] * r + c[1] / sqrt(r) + c[3] / pow(r, 3.5) + c[4] / pow(r, 5.0);
};

double DynamicCone::farFieldShapeDr(double r, const double (&c)[5])
{
    return c[0] - (5.0 * c[4]) / pow(r, 6.0) - (7 * c[3]) / (2. * pow(r, 4.5)) - c[1] / (2. * pow(r, 1.5));
};

double DynamicCone::farFieldShapeArc(double r0, double r1, const double (&c)[5])
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