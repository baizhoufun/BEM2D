#include <cmath>
#include <iostream>
#include "boundaryElement.hpp"
#include "numericTools/quadratureRules.hpp"
#include "numericTools/ellipticIntegral.hpp"

namespace bem2D
{

void BoundaryElement::initializeElement()
{
    nElement(sp().node().rows() - 1);
    _element.resize(nElement());
    for (std::size_t i = 0; i < element().size(); i++)
        _element[i].initialize(sp(), i, elementOrder(), quadratureOrder());
};

void BoundaryElement::initializeNode()
{
    const int nE = nElement();
    const int o = elementOrder();
    const int nN = nE * o + 1;

    _node.resize(nN);

    for (int idElement = 0; idElement < nE; idElement++)
    {
        const std::vector<double> &t = element()[idElement].t(); // grab intrinsic coordinates tk of an element
        for (int k = 0; k < o; k++)
        { // evaluate r(t), r'(t) and r''(t) at tk
            int idNode = idElement * o + k;
            auto xDiff = sp().d(sp().x(), idElement, t[k]);
            auto yDiff = sp().d(sp().y(), idElement, t[k]);
            auto &tmp = _node[idNode];
            tmp.x = xDiff(0);
            tmp.dx = xDiff(1);
            tmp.ddx = xDiff(2);
            tmp.y = yDiff(0);
            tmp.dy = yDiff(1);
            tmp.ddy = yDiff(2);
        }
    }
    int idNode = _node.size() - 1;
    int idElement = nE - 1;

    auto xDiff = sp().d(sp().x(), idElement + 1, 0);
    auto yDiff = sp().d(sp().y(), idElement + 1, 0.0);
    auto &tmp = _node[idNode];
    tmp.x = xDiff(0);
    tmp.dx = xDiff(1);
    tmp.ddx = xDiff(2);
    tmp.y = yDiff(0);
    tmp.dy = yDiff(1);
    tmp.ddy = yDiff(2);
};

void BoundaryElement::bindNodeToElement()
{
    const int nE = nElement();
    const int o = elementOrder();

    for (int idElement = 0; idElement < nE; idElement++)
    {
        for (int k = 0; k <= o; k++)
        {
            int idNode = idElement * o + k;
            _element[idElement].setNodeIndex()[k] = idNode;
        }
    }
}

void BoundaryElement::initialize()
{
    initializeElement();
    initializeNode();
    bindNodeToElement();
}

// ===== single/double integration =====

const BoundaryElement::Array6d BoundaryElement::regular(double rp, double zp, int idElement) const
{
    Array6d output;
    output.setZero();
    const auto &e = element()[idElement];
    const auto &nqd = e.x().size();
    const auto &o = e.elementOrder();
    const auto &basis = e.basis();
    const auto &r = e.x();
    const auto &dr = e.dx();
    const auto &z = e.y();
    const auto &dz = e.dy();
    const auto &J = e.J();
    const double *w;
    w = numericTools::QuadratureRules::weight(nqd, numericTools::QuadratureType::GAUSS_LEGENDRE);

    double a, b, m, K, E, f_single_K, f_double_K, f_double_E;

    for (int k = 0; k < static_cast<int>(nqd); k++)
    { // cumulative sum of gauss-legendre quadrature
        const double rk = r[k], zk = z[k], drk = dr[k], dzk = dz[k], Jk = J[k], wk = w[k];

        auxFunction_abm(rp, zp, rk, zk, a, b, m);
        numericTools::EllipticIntegral::ellipticKE(m, K, E);
        auxFunction_fKE(rp, zp, rk, zk, drk, dzk, Jk, a, b, f_single_K, f_double_K, f_double_E);

        for (int j = 0; j <= o; j++)
        { // iterate through all Lagrange basis functions of order o
            const double &N = basis[j][k];
            output[j] += f_single_K * K * N * wk;                            // equation 7.69
            output[o + 1 + j] += (f_double_K * K + f_double_E * E) * N * wk; // equation 7.74
        }
    }
    return output;
}

const BoundaryElement::Array6d BoundaryElement::singularFirstOrder(double tau, int idElement) const
{
    Array6d output;
    output.setZero();

    double rp, zp;
    setSourcePoint(tau, idElement, rp, zp);
    const int nqdRegular = 20;  // settings.qdOrder() * 2;
    const int nqdSingular = 20; //settings.qdOrder() * 2;
    double a, b, m, K, E, PK, QK, PE, QE, RK, RE, f_single_K, f_double_K, f_double_E;
    Element e;
    const double *ab;
    const double *w;

    ab = numericTools::QuadratureRules::abascissa(nqdRegular, numericTools::QuadratureType::GAUSS_LEGENDRE);
    w = numericTools::QuadratureRules::weight(nqdRegular, numericTools::QuadratureType::GAUSS_LEGENDRE);
    e.initialize(sp(), idElement, 1, nqdRegular, ab);
    auto o = e.elementOrder();
    auto basis = e.basis();
    auto r = e.x();
    auto dr = e.dx();
    auto z = e.y();
    auto dz = e.dy();
    auto J = e.J();

    for (int k = 0; k < static_cast<int>(nqdRegular); k++)
    {
        const double rk = r[k], zk = z[k], drk = dr[k], dzk = dz[k], Jk = J[k], wk = w[k], abk = ab[k];
        auxFunction_abm(rp, zp, rk, zk, a, b, m);
        numericTools::EllipticIntegral::ellipticKEPQ(m, K, E, PK, QK, PE, QE);
        RK = auxFunction_RKE(PK, QK, m, abk, tau);
        RE = auxFunction_RKE(PE, QE, m, abk, tau);
        auxFunction_fKE(rp, zp, rk, zk, drk, dzk, Jk, a, b, f_single_K, f_double_K, f_double_E);
        for (int j = 0; j <= o; j++)
        { // iterate through all Lagrange basis functions of order o
            const double &N = basis[j][k];
            output[j] += f_single_K * RK * N * wk;                             // equation 7.69
            output[o + 1 + j] += (f_double_K * RK + f_double_E * RE) * N * wk; // equation 7.74
        }
    }

    const double *abSingular = numericTools::QuadratureRules::abascissa(nqdSingular, numericTools::QuadratureType::GAUSS_LEGENDRE_LOG);
    w = numericTools::QuadratureRules::weight(nqdSingular, numericTools::QuadratureType::GAUSS_LEGENDRE_LOG);

    double *abb = new double[nqdSingular];
    for (int l = 0; l < nqdSingular; l++)
    { // if singularity tau = 0 then use abscissa, if tau = 1, use 1 - ab
        double abl = numericTools::QuadratureRules::abascissa(nqdSingular, numericTools::QuadratureType::GAUSS_LEGENDRE_LOG)[l];
        abb[l] = (1. - tau) * abl + tau * (1. - abl); // equation (7.36)
    }
    ab = abb;
    e.initialize(sp(), idElement, 1, nqdSingular, ab);

    o = e.elementOrder();
    basis = e.basis();
    r = e.x();
    dr = e.dx();
    z = e.y();
    dz = e.dy();
    J = e.J();

    for (int k = 0; k < static_cast<int>(nqdSingular); k++)
    {
        const double rk = r[k], zk = z[k], drk = dr[k], dzk = dz[k], Jk = J[k], wk = w[k], abk = ab[k];
        auxFunction_abm(rp, zp, rk, zk, a, b, m);
        numericTools::EllipticIntegral::ellipticKEPQ(m, K, E, PK, QK, PE, QE);

        auxFunction_fKE(rp, zp, rk, zk, drk, dzk, Jk, a, b, f_single_K, f_double_K, f_double_E);
        for (int j = 0; j <= o; j++)
        { // iterate through all Lagrange basis functions of order o
            const double &N = basis[j][k];
            output[j] += 2. * f_single_K * QK * N * wk;                             // equation 7.69
            output[o + 1 + j] += 2. * (f_double_K * QK + f_double_E * QE) * N * wk; // equation 7.74
        }
    }
    delete[] abb;
    return output;
};

const BoundaryElement::Array6d BoundaryElement::singularHigherOrder(double tau, int idElement) const
{
    Array6d output;
    output.setZero();
    double rp, zp;
    setSourcePoint(tau, idElement, rp, zp);
    const int nqdRegular = 12;  // settings.qdOrder() * 2;
    const int nqdSingular = 16; //settings.qdOrder() * 2;
    double a, b, m, K, E, PK, QK, PE, QE, RK, RE, f_single_K, f_double_K, f_double_E;
    Element e;
    const double *ab;
    const double *w;
    double *abb;

    //=================================================

    w = numericTools::QuadratureRules::weight(nqdRegular, numericTools::QuadratureType::GAUSS_LEGENDRE);
    abb = new double[nqdRegular];
    for (int l = 0; l < nqdRegular; l++)
    {
        double abl = numericTools::QuadratureRules::abascissa(nqdRegular, numericTools::QuadratureType::GAUSS_LEGENDRE)[l];
        abb[l] = tau * abl; // equation (7.36)
    }
    ab = abb;

    e.initialize(sp(), idElement, element()[idElement].elementOrder(), nqdRegular, ab);
    auto o = e.elementOrder();
    auto basis = e.basis();
    auto r = e.x();
    auto dr = e.dx();
    auto z = e.y();
    auto dz = e.dy();
    auto J = e.J();

    for (int k = 0; k < static_cast<int>(nqdRegular); k++)
    {

        const double rk = r[k], zk = z[k], drk = dr[k], dzk = dz[k], Jk = J[k], wk = w[k], abk = ab[k];
        auxFunction_abm(rp, zp, rk, zk, a, b, m);
        numericTools::EllipticIntegral::ellipticKEPQ(m, K, E, PK, QK, PE, QE);

        RK = auxFunction_RKE(PK, QK, m, abk, tau);
        RE = auxFunction_RKE(PE, QE, m, abk, tau);
        auxFunction_fKE(rp, zp, rk, zk, drk, dzk, Jk, a, b, f_single_K, f_double_K, f_double_E);

        for (int j = 0; j <= o; j++)
        { // iterate through all Lagrange basis functions of order o
            const double &N = basis[j][k];
            output[j] += f_single_K * RK * N * wk * tau - 2. * f_single_K * QK * N * wk * auxFunction_xLogX(tau); // equation 7.69
            output[o + 1 + j] += (f_double_K * RK + f_double_E * RE) * N * wk * tau - 2. * (f_double_K * QK + f_double_E * QE) * N * wk * auxFunction_xLogX(tau);
        }
    }
    delete[] abb;

    //=================================================

    abb = new double[nqdRegular];
    for (int l = 0; l < nqdRegular; l++)
    {
        double abl = numericTools::QuadratureRules::abascissa(nqdRegular, numericTools::QuadratureType::GAUSS_LEGENDRE)[l];
        abb[l] = tau + (1 - tau) * abl; // equation (7.36)
    }
    ab = abb;

    e.initialize(sp(), idElement, element()[idElement].elementOrder(), nqdRegular, ab);
    o = e.elementOrder();
    basis = e.basis();
    r = e.x();
    dr = e.dx();
    z = e.y();
    dz = e.dy();
    J = e.J();
    for (int k = 0; k < static_cast<int>(nqdRegular); k++)
    {
        const double rk = r[k], zk = z[k], drk = dr[k], dzk = dz[k], Jk = J[k], wk = w[k], abk = ab[k];
        auxFunction_abm(rp, zp, rk, zk, a, b, m);
        numericTools::EllipticIntegral::ellipticKEPQ(m, K, E, PK, QK, PE, QE);
        RK = auxFunction_RKE(PK, QK, m, abk, tau);
        RE = auxFunction_RKE(PE, QE, m, abk, tau);
        auxFunction_fKE(rp, zp, rk, zk, drk, dzk, Jk, a, b, f_single_K, f_double_K, f_double_E);
        for (int j = 0; j <= o; j++)
        { // iterate through all Lagrange basis functions of order o
            const double &N = basis[j][k];
            output[j] += f_single_K * RK * N * wk * (1. - tau) - 2. * f_single_K * QK * N * wk * auxFunction_xLogX(1. - tau); // equation 7.69
            output[o + 1 + j] += (f_double_K * RK + f_double_E * RE) * N * wk * (1. - tau) - 2. * (f_double_K * QK + f_double_E * QE) * N * wk * auxFunction_xLogX(1. - tau);
        }
    }

    delete[] abb;
    //=================================================

    w = numericTools::QuadratureRules::weight(nqdSingular, numericTools::QuadratureType::GAUSS_LEGENDRE_LOG);
    abb = new double[nqdSingular];
    for (int l = 0; l < nqdSingular; l++)
    {
        double abl = numericTools::QuadratureRules::abascissa(nqdSingular, numericTools::QuadratureType::GAUSS_LEGENDRE_LOG)[l];
        abb[l] = tau * (1 - abl); // equation (7.36)
    }
    ab = abb;

    e.initialize(sp(), idElement, element()[idElement].elementOrder(), nqdSingular, ab);
    o = e.elementOrder();
    basis = e.basis();
    r = e.x();
    dr = e.dx();
    z = e.y();
    dz = e.dy();
    J = e.J();
    for (int k = 0; k < static_cast<int>(nqdSingular); k++)
    {
        const double rk = r[k], zk = z[k], drk = dr[k], dzk = dz[k], Jk = J[k], wk = w[k], abk = ab[k];
        auxFunction_abm(rp, zp, rk, zk, a, b, m);
        numericTools::EllipticIntegral::ellipticKEPQ(m, K, E, PK, QK, PE, QE);
        auxFunction_fKE(rp, zp, rk, zk, drk, dzk, Jk, a, b, f_single_K, f_double_K, f_double_E);
        for (int j = 0; j <= o; j++)
        { // iterate through all Lagrange basis functions of order o
            const double &N = basis[j][k];
            output[j] += 2. * f_single_K * QK * N * wk * tau; // equation 7.69
            output[o + 1 + j] += 2. * (f_double_K * QK + f_double_E * QE) * N * wk * tau;
        }
    }

    delete[] abb;

    //=================================================

    //w = numericTools::QuadratureRules::weight(nqdSingular, numericTools::QuadratureType::GAUSS_LEGENDRE_LOG);
    abb = new double[nqdSingular];
    for (int l = 0; l < nqdSingular; l++)
    {
        double abl = numericTools::QuadratureRules::abascissa(nqdSingular, numericTools::QuadratureType::GAUSS_LEGENDRE_LOG)[l];
        abb[l] = tau + (1 - tau) * abl; // equation (7.36)
    }
    ab = abb;

    e.initialize(sp(), idElement, element()[idElement].elementOrder(), nqdSingular, ab);
    o = e.elementOrder();
    basis = e.basis();
    r = e.x();
    dr = e.dx();
    z = e.y();
    dz = e.dy();
    J = e.J();
    for (int k = 0; k < static_cast<int>(nqdSingular); k++)
    {
        const double rk = r[k], zk = z[k], drk = dr[k], dzk = dz[k], Jk = J[k], wk = w[k], abk = ab[k];
        auxFunction_abm(rp, zp, rk, zk, a, b, m);
        numericTools::EllipticIntegral::ellipticKEPQ(m, K, E, PK, QK, PE, QE);
        auxFunction_fKE(rp, zp, rk, zk, drk, dzk, Jk, a, b, f_single_K, f_double_K, f_double_E);
        for (int j = 0; j <= o; j++)
        { // iterate through all Lagrange basis functions of order o
            const double &N = basis[j][k];
            output[j] += 2. * f_single_K * QK * N * wk * (1. - tau); // equation 7.69
            output[o + 1 + j] += 2. * (f_double_K * QK + f_double_E * QE) * N * wk * (1. - tau);
        }
    }
    delete[] abb;

    return output;
};

const BoundaryElement::Array6d BoundaryElement::singular(double tau, int idElement) const
{
    const auto &e = element()[idElement];
    if (e.elementOrder() == 1)
        return singularFirstOrder(tau, idElement);
    else
        return singularHigherOrder(tau, idElement);
}

const BoundaryElement::Array6d BoundaryElement::axis(double zp, int idElement) const
{
    Array6d output;
    const auto &e = element()[idElement];
    const auto &nqd = e.x().size();
    const auto &o = e.elementOrder();
    const auto &basis = e.basis();
    const auto &r = e.x();
    const auto &dr = e.dx();
    const auto &z = e.y();
    const auto &dz = e.dy();
    const auto &J = e.J();
    const double *w;
    double f_single_axis, f_double_axis;

    for (int k = 0; k < static_cast<int>(nqd); k++)
    {
        w = numericTools::QuadratureRules::weight(nqd, numericTools::QuadratureType::GAUSS_LEGENDRE);
        const double rk = r[k], zk = z[k], drk = dr[k], dzk = dz[k], Jk = J[k], wk = w[k];

        f_single_axis = Jk * rk / 2. / sqrt(rk * rk + (zk - zp) * (zk - zp));
        f_double_axis = rk / 2. * (dzk * rk + drk * (zp - zk)) / pow(rk * rk + (zk - zp) * (zk - zp), 1.5);

        for (int j = 0; j <= o; j++)
        {
            const double &N = basis[j][k];
            output[j] += f_single_axis * N * wk;
            output[o + 1 + j] += f_double_axis * N * wk;
        }
    }
    return output;
};

BoundaryElement::BoundaryElement(){};

BoundaryElement::~BoundaryElement(){};

// ================ getter and setter ==================//
const int &BoundaryElement::nElement() const
{
    return _nElement;
};
void BoundaryElement::nElement(int i)
{
    _nElement = i;
};

const int &BoundaryElement::elementOrder() const
{
    return _elementOrder;
};
void BoundaryElement::elementOrder(int i)
{
    _elementOrder = i;
};

const int &BoundaryElement::indexShift() const
{
    return _indexShift;
};
void BoundaryElement::indexShift(int i)
{
    _indexShift = i;
};

const int &BoundaryElement::quadratureOrder() const
{
    return _quadratureOrder;
};
void BoundaryElement::quadratureOrder(int i)
{
    _quadratureOrder = i;
};

const spline::Quintic &BoundaryElement::sp() const
{
    return _sp;
};
spline::Quintic &BoundaryElement::sp()
{
    return _sp;
};

const std::vector<Element> &BoundaryElement::element() const
{
    return _element;
};

const std::vector<BoundaryElement::Node> &BoundaryElement::node() const
{
    return _node;
};

// ========== aux function
void BoundaryElement::auxFunction_abm(double rp, double zp, double r, double z, double &a, double &b, double &m)
{
    a = r * r + rp * rp + (zp - z) * (zp - z);
    b = 2. * rp * r;
    m = 2. * b / (a + b);
};

double BoundaryElement::auxFunction_RKE(double P, double Q, double m, double t, double tp)
{

    if (std::abs(1.0 - m) < 1e-13)
        return 0;
    else
        return P - Q * log((1. - m) / (t * t - 2 * tp * t + tp * tp));
};

void BoundaryElement::auxFunction_fKE(double rp, double zp, double r, double z, double dr, double dz, double J, double a, double b,
                                      double &f_single_K, double &f_double_K, double &f_double_E)
{
    f_single_K = r * J / M_PI / sqrt(a + b);   // equation (7.68)
    f_double_K = dz / 2. / M_PI / sqrt(a + b); // equation (7.73)

    if (std::abs(a - b) < 1e-13)
        f_double_E = 0;
    else
        f_double_E = ((dr * (zp - z) - dz * (rp - r)) * r / (a - b) - dz / 2.) / M_PI / sqrt(a + b);
};

void BoundaryElement::setSourcePoint(double tau, int idElement, double &rp, double &zp) const
{
    const auto &e = element()[idElement];
    const auto &o = e.elementOrder();
    double epsSnap = 1e-13;
    if (std::abs(tau) < epsSnap)
    {
        rp = node()[idElement * o].x;
        zp = node()[idElement * o].y;
    }
    else if (std::abs(1.0 - tau) < epsSnap)
    {
        rp = node()[(idElement + 1) * o].x;
        zp = node()[(idElement + 1) * o].y;
    }
    else
    {
        rp = sp().d(sp().x(), idElement, tau)(0);
        zp = sp().d(sp().y(), idElement, tau)(0);
    }

    //std::cout << rp << "\t" << zp << "\n";
};
double BoundaryElement::auxFunction_xLogX(double x)
{
    if (std::abs(x) < 1e-14)
    {
        return 0;
    }
    else
        return x * std::log(x);
};

} // namespace bem2D