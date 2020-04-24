#include <cmath>
#include <iostream>
#include "boundaryElement.hpp"
#include "numericTools/quadratureRules.hpp"
#include "numericTools/ellipticIntegral.hpp"
#include "numericTools/geometry2D.hpp"
#include "io/ioEigen.hpp"

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

void BoundaryElement::swapSDLR(const Eigen::MatrixXd &S, const Eigen::MatrixXd &D, int nSwap, Eigen::MatrixXd &L, Eigen::MatrixXd &R)
{
    int m = S.rows() - nSwap;
    L.topLeftCorner(nSwap, nSwap) = D.topLeftCorner(nSwap, nSwap);
    L.topRightCorner(nSwap, m) = -S.topRightCorner(nSwap, m);
    L.bottomLeftCorner(m, nSwap) = D.bottomLeftCorner(m, nSwap);
    L.bottomRightCorner(m, m) = -S.bottomRightCorner(m, m);

    R.topLeftCorner(nSwap, nSwap) = S.topLeftCorner(nSwap, nSwap);
    R.topRightCorner(nSwap, m) = -D.topRightCorner(nSwap, m);
    R.bottomLeftCorner(m, nSwap) = S.bottomLeftCorner(m, nSwap);
    R.bottomRightCorner(m, m) = -D.bottomRightCorner(m, m);
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

const BoundaryElement::Array6d BoundaryElement::regularGrad(double rp, double zp, int idElement, BoundaryElement::GradType type) const
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

    double a, b, m, K, E, f_single_K, f_single_E, f_double_K, f_double_E;

    for (int k = 0; k < static_cast<int>(nqd); k++)
    { // cumulative sum of gauss-legendre quadrature
        const double rk = r[k], zk = z[k], drk = dr[k], dzk = dz[k], Jk = J[k], wk = w[k];
        //double &f_single_K, double &f_double_K, double &f_single_E, double &f_double_E

        auxFunction_fKE(rp, zp, rk, zk, drk, dzk, Jk, type, f_single_K, f_double_K, f_single_E, f_double_E);
        auxFunction_abm(rp, zp, rk, zk, a, b, m);
        numericTools::EllipticIntegral::ellipticKE(m, K, E);

        for (int j = 0; j <= o; j++)
        { // iterate through all Lagrange basis functions of order o
            const double &N = basis[j][k];
            output[j] += (f_single_K * K + f_single_E * E) * Jk * N * wk;    // equation 7.69
            output[o + 1 + j] += (f_double_K * K + f_double_E * E) * N * wk; // equation 7.74
        }
    }
    return output;
}

const BoundaryElement::Array6d BoundaryElement::singular(double tau, int idElement) const
{
    Array6d output;
    output.setZero();
    double rp, zp;
    setSourcePoint(tau, idElement, rp, zp);
    const int nqdRegular = 16;  // settings.qdOrder() * 2;
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
        const double rk = r[k], zk = z[k], drk = dr[k], dzk = dz[k], Jk = J[k], wk = w[k]; //, abk = ab[k];
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
        const double rk = r[k], zk = z[k], drk = dr[k], dzk = dz[k], Jk = J[k], wk = w[k]; //, abk = ab[k];
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

// const BoundaryElement::Array6d BoundaryElement::singular(double tau, int idElement) const
// {
//     //const auto &e = element()[idElement];
//     //if (e.elementOrder() == 1)
//     //  return singularFirstOrder(tau, idElement);
//     //else
//     return singularHigherOrder(tau, idElement);
// }

const BoundaryElement::Array6d BoundaryElement::axis(double zp, int idElement) const
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

void BoundaryElement::assembleMatrix(const BoundaryElement &bem0,
                                     const BoundaryElement &bem1, Eigen::MatrixXd &S, Eigen::MatrixXd &D)
{
    const int nNode = bem0.node().size();
    const int shift0 = bem0.indexShift();
    //  const Eigen::VectorXd &r = bem0.node().r.col(0);
    //    const Eigen::VectorXd &z = bem0.node().z.col(0);
    const int shift1 = bem1.indexShift();
    const int o = bem1.elementOrder();
    const int nElem = bem1.element().size();
    const auto check = checkBoundaryRelation(bem0, bem1);

#pragma omp parallel for
    for (int i = 0; i < nNode; i++)
    { // loop through all field points of bem0
        double ri = bem0.node()[i].x, zi = bem0.node()[i].y;
        if (numericTools::Geometry2D::isOnZAxis(ri)) // field point at symmetry axis r = 0
        {
            for (int j = 0; j < nElem; j++)
            {
                const Array6d axisIntegral = bem1.axis(zi, j);
                for (int k = 0; k <= o; k++)
                {
                    S(i + shift0, shift1 + o * j + k) += axisIntegral[k];
                    D(i + shift0, shift1 + o * j + k) += axisIntegral[o + 1 + k];
                }
            }
        }
        else
        {
            switch (check)
            {
            case BoundaryRelationType::Identical:
            {
                for (int j = 0; j < i / o - 1 + (i % o); j++)
                { // regular integral
                    const Array6d regularIntegral = bem1.regular(ri, zi, j);
                    for (int k = 0; k <= o; k++)
                    {
                        S(i + shift0, o * j + k + shift1) += regularIntegral[k];
                        D(i + shift0, o * j + k + shift1) += regularIntegral[o + 1 + k];
                    }
                }
                if (i % o == 0)
                {                      // two elements that contains singularity: m - 1 and m
                    int j = i / o - 1; // element m - 1
                    if (j >= 0 && j < nElem)
                    {
                        const Array6d singularIntegral = bem1.singular(1, j); // singularit at local tau = 1
                        for (int k = 0; k <= o; k++)
                        {
                            S(i + shift0, o * j + k + shift1) += singularIntegral[k];
                            D(i + shift0, o * j + k + shift1) += singularIntegral[o + 1 + k];
                        }
                    }

                    j = i / o; // element m
                    if (j >= 0 && j < nElem)
                    {
                        const Array6d singularIntegral = bem1.singular(0, j); // singularit at local tau = 0
                        for (int k = 0; k <= o; k++)
                        {
                            S(i + shift0, o * j + k + shift1) += singularIntegral[k];
                            D(i + shift0, o * j + k + shift1) += singularIntegral[o + 1 + k];
                        }
                    }
                }
                else
                {                                                                                // only one spline that contains singularity: m
                    int j = i / o;                                                               // m
                    const Array6d singularIntegral = bem1.singular(bem1.element()[j].t()[1], j); // singularit at t where half arclength happens
                    for (int k = 0; k <= o; k++)
                    {
                        S(i + shift0, o * j + k + shift1) += singularIntegral[k];
                        D(i + shift0, o * j + k + shift1) += singularIntegral[o + 1 + k];
                    }
                }

                for (int j = i / o + 1; j < nElem; j++)
                { // regular integral
                    const Array6d regularIntegral = bem1.regular(ri, zi, j);
                    for (int k = 0; k <= o; k++)
                    {
                        S(i + shift0, o * j + k + shift1) += regularIntegral[k];
                        D(i + shift0, o * j + k + shift1) += regularIntegral[o + 1 + k];
                    }
                }
                break;
            }
            case BoundaryRelationType::JoinedBefore:
            {
                for (int j = 0; j < nElem; j++)
                {
                    if (i == 0 && j == nElem - 1)
                    { // singular: 1st field point of bem0 meets last source element of bem1
                        const Array6d singularIntegral = bem1.singular(1.0, j);
                        for (int k = 0; k <= o; k++)
                        {
                            S(i + shift0, o * j + k + shift1) += singularIntegral[k];
                            D(i + shift0, o * j + k + shift1) += singularIntegral[o + 1 + k];
                        }
                    }
                    else
                    { // most field-source pair are disjoint hence regular integral
                        const Array6d regularIntegral = bem1.regular(ri, zi, j);
                        for (int k = 0; k <= o; k++)
                        {
                            S(i + shift0, o * j + k + shift1) += regularIntegral[k];
                            D(i + shift0, o * j + k + shift1) += regularIntegral[o + 1 + k];
                        }
                    }
                }
                break;
            }
            case BoundaryRelationType::JoinedAfter:
            {
                for (int j = 0; j < nElem; j++)
                {
                    if (i == nNode - 1 && j == 0)
                    { // singular: last field point of bem0 meets 1st source element of bem1
                        const Array6d singularIntegral = bem1.singular(0, j);
                        for (int k = 0; k <= o; k++)
                        {
                            S(i + shift0, o * j + k + shift1) += singularIntegral[k];
                            D(i + shift0, o * j + k + shift1) += singularIntegral[o + 1 + k];
                        }
                    }
                    else
                    { // most field-source pair are disjoint hence regular integral
                        const Array6d regularIntegral = bem1.regular(ri, zi, j);
                        for (int k = 0; k <= o; k++)
                        {
                            S(i + shift0, o * j + k + shift1) += regularIntegral[k];
                            D(i + shift0, o * j + k + shift1) += regularIntegral[o + 1 + k];
                        }
                    }
                }
                break;
            }
            case BoundaryRelationType::Disjoint:
            {
                for (int j = 0; j < nElem; j++)
                { // only need regular integral
                    const Array6d regularIntegral = bem1.regular(ri, zi, j);
                    for (int k = 0; k <= o; k++)
                    {
                        S(i + shift0, o * j + k + shift1) += regularIntegral[k];
                        D(i + shift0, o * j + k + shift1) += regularIntegral[o + 1 + k];
                    }
                }
                break;
            }
            default:
                break;
            }
        }
    }
};

Eigen::Vector3d BoundaryElement::interiorField(double rp, double zp, const BoundaryElement &bem, const Eigen::VectorXd &q, const Eigen::VectorXd &p)
{
    const int shift1 = bem.indexShift();
    const int o = bem.elementOrder();
    const int nElem = bem.nElement();

    Eigen::VectorXd S, D, Sr, Dr, Sz, Dz;
    S.setZero(q.size());
    D.setZero(p.size());
    Sr.setZero(q.size());
    Dr.setZero(p.size());
    Sz.setZero(q.size());
    Dz.setZero(p.size());

    if (numericTools::Geometry2D::isOnZAxis(rp))
    {
        for (int j = 0; j < nElem; j++)
        {
            const Array6d axisIntegral = bem.axis(zp, j);
            for (int k = 0; k <= o; k++)
            {
                S(shift1 + o * j + k) += axisIntegral[k];
                D(shift1 + o * j + k) += axisIntegral[o + 1 + k];
            }
        }
    }
    else
    {
        for (int j = 0; j < nElem; j++)
        {
            const Array6d regularIntegral = bem.regular(rp, zp, j);
            for (int k = 0; k <= o; k++)
            {
                S(o * j + k + shift1) += regularIntegral[k];
                D(o * j + k + shift1) += regularIntegral[o + 1 + k];
            }
        }
    }

    for (int j = 0; j < nElem; j++)
    {
        const Array6d regularIntegralDr = bem.regularGrad(rp, zp, j, GradType::Dr);
        const Array6d regularIntegralDz = bem.regularGrad(rp, zp, j, GradType::Dz);
        for (int k = 0; k <= o; k++)
        {
            Sr(o * j + k + shift1) += regularIntegralDr[k];
            Dr(o * j + k + shift1) += regularIntegralDr[o + 1 + k];
            Sz(o * j + k + shift1) += regularIntegralDz[k];
            Dz(o * j + k + shift1) += regularIntegralDz[o + 1 + k];
        }
    }

    Eigen::Vector3d tmp;
    tmp.setZero();
    tmp(0) = S.dot(q) - D.dot(p);
    tmp(1) = Sr.dot(q) - Dr.dot(p);
    tmp(2) = Sz.dot(q) - Dz.dot(p);
    return tmp;
};

// ================ rule of five  ==================//

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
    const auto &o = e.elementOrder(); //std::cout << rp << "\t" << zp << "\n";
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

void BoundaryElement::auxFunction_fKE(double rp, double zp, double r, double z, double dr, double dz, double J, BoundaryElement::GradType type,
                                      double &f_single_K, double &f_double_K, double &f_single_E, double &f_double_E)
{
    switch (type)
    {
    case BoundaryElement::GradType::Dr:
    {
        f_single_K = -r / (2. * M_PI * rp * sqrt(pow(r + rp, 2) + pow(z - zp, 2)));
        f_single_E = (pow(r, 3) - r * (rp + z - zp) * (rp - z + zp)) / (2. * M_PI * rp * (pow(r - rp, 2) + pow(z - zp, 2)) * sqrt(pow(r + rp, 2) + pow(z - zp, 2)));
        f_double_K = (-(dz * (pow(r, 4) + pow(r, 2) * (-2 * pow(rp, 2) + pow(z - zp, 2)) + pow(rp, 2) * (pow(rp, 2) + pow(z - zp, 2)))) + dr * r * (pow(r, 2) - pow(rp, 2) + pow(z - zp, 2)) * (z - zp)) / (2. * M_PI * rp * (pow(r - rp, 2) + pow(z - zp, 2)) * pow(pow(r + rp, 2) + pow(z - zp, 2), 1.5));
        f_double_E = (dz * (pow(r, 6) - pow(r, 4) * (pow(rp, 2) - 2 * pow(z - zp, 2)) + pow(rp, 2) * pow(pow(rp, 2) + pow(z - zp, 2), 2) + pow(r, 2) * (-pow(rp, 4) - 12 * pow(rp, 2) * pow(z - zp, 2) + pow(z - zp, 4))) - dr * r * (pow(r, 4) - 7 * pow(rp, 4) + 2 * pow(r, 2) * (3 * pow(rp, 2) + pow(z - zp, 2)) - 6 * pow(rp, 2) * pow(z - zp, 2) + pow(z - zp, 4)) * (z - zp)) / (2. * M_PI * rp * pow(pow(r - rp, 2) + pow(z - zp, 2), 2) * pow(pow(r + rp, 2) + pow(z - zp, 2), 1.5));
        break;
    }
    case BoundaryElement::GradType::Dz:
    {
        f_single_K = 0;
        f_single_E = (r * (z - zp)) / (M_PI * (pow(r - rp, 2) + pow(z - zp, 2)) * sqrt(pow(r + rp, 2) + pow(z - zp, 2)));
        f_double_K = ((dz * (-pow(r, 2) + pow(rp, 2) + pow(z - zp, 2)) + 2 * dr * r * (z - zp)) * (z - zp)) / (2. * M_PI * (pow(r - rp, 2) + pow(z - zp, 2)) * pow(pow(r + rp, 2) + pow(z - zp, 2), 1.5));
        f_double_E = (2 * dr * r * (pow(r, 4) - 2 * pow(r, 2) * (pow(rp, 2) + pow(z - zp, 2)) + (pow(rp, 2) - 3 * pow(z - zp, 2)) * (pow(rp, 2) + pow(z - zp, 2))) - dz * (z - zp) * (-7 * pow(r, 4) + pow(pow(rp, 2) + pow(z - zp, 2), 2) + 6 * pow(r, 2) * (rp + z - zp) * (rp - z + zp))) / (2. * M_PI * pow(pow(r - rp, 2) + pow(z - zp, 2), 2) * pow(pow(r + rp, 2) + pow(z - zp, 2), 1.5));
        break;
    }
    default:
        break;
    };
};

const BoundaryElement::BoundaryRelationType BoundaryElement::checkBoundaryRelation(const BoundaryElement &bem0, const BoundaryElement &bem1)
{
    int s0 = bem0.indexShift();
    int s1 = bem1.indexShift();
    int ds = s1 - s0;

    if (ds == 0) // same Bem object
        return BoundaryRelationType::Identical;
    if (ds > 0 && std::abs(ds) == bem0.node().size()) // bem1 is adjacent to and "after" bem0
        return BoundaryRelationType::JoinedAfter;
    if (ds < 0 && std::abs(ds) == bem1.node().size()) // bem1 is adjacent to and "before" bem0
        return BoundaryRelationType::JoinedBefore;

    return BoundaryRelationType::Disjoint; // two Bem objects are not adjacent
};

void BoundaryElement::scan(const BoundaryElement &bem0, const BoundaryElement &bem1,
                           const Eigen::VectorXd &q, const Eigen::VectorXd &p, const std::string &inputFile, const std::string &outputFile)
{
    Eigen::MatrixXd inputMat = io::IOEigen::readMatrix(inputFile.c_str(), 1e6);
    Eigen::MatrixXd outputMat;
    outputMat.setZero(inputMat.rows() + p.size(), 5);
    // evaluate phi, dphidr, dphidz at these points
#pragma omp parallel for
    for (Eigen::Index row = 0; row < inputMat.rows(); row++)
    {
        double rp = inputMat(row, 0);
        double zp = inputMat(row, 1);
        Eigen::Vector3d field;
        field.setZero();
        field = BoundaryElement::interiorField(rp, zp, bem0, q, p);
        field += BoundaryElement::interiorField(rp, zp, bem1, q, p);
        outputMat(row, 0) = rp;
        outputMat(row, 1) = zp;
        outputMat(row, 2) = field(0);
        outputMat(row, 3) = field(1);
        outputMat(row, 4) = field(2);
    }

    // export field along domain boundary
    int n0 = bem0.node().size(), n1 = bem1.node().size();
    int nTotal = n0 + n1;

    for (int k = 0; k < nTotal; k++)
    {
        double rp, zp, field;
        if (k < n0)
        { // boundary 0
            rp = bem0.node()[k].x;
            zp = bem0.node()[k].y;
        }
        else
        {                               // boundary 1
            rp = bem1.node()[k - n0].x; //.r(k - n0, 0);
            zp = bem1.node()[k - n0].y; //.z(k - n0, 0);
        }
        field = p(k);
        int row = inputMat.rows() + k;
        outputMat(row, 0) = rp;
        outputMat(row, 1) = zp;
        outputMat(row, 2) = field;
    }

    io::IOEigen::write(outputFile, outputMat);
};

} // namespace bem2D