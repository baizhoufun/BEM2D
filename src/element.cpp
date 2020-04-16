#include "element.hpp"
#include "numericTools/quadratureRules.hpp"
#include "numericTools/lagrangeBasis.hpp"

using namespace numericTools;

namespace bem2D
{

Element::Element()
{
}

Element::~Element()
{
}
Element::Element(const Element &that)
{
    _arc = that.arc();
    _elementOrder = that.elementOrder();
    _quadratureOrder = that.quadratureOrder();
    _t = that.t();
    _x = that.x();
    _y = that.y();
    _dx = that.dx();
    _dy = that.dy();
    _J = that.J();
    _xi = that.xi();
    _basis = that.basis();
};

void Element::initialize(const spline::Quintic &sp, int iSplineSegment, int nElementOrder, int nQuadratureOrder, double const *abacissa)
{
    _arc = sp.localArc(iSplineSegment); // compute total arc of local spline sp
    _elementOrder = nElementOrder;      // element order
    _quadratureOrder = nQuadratureOrder;
    const int &i = iSplineSegment;
    const int &o = elementOrder();
    const int &nqd = quadratureOrder();
    const double *ab;
    if (abacissa == nullptr)
        ab = QuadratureRules::abascissa(nqd, QuadratureType::GAUSS_LEGENDRE);
    else
        ab = abacissa;
    // allocation
    _t.resize(o + 1); // nodes of element; linear t ={ 0, 1}
    _x.resize(nqd);   // various members evaluated at nqd-point quadrature abscissa
    _dx.resize(nqd);
    _y.resize(nqd);
    _dy.resize(nqd);
    _J.resize(nqd);
    _xi.resize(nqd);
    _basis.resize(o + 1); // o + 1 basis fuctions for o-th Lagrange basis
    for (auto &member : _basis)
        member.resize(nqd);
    // compute
    _t[0] = 0.0;
    _t[1] = sp.arc2t(i, 0.5 * arc()); // quadratic  t = {0, where xi(t) = 0.5,  1}
    _t[o] = 1.0;
    for (int k = 0; k < nqd; k++)
    { // evaluate data members at nqd-point quadrature abscissa
        const Eigen::Vector3d xDiff = sp.d(sp.x(), i, ab[k]);
        const Eigen::Vector3d yDiff = sp.d(sp.y(), i, ab[k]);

        _x[k] = xDiff(0);
        _y[k] = yDiff(0);
        _dx[k] = xDiff(1);
        _dy[k] = yDiff(1);
        _J[k] = sqrt(xDiff(1) * xDiff(1) + yDiff(1) * yDiff(1));
        _xi[k] = sp.localArc(i, ab[k]) / arc(); // arc length fraction - eqn (7.61)

        for (int j = 0; j <= o; j++)
            _basis[j][k] = LagrangeBasis::N[o][j](_xi[k]);
    }
}

const int &Element::elementOrder() const { return _elementOrder; };
const int &Element::quadratureOrder() const { return _quadratureOrder; };
const double &Element::arc() const { return _arc; };
const std::vector<double> &Element::t() const { return _t; };
const std::vector<double> &Element::x() const { return _x; };
const std::vector<double> &Element::y() const { return _t; };
const std::vector<double> &Element::dx() const { return _dx; };
const std::vector<double> &Element::dy() const { return _dy; };
const std::vector<double> &Element::J() const { return _J; };
const std::vector<double> &Element::xi() const { return _xi; };
const std::vector<std::vector<double>> &Element::basis() const { return _basis; };

} // namespace bem2D
