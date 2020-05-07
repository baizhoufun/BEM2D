#ifndef ELEMENT_HPP
#define ELEMENT_HPP
#include <vector>
#include "spline/quinticSpline.hpp"

namespace bem2D
{
class Element
{
private:
    double _arc = 0;
    int _elementOrder = 2;    // element order (linear or quadratic basis) 	//int _nqd = 6;
    int _quadratureOrder = 2; // element order (linear or quadratic basis) 	//int _nqd = 6;
    std::vector<double> _t, _x, _y, _dx, _dy, _J, _xi;
    std::vector<std::vector<double>> _basis;
    std::vector<int> _nodeIndex;

public:
    Element();
    Element(const Element &e);
    ~Element();
    void initialize(const spline::Quintic &sp, int iSplineSegment, int nElementOrder, int nQuadratureOrder = 6, double const *abacissa = nullptr);

    const double &arc() const;
    const int &elementOrder() const;
    const int &quadratureOrder() const;
    const std::vector<double> &t() const;
    const std::vector<double> &x() const;
    const std::vector<double> &y() const;
    const std::vector<double> &dx() const;
    const std::vector<double> &dy() const;
    const std::vector<double> &J() const;
    const std::vector<double> &xi() const;
    const std::vector<std::vector<double>> &basis() const;
    const std::vector<int> &nodeIndex() const;
    std::vector<int> &setNodeIndex();
};

} // namespace bem2D

#endif