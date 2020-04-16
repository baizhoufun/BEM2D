#ifndef BOUNDARYELEMENT_HPP
#define BOUNDARYELEMENT_HPP
#include <vector>
#include "spline/quinticSpline.hpp"
#include "element.hpp"

namespace bem2D
{
class BoundaryElement
{
public:
    struct Node
    {
        double x, y, dx, dy;
    };

private:
    int _nElement = 0;
    int _elementOrder = 2;
    int _indexShift = 0;
    int _quadratureOrder = 6;
    spline::Quintic _sp;
    std::vector<Element> _element;
    std::vector<Node> _node;

public:
    BoundaryElement();
    ~BoundaryElement();
    const int &nElement() const;
    const int &elementOrder() const;
    const int &indexShift() const;
    const int &quadratureOrder() const;
    const spline::Quintic &sp() const;
    const std::vector<Element> &element() const;
    const std::vector<Node> &node() const;

public:
    void nElement(int i);
    void elementOrder(int i);
    void indexShift(int i);
    void quadratureOrder(int i);
    void setSplineBC(int i, spline::BCType bc0, spline::BCType bc1, double a0 = 0, double b0 = 0, double a1 = 0, double b1 = 0);

private:
    size_t nodeToElement(size_t nodeIndex);
};

} // namespace bem2D

#endif