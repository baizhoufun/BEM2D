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
        double x, y, dx, dy, ddx, ddy;
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
    const int &quadratureOrder() const;
    const int &indexShift() const;

    void elementOrder(int i);
    void indexShift(int i);
    void quadratureOrder(int i);

    const spline::Quintic &sp() const;
    spline::Quintic &sp();
    const std::vector<Element> &element() const;
    const std::vector<Node> &node() const;

public:
    void initializeElement();
    void initializeNode();
    void bindNodeToElement();
    void initialize();

private:
    void nElement(int i); // never manually set element number, it's determined by spline segments
    size_t nodeToElement(size_t nodeIndex);
    void setSplineBC(int i, spline::BCType bc0, spline::BCType bc1, double a0 = 0, double b0 = 0, double a1 = 0, double b1 = 0);
};

} // namespace bem2D

#endif