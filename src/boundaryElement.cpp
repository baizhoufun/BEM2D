#include "boundaryElement.hpp"

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

} // namespace bem2D