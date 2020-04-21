#ifndef BOUNDARYELEMENT_HPP
#define BOUNDARYELEMENT_HPP
#include <vector>
#include <eigen3/Eigen/Core>
#include "spline/quinticSpline.hpp"
#include "element.hpp"

namespace bem2D
{
class BoundaryElement
{
public:
    typedef Eigen::Matrix<double, 6, 1> Array6d;
    enum class GradType
    {
        Dr,
        Dz
    };

private:
    struct Node
    {
        double x, y, dx, dy, ddx, ddy;
    };
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
    void initialize();
    const Array6d regular(double rp, double zp, int idElement) const;
    const Array6d regularGrad(double rp, double zp, int idElement, GradType type) const;
    //const Array6d regularDz(double rp, double zp, int idElement) const;
    const Array6d axis(double zp, int idElement) const;
    const Array6d singular(double tau, int idElement) const;
    static void assembleMatrix(const BoundaryElement &bem0, const BoundaryElement &bem1, Eigen::MatrixXd &S, Eigen::MatrixXd &D);
    static double interiorField(double rp, double zp, const BoundaryElement &bem, const Eigen::VectorXd &q, const Eigen::VectorXd &p);
    static void swapSDLR(const Eigen::MatrixXd &S, const Eigen::MatrixXd &D, int nSwap, Eigen::MatrixXd &L, Eigen::MatrixXd &R);

private:
    enum class BoundaryRelationType
    {
        Identical,
        JoinedBefore,
        JoinedAfter,
        Disjoint,
    };
    static const BoundaryRelationType checkBoundaryRelation(const BoundaryElement &bem0, const BoundaryElement &bem1);

    void nElement(int i); // never manually set element number, it's determined by spline segments
    void initializeElement();
    void initializeNode();
    void bindNodeToElement();
    void setSourcePoint(double tau, int idElement, double &rp, double &zp) const;

    static void auxFunction_abm(double rp, double zp, double r, double z, double &a, double &b, double &m);
    static double auxFunction_xLogX(double x);
    static double auxFunction_RKE(double P, double Q, double m, double t, double tp);
    static void auxFunction_fKE(double rp, double zp, double r, double z, double dr, double dz, double J, double a, double b,
                                double &f_single_K, double &f_double_K, double &f_double_E);
    static void auxFunction_fKE(double rp, double zp, double r, double z, double dr, double dz, double J, GradType type,
                                double &f_single_K, double &f_double_K, double &f_single_E, double &f_double_E);
};

} // namespace bem2D

#endif