#ifndef ELLIPTICINTEGRAL_HPP
#define ELLIPTICINTEGRAL_HPP

namespace numericTools
{

class EllipticIntegral
{
public:
    static const double ellipticK(double m);
    static const double ellipticE(double m);
    static const void ellipticKE(double m, double &K, double &E);
    static const void ellipticKEPQ(double m, double &K, double &E, double &KP, double &KQ, double &EP, double &EQ);

private:
    static const double KP[11], KQ[11], EP[11], EQ[11];
};

}; // namespace numericTools
#endif