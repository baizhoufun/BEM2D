#include <cmath>
#include "ellipticIntegral.hpp"

namespace numericTools
{

const double EllipticIntegral::ellipticK(double m)
{
    double kp = KP[0], kq = KQ[0];
    double m1 = 1. - m;
    for (int k = 1; k <= 10; k++)
    {
        double m1k = pow(m1, k);
        kp += KP[k] * m1k;
        kq += KQ[k] * m1k;
    }
    return kp - log(m1) * kq;
}

const double EllipticIntegral::ellipticE(double m)
{
    double ep = EP[0], eq = EQ[0];
    double m1 = 1. - m;
    for (int k = 1; k <= 10; k++)
    {
        double m1k = pow(m1, k);
        ep += EP[k] * m1k;
        eq += EQ[k] * m1k;
    }
    return ep - log(m1) * eq;
}

const void EllipticIntegral::ellipticKE(double m, double &K, double &E)
{
    double kp = KP[0], kq = KQ[0], ep = EP[0], eq = EQ[0];
    double m1 = 1. - m;
    for (int k = 1; k <= 10; k++)
    {
        double m1k = pow(m1, k);
        kp += KP[k] * m1k;
        kq += KQ[k] * m1k;
        ep += EP[k] * m1k;
        eq += EQ[k] * m1k;
    }
    K = kp - log(m1) * kq;
    E = ep - log(m1) * eq;
}

const void EllipticIntegral::ellipticKEPQ(double m, double &K, double &E, double &kp, double &kq, double &ep, double &eq)
{
    kp = KP[0], kq = KQ[0], ep = EP[0], eq = EQ[0];
    double m1 = 1. - m;
    for (int k = 1; k <= 10; k++)
    {
        double m1k = pow(m1, k);
        kp += KP[k] * m1k;
        kq += KQ[k] * m1k;
        ep += EP[k] * m1k;
        eq += EQ[k] * m1k;
    }
    K = kp - log(m1) * kq;
    E = ep - log(m1) * eq;
}

const double EllipticIntegral::KP[11] = {1.38629436111989062502E0, 9.65735902811690126535E-2, 3.08851465246711995998E-2, 1.49380448916805252718E-2, 8.79078273952743772254E-3, 6.18901033637687613229E-3, 6.87489687449949877925E-3, 9.85821379021226008714E-3, 7.97404013220415179367E-3, 2.28025724005875567385E-3, 1.37982864606273237150E-4};
const double EllipticIntegral::KQ[11] = {4.99999999999999999821E-1, 1.24999999999870820058E-1, 7.03124996963957469739E-2, 4.88280347570998239232E-2, 3.73774314173823228969E-2, 3.01204715227604046988E-2, 2.39089602715924892727E-2, 1.54850516649762399335E-2, 5.94058303753167793257E-3, 9.14184723865917226571E-4, 2.94078955048598507511E-5};
const double EllipticIntegral::EP[11] = {1.00000000000000000299E0, 4.43147180560990850618E-1, 5.68051945617860553470E-2, 2.18317996015557253103E-2, 1.15688436810574127319E-2, 7.58395289413514708519E-3, 7.77395492516787092951E-3, 1.07350949056076193403E-2, 8.68786816565889628429E-3, 2.50888492163602060990E-3, 1.53552577301013293365E-4};
const double EllipticIntegral::EQ[11] = {0.0, 2.49999999999888314361E-1, 9.37499997197644278445E-2, 5.85936634471101055642E-2, 4.27180926518931511717E-2, 3.34833904888224918614E-2, 2.61769742454493659583E-2, 1.68862163993311317300E-2, 6.50609489976927491433E-3, 1.00962792679356715133E-3, 3.27954898576485872656E-5};

} // namespace numericTools