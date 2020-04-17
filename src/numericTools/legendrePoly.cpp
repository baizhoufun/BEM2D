#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif // !M_PI
#include "legendrePoly.hpp"
#include "ellipticIntegral.hpp"

namespace numericTools
{

double LegendrePoly::legendreP(int l, double x, LegendrePolyType type)
{
    switch (type)
    {
    case LegendrePolyType::INTEGER:
    {
        switch (l)
        {
        case 0:
        {
            return 1;
            break;
        }
        case 1:
        {
            return x;
            break;
        }
        case 2:
        {
            return (-1. + 3. * (x * x)) / 2.;
            break;
        }
        case 3:
        {
            return (-3. * x + 5. * (x * x * x)) / 2.;
            break;
        }
        case 4:
        {
            return (3. - 30. * (x * x) + 35. * (x * x * x * x)) / 8.;
            break;
        }
        case 5:
        {
            return (15. * x - 70. * (x * x * x) + 63. * (x * x * x * x * x)) / 8.;
            break;
        }
        default:
        {
            break;
        }
        }
        break;
    }

    case LegendrePolyType::HALF_INTEGER:
    {
        double m = (1. - x) / 2.;
        double E, K;
        EllipticIntegral::ellipticKE(m, K, E);
        switch (l)
        {
        case 1:
        {
            double c = 2. / M_PI;
            double PE = 2.;
            double PK = -1.;
            return c * (PE * E + PK * K);
            break;
        }
        case 3:
        {
            double c = 1. / (3. * M_PI);
            double PE = 16. * x;
            double PK = -2. * (1. + 4. * x);
            return c * (PE * E + PK * K);
            break;
        }
        case 5:
        {
            double c = 1. / (15. * M_PI);
            double PE = 4. * (-9. + 32. * (x * x));
            double PK = -2. * (-9. + 8. * x + 32. * (x * x));
            return c * (PE * E + PK * K);
            break;
        }
        case 7:
        {
            double c = 1. / (105. * M_PI);
            double PE = 64. * x * (-13. + 24. * (x * x));
            double PK = 2. * (25. + 208. * x - 96. * (x * x) - 384. * (x * x * x));
            return c * (PE * E + PK * K);
            break;
        }
        case 9:
        {
            double c = 1. / (315. * M_PI);
            double PE = 588. - 6528. * (x * x) + 8192. * (x * x * x * x);
            double PK = -2. * (147. - 264. * x - 1632. * (x * x) + 512. * (x * x * x) + 2048. * (x * x * x * x));
            return c * (PE * E + PK * K);
            break;
        }
        default:
        {
            break;
        }
        }

        /* code */
        break;
    }

    default:
        return x;
        break;
    }
    return x;
}

} // namespace numericTools