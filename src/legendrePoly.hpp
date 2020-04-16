#ifndef LEGENDREPOLY_HPP
#define LEGENDREPOLY_HPP

namespace numericTools
{
enum class LegendrePolyType
{
    INTEGER,
    HALF_INTEGER
};

class LegendrePoly
{
public:
    static double legendreP(int l, double x, LegendrePolyType type); // integer Legendre poly
private:
};
} // namespace numericTools

#endif