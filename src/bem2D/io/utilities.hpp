#ifndef UTILITIES_HPP
#define UTILITIES_HPP
#include <string>

namespace bem2D
{
namespace io
{
class Utilities
{
public:
    static double startTime;
    static double endTime;
    static void waterMark();
    static void tic(bool output = false);
    static void toc(bool output = false);
    static double tictoc(bool output = false);
    static const std::string padZero(int a, int total);

}; // namespace io
} // namespace io
} // namespace bem2D

#endif