#include "utilities.hpp"
#include <iostream>
#include <omp.h>

namespace bem2D
{
namespace io
{
void Utilities::waterMark()
{
    std::cout << "  ________   ________   _________   ___   _____ ______    ________   ___           \n";
    std::cout << " |\\   __  \\ |\\   __  \\ |\\___   ___\\|\\  \\ |\\   _ \\  _   \\ |\\   __  \\ |\\  \\          \n";
    std::cout << " \\ \\  \\|\\  \\\\ \\  \\|\\  \\\\|___ \\  \\_|\\ \\  \\\\ \\  \\\\\\__\\ \\  \\\\ \\  \\|\\  \\\\ \\  \\         \n";
    std::cout << "  \\ \\  \\\\\\  \\\\ \\   ____\\    \\ \\  \\  \\ \\  \\\\ \\  \\\\|__| \\  \\\\ \\   __  \\\\ \\  \\        \n";
    std::cout << "   \\ \\  \\\\\\  \\\\ \\  \\___|     \\ \\  \\  \\ \\  \\\\ \\  \\    \\ \\  \\\\ \\  \\ \\  \\\\ \\  \\____   \n";
    std::cout << "    \\ \\_______\\\\ \\__\\         \\ \\__\\  \\ \\__\\\\ \\__\\    \\ \\__\\\\ \\__\\ \\__\\\\ \\_______\\ \n";
    std::cout << "     \\|_______| \\|__|          \\|__|   \\|__| \\|__|     \\|__| \\|__|\\|__| \\|_______| \n";
    std::cout << "                               _______    ___  ___   ___                        \n";
    std::cout << "                              |\\  ___ \\  |\\  \\|\\  \\ |\\  \\                       \n";
    std::cout << "                              \\ \\   __/| \\ \\  \\\\\\  \\\\ \\  \\                      \n";
    std::cout << "                               \\ \\  \\_|/__\\ \\   __  \\\\ \\  \\                     \n";
    std::cout << "                                \\ \\  \\_|\\ \\\\ \\  \\ \\  \\\\ \\  \\____                \n";
    std::cout << "                                 \\ \\_______\\\\ \\__\\ \\__\\\\ \\_______\\              \n";
    std::cout << "                                  \\|_______| \\|__|\\|__| \\|_______|              \n";
    std::cout << std::endl;
}

// MATLAB-like time counter
void Utilities::tic(bool output)
{
    if (output)
    {
        std::cout << "tic\n";
    };
    startTime = omp_get_wtime();
}
void Utilities::toc(bool output)
{
    endTime = omp_get_wtime();
    if (output)
    {
        std::cout << "toc ... " << endTime - startTime << "\n";
    };
}

const std::string Utilities::padZero(int a, int total)
{

    int dataLength = total;
    unsigned int number_of_digits = 0;
    do
    {
        ++number_of_digits;
        dataLength /= 10;
    } while (dataLength);
    int number_of_zeros = number_of_digits - std::to_string(a).length();
    return std::string(number_of_zeros, '0').append(std::to_string(a));
}

double Utilities::tictoc(bool output)
{

    if (output)
    {
        std::cout << "tictoc = " << endTime - startTime << std::endl;
    };
    return endTime - startTime;
    return 0;
}
double Utilities::startTime = 0.;
double Utilities::endTime = 0.;
} // namespace io
} // namespace bem2D
