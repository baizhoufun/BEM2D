#ifndef WRITEEIGEN_HPP
#define WRITEEIGEN_HPP

//#include <opencv2/opencv.hpp>
#include <eigen3/Eigen/Core>
#include <string>

namespace bem2D
{
namespace io
{
class IOEigen
{
public:
    static void waterMark();

    static Eigen::MatrixXd readMatrix(const char *filename, int MAXBUFFSIZE);

    static void write(const std::string &fileName, const Eigen::VectorXd &f);

    static void write(const std::string &fileName, const Eigen::MatrixXd &f);

    static void write(const std::string &fileName, const std::vector<Eigen::VectorXd> &fContainer, int k = 1);

    //    static void img2Mat(const cv::Mat &img, Eigen::VectorXd &b);
    //  static void mat2Img(const Eigen::VectorXd &b, int col, int row, float aspectRatio, float bmin = 0, float bmax = 1);
    static const Eigen::VectorXd std2Eigenvec(const std::vector<double> &b);
};
} // namespace io

} // namespace bem2D

#endif