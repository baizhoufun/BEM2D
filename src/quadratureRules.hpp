#ifndef QUADRATURE_HPP
#define QUADRATURE_HPP

namespace numericTools
{
enum class QuadratureType
{
	GAUSS_LEGENDRE,
	GAUSS_LEGENDRE_LOG
};

class QuadratureRules
{
public:
	static const double *abascissa(int order, QuadratureType type);
	static const double *weight(int order, QuadratureType type);

private:
	static const double *qd_GL_x[21];  // Gauss-Legendre abscissa
	static const double *qd_GL_w[21];  // Gauss-Legendre weights
	static const double *qd_LOG_x[21]; // logarithmic-weighted abscissa
	static const double *qd_LOG_w[21]; // logarithmic-weighted weights
	static const double qd_GL_x01[1], qd_GL_x02[2], qd_GL_x03[3], qd_GL_x04[4], qd_GL_x05[5], qd_GL_x06[6], qd_GL_x07[7], qd_GL_x08[8], qd_GL_x09[9], qd_GL_x10[10], qd_GL_x11[11], qd_GL_x12[12], qd_GL_x13[13], qd_GL_x14[14], qd_GL_x15[15], qd_GL_x16[16], qd_GL_x17[17], qd_GL_x18[18], qd_GL_x19[19], qd_GL_x20[20];
	static const double qd_GL_w01[1], qd_GL_w02[2], qd_GL_w03[3], qd_GL_w04[4], qd_GL_w05[5], qd_GL_w06[6], qd_GL_w07[7], qd_GL_w08[8], qd_GL_w09[9], qd_GL_w10[10], qd_GL_w11[11], qd_GL_w12[12], qd_GL_w13[13], qd_GL_w14[14], qd_GL_w15[15], qd_GL_w16[16], qd_GL_w17[17], qd_GL_w18[18], qd_GL_w19[19], qd_GL_w20[20];
	static const double qd_LOG_x01[1], qd_LOG_x02[2], qd_LOG_x03[3], qd_LOG_x04[4], qd_LOG_x05[5], qd_LOG_x06[6], qd_LOG_x07[7], qd_LOG_x08[8], qd_LOG_x09[9], qd_LOG_x10[10], qd_LOG_x11[11], qd_LOG_x12[12], qd_LOG_x13[13], qd_LOG_x14[14], qd_LOG_x15[15], qd_LOG_x16[16], qd_LOG_x17[17], qd_LOG_x18[18], qd_LOG_x19[19], qd_LOG_x20[20];
	static const double qd_LOG_w01[1], qd_LOG_w02[2], qd_LOG_w03[3], qd_LOG_w04[4], qd_LOG_w05[5], qd_LOG_w06[6], qd_LOG_w07[7], qd_LOG_w08[8], qd_LOG_w09[9], qd_LOG_w10[10], qd_LOG_w11[11], qd_LOG_w12[12], qd_LOG_w13[13], qd_LOG_w14[14], qd_LOG_w15[15], qd_LOG_w16[16], qd_LOG_w17[17], qd_LOG_w18[18], qd_LOG_w19[19], qd_LOG_w20[20];
};
} // namespace numericTools

#endif