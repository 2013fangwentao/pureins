#ifndef IMUSTRUCT_H_
#define IMUSTRUCT_H_
#pragma warning(disable:4996)

#ifndef DLL_API
#ifdef WIN32
#define DLL_API __declspec(dllexport)
#else
#define DLL_API
#endif
#endif

#include <valarray>
#include <vector>
#include <sstream>
#include <tuple>
#include "cmatrix"

template <typename _Ty>
inline _Ty POW3(const _Ty& _t)
{
	return _t*_t*_t;
}
template <typename _Ty>
inline _Ty POW2(const _Ty& _t)
{
	return _t*_t;
}
const double  PI(3.141592653589793238462643383279);	
const double  EPSILON(1e-15);
const double  e_5(1.0e-5);
const double  e_6(1.0e-6);
const double  e5(1.0e5);
const double  e6(1.0e6);
const double  deg2rad(PI / 180.0);
const double  rad2deg(180.0 / PI);
const double  dh2rs(PI / 180.0 / 3600.0);
const double  rs2dh(180.0 / PI*3600.0);
const double  oneH(3600.0);
const double  J2(0.00108263);
const double  J4(-2.37091222e-6);
const double  J6(6.08347e-9);
const double  WGS84_AngleRate(7.2921151467E-5);
const double  WGS84_GM(3.986005e14);
const double  WGS84_a(6378137.0);
const double  WGS84_b(6356752.3142);
const double  WGS84_e1(sqrt(POW2(WGS84_a) - POW2(WGS84_b)) / WGS84_a);
const double  WGS84_e2(sqrt(POW2(WGS84_a) - POW2(WGS84_b)) / WGS84_b);
const double  WGS84_f(1 / 298.257223563);

using Vector = std::valarray<double>;
using Matrix = techsoft::matrix<double>;
using Euler  = std::tuple<double, double, double>;//roll pitch yaw
const Vector WIE = { 0, 0, WGS84_AngleRate };



inline double RM(double B)
{
	return WGS84_a*(1 - WGS84_e1*WGS84_e1) / POW3(sqrt(1 - POW2(WGS84_e1*sin(B))));
}

inline double RN(double B)
{
	return WGS84_a / sqrt(1 - POW2(WGS84_e1*sin(B)));
}

inline Vector Wien(double B)
{
	Vector res(3);
	res[0] = WGS84_AngleRate*cos(B);
	res[1] = 0;
	res[2] = -WGS84_AngleRate*sin(B);
	return res;
}
inline Vector Wenn(const Vector& pos,const Vector& vel)
{
	Vector res(3);
	res[0] = vel[1]/(RN(pos[0])+pos[2]);
	res[1] = -vel[0] / (RM(pos[0]) + pos[2]);
	res[2] = -vel[1]*tan(pos[0]) / (RN(pos[0]) + pos[2]);
	return res;
}

/*ÔöÁ¿ÐÍ*/
struct IMU_DATA
{
	double		time;
	double		gyro[3];
	double		acce[3];
	IMU_DATA()
	{
		time = 0.0;
		memset(gyro, 0x0, sizeof(double) * 3);
		memset(acce, 0x0, sizeof(double) * 3);
	}
};


class DLL_API Quaternion
{

public:
	double      m_q[4];

public:
	Quaternion(double q[4]);
	Quaternion();
	Quaternion(double scalar, double v1, double v2, double v3);
	Quaternion(double m_roll, double m_pitch, double m_yaw);
	~Quaternion() {}

public:
	const Quaternion  operator*(const Quaternion& q1) const;
	const Quaternion  operator*(double k) const;
	const Quaternion  operator+(double k) const;
	Quaternion  operator*(const Quaternion& q1);
	Quaternion  operator*(double k);
	Quaternion  operator+(double k);
	Matrix      Quat2DCM();
	Euler       Quat2Euler();
	void        Normalization();
};
DLL_API Vector  CalculateGravity(const Vector& pos);

DLL_API Vector VectorCross(const Vector& vector1, const Vector& vector2);

DLL_API Matrix CroessProduct(const Vector& vector);

DLL_API double ModelLength(const Vector& vector);

DLL_API Euler DCM2Euler(const Matrix& DCM);

struct NAVIGATION_INFO
{
	double  timeofweek;                    /*the time of navigation information*/
	Vector  XYZ;                           /*the body pos in ECEF */
	Vector  BLH;                           /*the body pos in Geographic coordinate system */
	Vector  vXYZ;                          /*the body vel in ECEF */
	Vector  vNEU;                          /*the body vel in navigation  coordinate*/
	Matrix  Cbn;                           /*the DCM of body to navigation coordinate*/
	Euler   Attitude;                      /*the attitude of body in navigation coordinate*/
	Quaternion qbn;                        /*the quaternion of body to navigation coordinate*/

	NAVIGATION_INFO()
	{
		timeofweek = 0.0;
		XYZ.resize(3, 0);
		BLH.resize(3, 0);
		vXYZ.resize(3, 0);
		vNEU.resize(3, 0);
		Cbn.resize(3, 3, 0);
	}
};
#endif
