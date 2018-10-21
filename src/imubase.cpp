#include "imustruct.h"
#include "cmath"
Quaternion::Quaternion()
{
	m_q[0] = m_q[1] = m_q[2] = m_q[3] = 0;
}

Quaternion::Quaternion(double q[4])
{
	m_q[0] = q[0];
	m_q[1] = q[1];
	m_q[2] = q[2];
	m_q[3] = q[3];
}

Quaternion::Quaternion(double scalar, double v1, double v2, double v3)
{
	m_q[0] = scalar;
	m_q[1] = v1;
	m_q[2] = v2;
	m_q[3] = v3;
}

Quaternion::Quaternion(double m_roll,double m_pitch,double m_yaw)
{
	m_q[0] = cos(m_roll / 2)*cos(m_pitch / 2)*cos(m_yaw / 2) + sin(m_roll / 2)*sin(m_pitch / 2)*sin(m_yaw / 2);
	m_q[1] = sin(m_roll / 2)*cos(m_pitch / 2)*cos(m_yaw / 2) - cos(m_roll / 2)*sin(m_pitch / 2)*sin(m_yaw / 2);
	m_q[2] = cos(m_roll / 2)*sin(m_pitch / 2)*cos(m_yaw / 2) + sin(m_roll / 2)*cos(m_pitch / 2)*sin(m_yaw / 2);
	m_q[3] = cos(m_roll / 2)*cos(m_pitch / 2)*sin(m_yaw / 2) - sin(m_roll / 2)*sin(m_pitch / 2)*cos(m_yaw / 2);
}

Euler Quaternion::Quat2Euler()
{
	return DCM2Euler(Quat2DCM());
}

const Quaternion Quaternion::operator*(const Quaternion& q1) const
{
	Quaternion q_res;
	Matrix v1(3, 1, m_q + 1);
	Matrix v2(3, 1, q1.m_q + 1);
	Vector vec1(m_q + 1,3);
	Vector vec2(q1.m_q + 1,3);
	q_res.m_q[0] = *m_q*q1.m_q[0] - ((~v1)*v2)(0, 0);
	Vector iv_cross_res, iv_vector1, iv_vector2;
	iv_cross_res = VectorCross(vec1, vec2);
	iv_vector1 = m_q[0]*vec2;
	iv_vector2 = vec1*q1.m_q[0];
	Vector m_vector = iv_vector1 + iv_vector2 + iv_cross_res;
	q_res.m_q[1] =m_vector[0];
	q_res.m_q[2] =m_vector[1];
	q_res.m_q[3] =m_vector[2];
	return q_res;
}

Quaternion Quaternion::operator*(const Quaternion& q1)
{
	return *(const_cast<Quaternion*>(&(static_cast<const Quaternion&>(*this)*(q1))));
}

const Quaternion Quaternion::operator*(double k) const
{
	Quaternion quat_res;
	quat_res.m_q[0] = m_q[0] * k;
	quat_res.m_q[1] = m_q[1] * k;
	quat_res.m_q[2] = m_q[2] * k;
	quat_res.m_q[3] = m_q[3] * k;
	return quat_res;
}

Quaternion Quaternion::operator*(double k)
{
	return *const_cast<Quaternion*>(&(static_cast<const Quaternion&>(*this)*k));
}

const Quaternion Quaternion::operator+(double k) const
{
	Quaternion quat_res;
	quat_res.m_q[0] = m_q[0] + k;
	quat_res.m_q[1] = m_q[1] + k;
	quat_res.m_q[2] = m_q[2] + k;
	quat_res.m_q[3] = m_q[3] + k;
	return quat_res;
}

Quaternion Quaternion::operator+(double k)
{
	return *const_cast<Quaternion*>(&(static_cast<const Quaternion&>(*this) + k));
}

void Quaternion::Normalization()
{
	double length = 1.0 / sqrt(m_q[0] * m_q[0] + m_q[1] * m_q[1] + m_q[2] * m_q[2] + m_q[3] * m_q[3]);
	*this = (*this)*length;
}

Matrix Quaternion::Quat2DCM()
{
	Matrix dcm_res(3,3);
	dcm_res(0,0) = m_q[0] * m_q[0] + m_q[1] * m_q[1] - m_q[2] * m_q[2] - m_q[3] * m_q[3];
	dcm_res(0,1) = 2 * (m_q[1] * m_q[2] - m_q[0] * m_q[3]);
	dcm_res(0,2) = 2 * (m_q[1] * m_q[3] + m_q[0] * m_q[2]);
	dcm_res(1,0) = 2 * (m_q[1] * m_q[2] + m_q[0] * m_q[3]);
	dcm_res(1,1) = m_q[0] * m_q[0] - m_q[1] * m_q[1] + m_q[2] * m_q[2] - m_q[3] * m_q[3];
	dcm_res(1,2) = 2 * (m_q[2] * m_q[3] - m_q[0] * m_q[1]);
	dcm_res(2,0) = 2 * (m_q[1] * m_q[3] - m_q[0] * m_q[2]);
	dcm_res(2,1) = 2 * (m_q[2] * m_q[3] + m_q[0] * m_q[1]);
	dcm_res(2,2) = m_q[0] * m_q[0] - m_q[1] * m_q[1] - m_q[2] * m_q[2] + m_q[3] * m_q[3];
	return dcm_res;
}

Vector VectorCross(const Vector& vector1, const Vector& vector2)
{
	Vector res(3);
	res[0] = vector1[1] * vector2[2] - vector1[2] * vector2[1];
	res[1] = vector1[2] * vector2[0] - vector1[0] * vector2[2];
	res[2] = vector1[0] * vector2[1] - vector1[1] * vector2[0];
	return res;
}

Matrix CroessProduct(const Vector& vector)
{
	Matrix dcm1(3,3);
	dcm1(0, 0) = 0;
	dcm1(0, 1) = -vector[2];
	dcm1(0, 2) = vector[1];

	dcm1(1, 0) = vector[2];
	dcm1(1, 1) = 0;
	dcm1(1, 2) = -vector[0];

	dcm1(2, 0) = -vector[1];
	dcm1(2, 1) = vector[0];
	dcm1(2, 2) = 0;
	return dcm1;
}

double ModelLength(const Vector& vector)
{
	double sumsq = 0;
	for (size_t i = 0; i < vector.size(); i++)
	{
		sumsq += vector[i] * vector[i];
	}
	return sqrt(sumsq);
}

Vector  CalculateGravity(const Vector& pos)
{
	double gn = 9.7803267715*(1 + 0.0052790414*sin(pos[0])*sin(pos[0]) + 0.0000232719*POW3(sin(pos[0]))*sin(pos[0]));
	gn += (-0.0000030876910891 + 0.0000000043977311*sin(pos[0])*sin(pos[0]))*pos[2];
	gn += 0.0000000000007211*pos[2] * pos[2];
	Vector gn_vec = { 0, 0, gn };
	return gn_vec;
}

Euler DCM2Euler(const Matrix& DCM)
{
	Euler euler1;
	std::get<1>(euler1) = atan2(-DCM(2, 0) , (sqrt(DCM(2, 2) * DCM(2, 2) + DCM(2, 1) * DCM(2, 1))));
	std::get<0>(euler1) = atan2(DCM(2, 1) , DCM(2, 2));
	std::get<2>(euler1)= atan2(DCM(1, 0) , DCM(0, 0));
	return euler1;
}