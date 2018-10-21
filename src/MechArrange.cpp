#include "MechArrange.h"


NAVIGATION_INFO MechArrange::f_InitializeNav(double time,const Vector & pos, const Vector & vel, const Euler & euler)
{
	NAVIGATION_INFO m_nav;
	m_nav.BLH = pos;
	m_nav.vNEU = vel;
	m_nav.Attitude = euler;
	Quaternion temp(std::get<0>(euler), std::get<1>(euler), std::get<2>(euler));
	m_nav.qbn = temp;
	m_nav.Cbn = temp.Quat2DCM();
	m_nav.timeofweek = time;
	return m_nav;
}

bool MechArrange::f_Mechanization(const IMU_DATA & pre_imuobs, const IMU_DATA & curr_imuobs, const NAVIGATION_INFO & pre_nav, NAVIGATION_INFO & curr_nav)
{
	if (curr_imuobs.time-pre_imuobs.time<1e-6)
	{
		curr_nav = pre_nav;
		return true;
	}
	f_VelocityUpdate(pre_imuobs, curr_imuobs, pre_nav);
	f_PositionUpdate(pre_nav, (-pre_imuobs.time + curr_imuobs.time));
	f_AttitudeUpdate(pre_imuobs, curr_imuobs, pre_nav);
	temp_current_nav.timeofweek = curr_imuobs.time;
	curr_nav = temp_current_nav;
	return true;
}

bool MechArrange::f_AttitudeUpdate(const IMU_DATA & pre_imuobs, const IMU_DATA & curr_imuobs, const NAVIGATION_INFO & pre_nav)
{
	/*---------------------------------------------------calculate the rotate vector----------------------------------------------------------*/
	Vector curr_phi = { 0, 0, 0 };
	Vector delta_thet_k = { curr_imuobs.gyro[0], curr_imuobs.gyro[1], curr_imuobs.gyro[2] };
	Vector delta_thet_k_1 = { pre_imuobs.gyro[0], pre_imuobs.gyro[1], pre_imuobs.gyro[2] };
	Vector temp;
	temp = VectorCross(delta_thet_k_1, delta_thet_k);
	temp = temp*(1.0 / 12.0);
	curr_phi = delta_thet_k + temp;

	/*-------------------------------------------calculate the quaternion for k-1 to k in body coorinate--------------------------------------*/
	Quaternion qbk_bk_1;
	temp = curr_phi*0.5;
	qbk_bk_1.m_q[0] = cos(ModelLength(temp));
	double coeff_v = sin(ModelLength(temp)) / ModelLength(temp);
	Vector m_vector = curr_phi*0.5*coeff_v;
	qbk_bk_1.m_q[1] = m_vector[0];
	qbk_bk_1.m_q[2] = m_vector[1];
	qbk_bk_1.m_q[3] = m_vector[2];

	/*--------------------------------------------calculate the quaternion for k-1 to k in navigtion coorinate------------------------------------------------*/
	Vector current_Vb = { curr_imuobs.acce[0], curr_imuobs.acce[1], curr_imuobs.acce[2] };
	NAVIGATION_INFO predict;
	predict.vNEU = 0.5*(temp_current_nav.vNEU + pre_nav.vNEU);
	predict.BLH = 0.5*(temp_current_nav.BLH + pre_nav.BLH);
	Vector zeta =(Wien(predict.BLH[0])+ Wenn(predict.BLH, predict.vNEU))*(curr_imuobs.time - pre_imuobs.time);
	Quaternion qek_ek_1;
	temp = zeta*0.5;
	qek_ek_1.m_q[0] = cos(ModelLength(temp));
	coeff_v = sin(ModelLength(temp)) / ModelLength(temp);
	m_vector = temp*(-coeff_v);
	qek_ek_1.m_q[1] = m_vector[0];
	qek_ek_1.m_q[2] = m_vector[1];
	qek_ek_1.m_q[3] = m_vector[2];

	/*-----------------------------------------calculate quaternion from body to ecef currently------------------------------------------------*/
	Quaternion qbk_ek_1 = pre_nav.qbn*qbk_bk_1;
	Quaternion qbk_ek = qek_ek_1*qbk_ek_1;
	qbk_ek.Normalization();

	/*--------------------------------------------------update attitude information-------------------------------------------------------------*/
	temp_current_nav.qbn = qbk_ek;
	temp_current_nav.Cbn = qbk_ek.Quat2DCM();
	return false;
}

bool MechArrange::f_VelocityUpdate(const  IMU_DATA& preimuobs, const  IMU_DATA& currimuobs, const NAVIGATION_INFO& pre_nav)
{
	Vector pre_deltaVb ={ preimuobs.acce[0], preimuobs.acce[1], preimuobs.acce[2] };
	Vector current_Vb = { currimuobs.acce[0], currimuobs.acce[1], currimuobs.acce[2] };
	/*--------------------------------------------------------------calculate the gravity and Coriolis correction-------------------------------------*/
	double dt = currimuobs.time - preimuobs.time;
	NAVIGATION_INFO predict = pre_nav;

	Vector gn_vec = CalculateGravity(predict.BLH);
	Vector wien = Wien(predict.BLH[0]);
	Vector wenn = Wenn(predict.BLH, predict.vNEU);
	Vector omega_n = (wien* 2.0 + wenn);
	Vector delta_Vgcor = (gn_vec-VectorCross(omega_n, predict.vNEU))*dt;

	/*-------------------------------------------------------------------calculate the C from time of k-1 to k---------------------------------------*/
	Vector 	temp = (wenn + wien)*0.5*dt;
	Matrix Cek_1_ek = (CroessProduct(temp))*-1.0;
	Cek_1_ek(0, 0) += 1;
	Cek_1_ek(1, 1) += 1;
	Cek_1_ek(2, 2) += 1;

	/*----------------------------------------------------calculate rotational correction and sculling correction-------------------------------------*/
	Vector Vrot_vec = f_CalVrot(preimuobs, currimuobs, current_Vb);
	Vector Vscul_vec = f_CalVscul(preimuobs, currimuobs, current_Vb, pre_deltaVb);

	/*-----------------------------------------------------calculate the increment of vel in ecef-----------------------------------------------------*/
	Vector delta_Vb = current_Vb + Vrot_vec + Vscul_vec;
	Vector delta_Ve = Cek_1_ek*pre_nav.Cbn*delta_Vb;

	/*----------------------------------------------------------update velocity information-----------------------------------------------------------*/
	temp_current_nav.vNEU = pre_nav.vNEU + delta_Ve + delta_Vgcor;
	return true;
}

bool MechArrange::f_PositionUpdate(const NAVIGATION_INFO & pre_nav, double dt)
{
	temp_current_nav.BLH[2] = pre_nav.BLH[2] - 0.5*(pre_nav.vNEU[2] + temp_current_nav.vNEU[2])*dt;
	double aveh = 0.5*(temp_current_nav.BLH[2] + pre_nav.BLH[2]);
	temp_current_nav.BLH[0] = pre_nav.BLH[0] + 0.5*(pre_nav.vNEU[0] + temp_current_nav.vNEU[0])*dt / (RM(pre_nav.BLH[0]) + aveh);
	double aveB = 0.5*(temp_current_nav.BLH[0] + pre_nav.BLH[0]);
	temp_current_nav.BLH[1] = pre_nav.BLH[1] + 0.5*(pre_nav.vNEU[1] + temp_current_nav.vNEU[1])*dt / ((RN(aveB)+aveh)*cos(aveB));
	return true;
}

IMU_DATA MechArrange::f_InterpolateImu(const IMU_DATA& _preimudata, IMU_DATA& _currimudata, double _time)
{
	IMU_DATA inter_imu_data=_currimudata;
	//inter_imu_data.time = -1;
	if ((_currimudata.time - _time  > (1e-6)) && (_time - _preimudata.time>1e-6))
	{
		inter_imu_data.time = _time;
		double t_interp = (inter_imu_data.time - _preimudata.time) / (_currimudata.time - _preimudata.time);
		{
			for (int i = 0; i < 3; i++)
			{
				inter_imu_data.gyro[i] = _currimudata.gyro[i] * t_interp;
				inter_imu_data.acce[i] = _currimudata.acce[i] * t_interp;
				_currimudata.gyro[i] -= inter_imu_data.gyro[i];
				_currimudata.acce[i] -= inter_imu_data.acce[i];

			}
		}
	}
	return inter_imu_data;
}

Vector MechArrange::f_CalVrot(const  IMU_DATA& preimuobs, const  IMU_DATA& currimuobs, const Vector& Vb_k)
{
	Vector delta_thet = { currimuobs.gyro[0], currimuobs.gyro[1], currimuobs.gyro[2] };
	Vector delta_Vrot;
	delta_Vrot = VectorCross(delta_thet, Vb_k);
	delta_Vrot = delta_Vrot*0.5;
	return delta_Vrot;
}

Vector MechArrange::f_CalVscul(const  IMU_DATA& preimuobs, const  IMU_DATA& currimuobs, const Vector& Vb_k, const Vector& Vb_k_1)
{
	Vector thet_k_1 = { preimuobs.gyro[0], preimuobs.gyro[1], preimuobs.gyro[2] };
	Vector thet_k = { currimuobs.gyro[0], currimuobs.gyro[1], currimuobs.gyro[2] };
	Vector temp1, temp2, temp3;
	temp1 = VectorCross(thet_k_1, Vb_k);
	temp2 = VectorCross(Vb_k_1, thet_k);
	temp3 = temp1 + temp2;
	Vector delta_Vscul = temp3*(1.0 / 12.0);
	return delta_Vscul;
}

//NAVIGATION_INFO MechArrange::f_PredictNav(const NAVIGATION_INFO & pre_nav, double dt, const Vector & Vb_k)
//{
//	NAVIGATION_INFO  predict_nav;
//	Vector gn_vec = CalculateGravity(pre_nav.BLH);
//	Vector wien = Wien(pre_nav.BLH[0]);
//	Vector wenn = Wenn(pre_nav.BLH, pre_nav.vNEU);
//	Vector omega_n = (wien* 2.0 + wenn);
//	Vector delta_Vgcor = (gn_vec - VectorCross(omega_n, pre_nav.vNEU))*dt;
//	dt *= 0.5;
//	predict_nav.vNEU = pre_nav.vNEU + pre_nav.Cbn*(Vb_k*0.5) + delta_Vgcor*0.5;
//	predict_nav.BLH[2] = pre_nav.BLH[2] - (pre_nav.vNEU[2])*dt;
//	double aveh = predict_nav.BLH[2];
//	predict_nav.BLH[0] = pre_nav.BLH[0] + (pre_nav.vNEU[0])*dt / (RM(pre_nav.BLH[0]) + aveh);
//	double aveB = pre_nav.BLH[0];
//	predict_nav.BLH[1] = pre_nav.BLH[1] + (pre_nav.vNEU[1])*dt / ((RN(aveB) + aveh)*cos(aveB));
//	return predict_nav;
//}
