#ifndef MECH_H_
#define MECH_H_
#include "imustruct.h"

class MechArrange
{
public:
	MechArrange() {}
	~MechArrange() {}

public:
	NAVIGATION_INFO f_InitializeNav(double time, const Vector& pos,const Vector& vel,const Euler& euler);
	bool f_Mechanization(const IMU_DATA& pre_imuobs, const IMU_DATA& curr_imuobs, const NAVIGATION_INFO& pre_nav, NAVIGATION_INFO& curr_nav);
	IMU_DATA f_InterpolateImu(const IMU_DATA& _preimudata, IMU_DATA& _currimudata, double _time);

private:
	bool	f_AttitudeUpdate(const IMU_DATA& pre_imuobs, const IMU_DATA& curr_imuobs, const NAVIGATION_INFO& pre_nav);
	bool	f_VelocityUpdate(const IMU_DATA& preimuobs, const  IMU_DATA& currimuobs, const NAVIGATION_INFO& pre_nav);
	bool    f_PositionUpdate(const NAVIGATION_INFO& pre_nav, double dt);
	Vector  f_CalVrot(const  IMU_DATA& preimuobs, const  IMU_DATA& currimuobs, const Vector& Vb_k);   /*Calculate the rotational motion of vel*/
	Vector  f_CalVscul(const  IMU_DATA& preimuobs, const  IMU_DATA& currimuobs, const Vector& Vb_k, const Vector& Vb_k_1);   /*Calculate the sculling motion of vel*/
	//NAVIGATION_INFO f_PredictNav(const NAVIGATION_INFO & pre_nav, double dt, const Vector & Vb_k);

private:
	NAVIGATION_INFO temp_current_nav;
};

#endif
