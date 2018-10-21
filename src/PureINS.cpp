#include "MechArrange.h"
#include "imustruct.h"

//		91620.00500000000	  �Ӳο��ļ�ȡ��     
//		23.13739500000	      113.37136499999	        2.17499937627	        
//		0.00017423244	       -0.00032699472	        0.00024949312	        
//		0.01083286617	       -2.14248721480	      -75.74984266986	
const double rad2m(WGS84_a);
const double starttime = 91620.005;
const double pos[3] = { 23.1373950*deg2rad,  113.37136499999*deg2rad,  2.17499937627 };
const double vel[3] = { 0.00017423244,-0.00032699472,0.00024949312 };
const double euler[3] = { 0.01083286617*deg2rad,	-2.14248721480*deg2rad, -75.74984266986*deg2rad };
struct result
{
	double data[10];
	result()
	{
		memset(data, 0x0, sizeof(double) * 10);
	}
};

void GetResult()
{
	FILE* imufile = fopen("E:\\˶ʿ\\�γ�\\�ߵ�\\��ҵ\\��ҵ3_�γ����_���ߵ�����_���� (1)\\��ҵ4_���ߵ�����\\Data1_PureINS.bin", "rb");
	FILE* fresult = fopen("./result.txt", "w");
	result res;
	while (!feof(imufile))
	{
		fread(&res, sizeof(result), 1, imufile);
		for (size_t i = 0; i < 10; i++)
		{
			fprintf(fresult, "%21.11f\t", res.data[i]);
		}
		fprintf(fresult, "\n");
	}
	fclose(imufile);
	fclose(fresult);
}


void f_ResultPrint(const result& res,const NAVIGATION_INFO& nav, FILE* m_fresult)
{
	double residual[10] = { 0 };
	if (nav.timeofweek!= res.data[0])
	{
		printf("Time Error\n");
	}
	Euler Attitude = DCM2Euler(nav.Cbn);
	for (size_t i = 0; i < 3; i++)
	{
		residual[i + 1] = (nav.BLH[i] - res.data[i + 1] * deg2rad)*rad2m;
		residual[i + 4] = (nav.vNEU[i] - res.data[i + 4]);
	}
	residual[3] = (nav.BLH[2] - res.data[3]);
	residual[7] = (std::get<0>(Attitude)*rad2deg - res.data[7]);
	residual[8] = (std::get<1>(Attitude)*rad2deg - res.data[8]);
	residual[9] = (std::get<2>(Attitude)*rad2deg - res.data[9]);

	fprintf(m_fresult, "%10.3f\t", res.data[0]);
	for (size_t i = 1; i < 10; i++)
	{
		fprintf(m_fresult, "%12.9f\t", residual[i]);
	}
	fprintf(m_fresult, "\n");
	printf( "%10.3f\n", res.data[0]);
	fflush(m_fresult);
}

int main()
{
	GetResult();
	Vector	Pos(pos,3);
	Vector	Vel(vel, 3);
	Euler	mEuler(euler[0], euler[1], euler[2]);
	Quaternion  q(euler[0], euler[1], euler[2]);
	Matrix cbn = q.Quat2DCM();
	Euler temp = DCM2Euler(cbn);
	result res;
	NAVIGATION_INFO m_nav,m_navpre;
	MechArrange mech;
	m_navpre = mech.f_InitializeNav(starttime, Pos, Vel, mEuler);
	IMU_DATA pre_imu, curr_imu;
	FILE* imufile = fopen("E:\\˶ʿ\\�γ�\\�ߵ�\\��ҵ\\��ҵ3_�γ����_���ߵ�����_���� (1)\\��ҵ4_���ߵ�����\\Data1.bin", "rb");
	FILE* refimufile = fopen("E:\\˶ʿ\\�γ�\\�ߵ�\\��ҵ\\��ҵ3_�γ����_���ߵ�����_���� (1)\\��ҵ4_���ߵ�����\\Data1_PureINS.bin", "rb");
	FILE* fresult = fopen("E:\\˶ʿ\\�γ�\\�ߵ�\\��ҵ\\��ҵ3_�γ����_���ߵ�����_���� (1)\\��ҵ4_���ߵ�����\\residual.txt", "w");
	fprintf(fresult,"time(s)\tdeltaB(m)\tdeltaL(m)\tdeltaH(m)\tdeltaN(m/s)\tdeltaE(m/s)\tdeltaD(m/s)\tdeltaR(deg)\tdeltaP(deg)\tdeltaH(deg)\n");
	fread(&pre_imu, sizeof(IMU_DATA), 1, imufile);
	/*����ʱ��*/
	while (!feof(imufile))
	{
		fread(&curr_imu, sizeof(IMU_DATA), 1, imufile);
		if (curr_imu.time>starttime)
		{			
			break;
		}
		pre_imu = curr_imu;
	}
	while (!feof(refimufile))
	{
		fread(&res, sizeof(result), 1, refimufile);
		if (res.data[0]>starttime)
		{
			break;
		}
	}
	/*����*/
	while (!feof(imufile))
	{
		mech.f_Mechanization(pre_imu, curr_imu, m_navpre, m_nav);
		f_ResultPrint(res,m_nav, fresult);
		pre_imu = curr_imu;
		fread(&curr_imu, sizeof(IMU_DATA), 1, imufile);
		fread(&res, sizeof(result), 1, refimufile);
		m_navpre = m_nav;
	}
	fclose(imufile);
	fclose(fresult);
}

