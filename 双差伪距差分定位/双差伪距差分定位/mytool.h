#pragma once

#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip> // ��������������ȵ�ͷ�ļ�
#include <fstream> //�ļ���
#include <string>
#include <sstream>
#define M_PI 3.14159265358979323846
using namespace std;

//�ο��������
struct referenceell {
	double miu;
	double omigae;//������ת����
	double a;
	double f;
};
//������������
struct ionospheric_corr {
	char corr_type[5];//����4λ��1�� A4��1X
	double parameters[4];//��12λ��С�������4λ  4D12.4
	char time_mark;//1X,A1
	int SV_ID;//��2λ 1X,I2
};
//ʱ���������
struct time_system_corr {
	char type[5]; //A4,1X
	double a0;//��17λ��С�������10λ  D17.10
	double a1;//��16λ��С�������9λ  D16.9
	int t;//�ο�ʱ�䣬��6λ  1XI6
	int w;//gps�ܣ���4λ  1XI4
	char S[6];  //1X,A5,1X
	int U;//��2λ  I2,1X
};
//��������
struct leaps {
	int leapnumber;//Current Number of leap seconds
	int fleap; //future leap second if the week and day numberare in the future
	int Rweek;//Respective week
	int Rday;//Respective day number
};
//GPS��������
struct GPSdata {
	char SatSys;//A1
	int SatNumber;// I2.2
	int time[6];// Satellite Time of Clock 1X,I4+5(1X,I2.2) 
	double a0;
	double a1;
	double a2;//3D19.12
	double IODE;
	double Crs;
	double Delta_n;
	double M0;//4X,4D19.12
	double Cuc;
	double e;
	double Cus;
	double sqrt_a;//4X,4D19.12
	double toe;
	double Cic;
	double OMEGA0;
	double Cis;//4X,4D19.12
	double i0;
	double Crc;
	double omega;
	double OMEGA_DOT;//4X,4D19.12
	double I_DOT;
	double CodeL2;
	double GPSWeek;
	double L2Flag;//4X,4D19.12
	double SV_accuracy;
	double SV_health;
	double TGD;//��������ΪTGD1
	double IODC;//4X,4D19.12  ��������ΪTGD2
	//double Trans_time;
	//double Fit_Interval;
	//double Spare1;
	//double Spare2;//4X,4D19.12
};
////BDS��������
//struct BDSdata {
//	char SatSys;//A1
//	int SatNumber;// I2.2
//	int time[6];// Satellite Time of Clock 1X,I4+5(1X,I2.2) 
//	double a0;
//	double a1;
//	double a2;//3D19.12
//	double AODE;
//	double Crs;
//	double Delta_n;
//	double M0;//4X,4D19.12
//	double Cuc;
//	double e;
//	double Cus;
//	double sqrt_a;//4X,4D19.12
//	double toe;
//	double Cic;
//	double OMEGA0;
//	double Cis;//4X,4D19.12
//	double i0;
//	double Crc;
//	double omega;
//	double OMEGA_DOT;//4X,4D19.12
//	double I_DOT;
//	double spare1;
//	double BDTWeek;
//	double spare2;//4X,4D19.12
//	double SV_accuracy;
//	double SatH1;
//	double TGD1;
//	double TGD2;//4X,4D19.12
//	double Trans_time;
//	double AODC;
//	double Spare3;
//	double Spare4;//4X,4D19.12
//};

//���ǵ����ļ���������
struct SatCoor
{
	char SatSys;//A1
	int SatNumber;// I2.2
	double X;
	double Y;
	double Z;
	double CC;//clock corr;
	double TGD1;
};
//���ڴ洢ͷ�ļ��м�¼�ĸ�ϵͳ�Ĺ۲�ֵ��Ŀ�����͵Ľṹ��
struct obs_types
{
	char SatSystem;//A1
	int ObsNum;//2X,I3
	std::vector<std::string>  ObsTypes;//N(1X,A3)
};//�ṹ������vector��������˳�ʼ����������
//�洢һ�����ݿ�Ľṹ��
struct obs_epoch
{
	int y;
	int m;
	int d;
	int h;
	int min;
	double sec;
	int DataState;
	int ObsSatNumber;
	vector<char> SatSystem;
	vector<int> SatNum;
	vector<vector<double>> obsdata;
};

static int Vc = 299792458;

/*
* ������ת��Ϊ������
* ���룺�ꡢ�¡��ա�Сʱ�����ӡ���
* �����������
*/
double Date2JD(const int& Y, const int& M, const int& D, const int& H, const int& Min, const double& Sec);

/*
* ������ת��ΪGPSʱ��
* ���룺������
* �����GPS�ܡ�GPS��
*/
void JD2GPStime(const double& JD, int& GPSweek, double& GPSsecond);
/*
* ������ת��Ϊ������
* ���룺������
* ������ꡢ�¡��ա�Сʱ
*/
void JD2date(const double& JD, double & hour, int& year, int& month, int& day);

//�������ļ��ĺ���
void ReadNavFile(const string& pathname, vector<ionospheric_corr>& ionos_corr,
	vector<time_system_corr>& time_corr, vector<vector<GPSdata>>& gpsdata, leaps& myLeaps);

//���۲��ļ��ĺ���
void ReadObsFile(const string& pathname, vector<obs_epoch>& ObsData,
	vector<double>& ApproxXYZ, vector<double>& AntDelta, vector<obs_types>& Obs_Types);

/*���ڽ����ǵĿռ�ֱ������תΪ��վ�ĳ��ֱ������
  ���룺���ǵĿռ�ֱ�����ꡢ��վ�Ŀռ�ֱ������
  ���أ����ǵĲ�վ���ֱ������*/
vector<double> kongjian2chidao(vector<double> array1, vector<double> array2);

/*
	���ڼ���ĳ���BLH
	���룺�ռ�ֱ������
	���أ��������
*/
vector<double> calculate_BLH(double a, double f, vector<double> array1);

/*
	�������
	���룺��������
	���أ���������ĳ˻�
*/
vector<vector<double>> multiplyMatrices(const vector<vector<double>>& mat1, const vector<vector<double>>& mat2);

// ��������ת��
vector<vector<double>> transposeMatrix(const vector<vector<double>>& mat);

/*
  ����������վ�ĵ�ƽ����ϵ�е�����
  ���룺���ǵ�վ�ĳ������ϵ���꣬��վ�Ĵ������
  ��������ǵ�վ�ĵ�ƽ����ϵ����
*/
vector<double> calculateNEU(const vector<double>& XYZ_TES, const vector<double>& BLH);

/*
  ����������Բ�վ�ĸ߶Ƚ��뷽λ��
  ���룺���ǵ�վ�ĵ�ƽ����ϵ����
  �����������Բ�վ�ĸ߶ȽǼ���λ�ǵ�
*/
vector<double> calculateRAEL(const vector<double>& NEU);

// ���������
vector<vector<double>> inverseMatrix(vector<vector<double>>& A);

////����������
//int day_of_year(int year, int month, int day);
//
///*UNB3ģ�ͼ����վ�������ӳ�
//���룺��վγ�ȡ���ظߡ������ա����Ƿ���߶Ƚ�
//*/
//double UNB3(double Phi, double H, double DOY, double Ei);