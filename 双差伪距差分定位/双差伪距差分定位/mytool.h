#pragma once

#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip> // 包含设置输出精度的头文件
#include <fstream> //文件流
#include <string>
#include <sstream>
#define M_PI 3.14159265358979323846
using namespace std;

//参考椭球参数
struct referenceell {
	double miu;
	double omigae;//地球自转参数
	double a;
	double f;
};
//电离层误差数据
struct ionospheric_corr {
	char corr_type[5];//读完4位空1格 A4，1X
	double parameters[4];//读12位，小数点后有4位  4D12.4
	char time_mark;//1X,A1
	int SV_ID;//读2位 1X,I2
};
//时间误差数据
struct time_system_corr {
	char type[5]; //A4,1X
	double a0;//读17位，小数点后有10位  D17.10
	double a1;//读16位，小数点后有9位  D16.9
	int t;//参考时间，读6位  1XI6
	int w;//gps周，读4位  1XI4
	char S[6];  //1X,A5,1X
	int U;//读2位  I2,1X
};
//跳秒数据
struct leaps {
	int leapnumber;//Current Number of leap seconds
	int fleap; //future leap second if the week and day numberare in the future
	int Rweek;//Respective week
	int Rday;//Respective day number
};
//GPS导航数据
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
	double TGD;//北斗卫星为TGD1
	double IODC;//4X,4D19.12  北斗卫星为TGD2
	//double Trans_time;
	//double Fit_Interval;
	//double Spare1;
	//double Spare2;//4X,4D19.12
};
////BDS导航数据
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

//卫星导航文件解算数据
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
//用于存储头文件中记录的各系统的观测值数目与类型的结构体
struct obs_types
{
	char SatSystem;//A1
	int ObsNum;//2X,I3
	std::vector<std::string>  ObsTypes;//N(1X,A3)
};//结构体中有vector容器，因此初始化会有问题
//存储一个数据块的结构体
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
* 年月日转化为儒略日
* 输入：年、月、日、小时、分钟、秒
* 输出：儒略日
*/
double Date2JD(const int& Y, const int& M, const int& D, const int& H, const int& Min, const double& Sec);

/*
* 儒略日转化为GPS时间
* 输入：儒略日
* 输出：GPS周、GPS秒
*/
void JD2GPStime(const double& JD, int& GPSweek, double& GPSsecond);
/*
* 儒略日转化为年月日
* 输入：儒略日
* 输出：年、月、日、小时
*/
void JD2date(const double& JD, double & hour, int& year, int& month, int& day);

//读导航文件的函数
void ReadNavFile(const string& pathname, vector<ionospheric_corr>& ionos_corr,
	vector<time_system_corr>& time_corr, vector<vector<GPSdata>>& gpsdata, leaps& myLeaps);

//读观测文件的函数
void ReadObsFile(const string& pathname, vector<obs_epoch>& ObsData,
	vector<double>& ApproxXYZ, vector<double>& AntDelta, vector<obs_types>& Obs_Types);

/*用于将卫星的空间直角坐标转为测站的赤道直角坐标
  输入：卫星的空间直角坐标、测站的空间直角坐标
  返回：卫星的测站赤道直角坐标*/
vector<double> kongjian2chidao(vector<double> array1, vector<double> array2);

/*
	用于计算某点的BLH
	输入：空间直角坐标
	返回：大地坐标
*/
vector<double> calculate_BLH(double a, double f, vector<double> array1);

/*
	矩阵相乘
	输入：两个矩阵
	返回：两个矩阵的乘积
*/
vector<vector<double>> multiplyMatrices(const vector<vector<double>>& mat1, const vector<vector<double>>& mat2);

// 计算矩阵的转置
vector<vector<double>> transposeMatrix(const vector<vector<double>>& mat);

/*
  计算卫星在站心地平坐标系中的坐标
  输入：卫星的站心赤道坐标系坐标，测站的大地坐标
  输出：卫星的站心地平坐标系坐标
*/
vector<double> calculateNEU(const vector<double>& XYZ_TES, const vector<double>& BLH);

/*
  计算卫星相对测站的高度角与方位角
  输入：卫星的站心地平坐标系坐标
  输出：卫星相对测站的高度角及方位角等
*/
vector<double> calculateRAEL(const vector<double>& NEU);

// 计算逆矩阵
vector<vector<double>> inverseMatrix(vector<vector<double>>& A);

////计算年中日
//int day_of_year(int year, int month, int day);
//
///*UNB3模型计算测站对流层延迟
//输入：测站纬度、大地高、年中日、卫星方向高度角
//*/
//double UNB3(double Phi, double H, double DOY, double Ei);