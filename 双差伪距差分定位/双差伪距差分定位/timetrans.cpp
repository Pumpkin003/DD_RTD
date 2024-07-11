#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip> // 包含设置输出精度的头文件
using namespace std;

/*
* 年月日转化为儒略日
* 输入：年、月、日、小时、分钟、秒
* 输出：儒略日
*/
double Date2JD(const int& Y, const int& M, const int& D, const int& H, const int& Min, const double& Sec)
{
	double m = 0, y = 0;
	if (M <= 2)
	{
		y = Y - 1;
		m = M + 12;
	}
	else
	{
		y = Y, m = M;
	}
	double JD = int(365.25 * y) + int(30.6001 * (m + 1)) + D + (H / 1.0 + Min / 60.0 + Sec / 3600.0) / 24 + 1720981.5;
	return JD;
}

/*
* 儒略日转化为年月日
* 输入：儒略日
* 输出：年、月、日、小时
*/
void JD2date(const double& JD, double& hour, int& year, int& month, int& day)
{
	int a = int(JD + 0.5);
	int b = a + 1537;
	int c = int((b - 122.1) / 365.25);
	int d = int(365.25 * c);
	int e = int((b - d) / 30.6001);
	double D = b - d - int(30.6001 * e) + (JD + 0.5 - a);
	int M = e - 1 - 12 * int(e / 14);
	int Y = c - 4715 - int((7 + M) / 10);

	// 计算小时部分
	hour = (JD + 0.5 - a) * 24;

	year = Y;
	month = M;
	day = int(D);
	return;
}

/*
* 儒略日转化为GPS时间
* 输入：儒略日
* 输出：GPS周、GPS秒
*/
void JD2GPStime(const double& JD, int& GPSweek, double& GPSsecond)
{
	GPSweek = int((JD - 2444244.5) / 7);
	GPSsecond = (JD - 2444244.5 - 7 * GPSweek) * 86400.0;
	return;
}