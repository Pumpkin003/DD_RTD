#include "mytool.h"

//WGS84和BDCS椭球的参数，将其设置为全局变量
referenceell WGS84 = { 3.9860050e+14,7.2921151467e-5,6378137.0,(1 / 298.257223563) };
referenceell BDCS = { 3.986004418e+14,7.2921150e-5,6378137.0,(1 / 298.257222101) };

leaps myLeaps;
//跳秒数据设置为全局变量

//选择参考时间
//输入为需要计算的时间、某个测站的在整天的观测数据数组
//输出为参考时间的观测数据
GPSdata ChoiceRef(const obs_epoch& ObsData, const vector<GPSdata>& SatData)
{
	double code = 0;
	for (int j = 0; j < ObsData.obsdata.size(); j++)
	{
		if (ObsData.SatNum[j] == SatData[0].SatNumber && ObsData.SatSystem[j] == SatData[0].SatSys)
		{
			code = ObsData.obsdata[j][0]; break;
		}
	}
	int Tw;
	double Ts;
	JD2GPStime(Date2JD(ObsData.y, ObsData.m, ObsData.d, ObsData.h, ObsData.min, ObsData.sec), Tw, Ts);
	Ts = Ts - code / Vc;
	int size = SatData.size();
	vector<GPSdata> RefData = SatData;
	for (int i = 0; i < size; i++)
	{
		int Toew;
		double Toes;
		JD2GPStime(Date2JD(RefData[i].time[0], RefData[i].time[1], RefData[i].time[2], RefData[i].time[3], RefData[i].time[4], RefData[i].time[5]), Toew, Toes);
		double Tcha = abs(Toes - RefData[i].toe);
		if (Tcha >= 1.0)//RefData[i].time[5] != 0 ||
		{
			RefData.erase(RefData.begin() + i);
			i--;
			size--;
			continue;
		}
	}
	int num = 0;
	double Tk = abs(Ts - RefData[num].toe);
	for (int k = 0; k < size; k++)
	{
		if (Tk > abs(Ts - RefData[k].toe))
		{
			Tk = abs(Ts - RefData[k].toe);
			num = k;
		}
	}

	return RefData[num];
}

//计算某一时刻特定卫星的位置
//输入：观测文件数据块、特定卫星在参考时间的导航文件、钟差、是否考虑相对论和地球自转的影响
//输出：此卫星在这一时刻的地心地固坐标以及钟差
SatCoor CalculateCoor(const obs_epoch& ObsData, const GPSdata& navdata, double dTs, bool Comprehensive)
{
	double miu;
	double omigae;
	double code = 0;
	for (int j = 0; j < ObsData.obsdata.size(); j++)
	{
		if (ObsData.SatNum[j] == navdata.SatNumber && ObsData.SatSystem[j] == navdata.SatSys)
		{
			code = ObsData.obsdata[j][0];
			break;
		}
	}
	int Tw;
	double Ts;
	JD2GPStime(Date2JD(ObsData.y, ObsData.m, ObsData.d, ObsData.h, ObsData.min, ObsData.sec), Tw, Ts);
	Ts = Ts - code / Vc - dTs;
	//JD2GPStime(Date2JD(currenttime[0], currenttime[1], currenttime[2], currenttime[3], currenttime[4], sec), Tw, Ts);
	if (navdata.SatSys == 'G')
	{
		miu = WGS84.miu;
		omigae = WGS84.omigae;
	}
	else
	{
		miu = BDCS.miu;
		omigae = BDCS.omigae;
		Ts -= 14;
	}

	double Toes = navdata.toe;

	double Tk = Ts - Toes;
	double a = navdata.sqrt_a * navdata.sqrt_a;
	double n0 = sqrt(miu / (a * a * a));
	double n = n0 + navdata.Delta_n;

	double Mk = navdata.M0 + n * Tk;
	double Ek;
	{
		double Ek1 = Mk, Ek0 = 0.0;
		int diedai = 0;
		while (abs(Ek1 - Ek0) > 1.0e-12)
		{
			Ek0 = Ek1;
			Ek1 = Mk + navdata.e * sin(Ek0);
			diedai++;
		}
		Ek = Ek1;
	}
	//迭代计算Ek

	double COSVK = (cos(Ek) - navdata.e) / (1 - navdata.e * cos(Ek));
	double Tanv2 = sqrt(1 - navdata.e * navdata.e) / (1 - navdata.e) * tan(Ek / 2);
	double Vk = 2 * atan(Tanv2);
	//double TanVk = 2 * Tanv2 / (1 - Tanv2 * Tanv2);
	//double Vk = atan2(TanVk*COSVK, COSVK);//这里的范围不太确定

	double Uk = Vk + navdata.omega;

	double DUk = navdata.Cuc * cos(2 * Uk) + navdata.Cus * sin(2 * Uk);
	double DRk = navdata.Crc * cos(2 * Uk) + navdata.Crs * sin(2 * Uk);
	double Dik = navdata.Cic * cos(2 * Uk) + navdata.Cis * sin(2 * Uk);

	double U = Uk + DUk;
	double R = a * (1 - navdata.e * cos(Ek)) + DRk;
	double i = navdata.i0 + Dik + navdata.I_DOT * Tk;
	double lamda;
	double X, Y, Z;
	if (abs(navdata.i0 - 0) < 0.2)
	{
		lamda = navdata.OMEGA0 + navdata.OMEGA_DOT * Tk - omigae * navdata.toe;
		double Xgk = R * (cos(U) * cos(lamda) - sin(U) * cos(i) * sin(lamda));
		double Ygk = R * (cos(U) * sin(lamda) + sin(U) * cos(i) * cos(lamda));
		double Zgk = R * sin(U) * sin(i);
		double jiao1 = omigae * Tk;
		double jiao2 = -(5.0 / 180.0 * M_PI);
		X = cos(jiao1) * Xgk + sin(jiao1) * cos(jiao2) * Ygk + sin(jiao1) * sin(jiao2) * Zgk;
		Y = -Xgk * sin(jiao1) + Ygk * cos(jiao1) * cos(jiao2) + sin(jiao2) * Zgk * cos(jiao1);
		Z = -Ygk * sin(jiao2) + Zgk * cos(jiao2);
	}
	else
	{
		lamda = navdata.OMEGA0 + (navdata.OMEGA_DOT - omigae) * Tk - omigae * navdata.toe;

		X = R * (cos(U) * cos(lamda) - sin(U) * cos(i) * sin(lamda));//单位为米
		Y = R * (cos(U) * sin(lamda) + sin(U) * cos(i) * cos(lamda));//单位为米
		Z = R * sin(U) * sin(i);//单位为米
	}
	double Xeci = X, Yeci = Y;
	dTs = navdata.a0 + navdata.a1 * Tk + navdata.a2 * Tk * Tk;
	if (Comprehensive)
	{
		dTs -= 2 * sqrt(miu) * navdata.e * navdata.sqrt_a * sin(Ek) / Vc / Vc;
		Xeci = X * cos(omigae * code / Vc) + Y * sin(omigae * code / Vc);
		Yeci = -X * sin(omigae * code / Vc) + Y * cos(omigae * code / Vc);
	}
	//dTs = dTs * 1000000;//单位为微秒
	return { navdata.SatSys,navdata.SatNumber,Xeci, Yeci, Z,dTs ,navdata.TGD };
}

//计算发射信号时刻各卫星的坐标
//输入：读入的导航文件、读入的观测文件中接收时间和各卫星的伪距值
//输出：发射信号时刻各卫星的地心地固坐标，每一行是观测文件中一个数据块里面各卫星的坐标
void GetSatsendCoor(const vector<vector<GPSdata>>& gpsdata, vector<obs_epoch>& MatchCode, vector<vector<SatCoor>>& Satsendcoor)
{
	for (int i = 0; i < MatchCode.size(); i++)
	{
		vector<SatCoor> currentSat;
		for (int j = 0; j < MatchCode[i].obsdata.size(); j++)
		{
			GPSdata navadata;
			navadata.SatNumber = 0;
			int k = 0;
			for (; k < gpsdata.size(); k++)
			{
				if (MatchCode[i].SatNum[j] == gpsdata[k][0].SatNumber && MatchCode[i].SatSystem[j] == gpsdata[k][0].SatSys)
				{
					navadata = ChoiceRef(MatchCode[i], gpsdata[k]);
					break;
				}
			}
			if (k == gpsdata.size() && navadata.SatNumber == 0)
			{
				MatchCode[i].SatSystem.erase(MatchCode[i].SatSystem.begin() + j);
				MatchCode[i].SatNum.erase(MatchCode[i].SatNum.begin() + j);
				MatchCode[i].obsdata.erase(MatchCode[i].obsdata.begin() + j);
				MatchCode[i].ObsSatNumber--;
				j--;
				continue;
			}
			SatCoor satcoordata;
			satcoordata.CC = 0.0;
			double Dt;
			do
			{
				Dt = satcoordata.CC;
				satcoordata = CalculateCoor(MatchCode[i], navadata, Dt, 0);
			} while (abs(Dt - satcoordata.CC) > 0.000001);//迭代计算发射时间
			currentSat.push_back(CalculateCoor(MatchCode[i], navadata, satcoordata.CC, 1));//考虑发射时间的影响
		}
		Satsendcoor.push_back(currentSat);
	}
	cout << "发射时刻坐标计算完毕" << endl;
	return;
}

//将观测到的伪距取出来的函数
//伪距观测值保存在CodeEpoch中
void GetCode(const vector<obs_epoch>& ObsData, const vector<obs_types>& Obs_Types, vector<obs_epoch>& CodeEpoch)
{
	for (int i = 0; i < ObsData.size(); i++)//遍历每一个数据块
	{
		obs_epoch data = ObsData[i];
		data.obsdata.clear();
		data.SatNum.clear();
		data.SatSystem.clear();
		vector<double> CodeData;
		for (int j = 0; j < ObsData[i].obsdata.size(); j++)//遍历每一颗卫星
		{
			if (ObsData[i].SatSystem[j] == 'G')
			{
				int k = 0;
				for (; k < Obs_Types.size(); k++)
				{
					if (Obs_Types[k].SatSystem == 'G')
						break;
				}
				int C1C = 0, C2W = 0;
				for (; C1C < Obs_Types[k].ObsTypes.size(); C1C++)
				{
					if (Obs_Types[k].ObsTypes[C1C] == "C1C")
						break;
				}
				for (; C2W < Obs_Types[k].ObsTypes.size(); C2W++)
				{
					if (Obs_Types[k].ObsTypes[C2W] == "C2W")
						break;
				}
				CodeData.push_back(ObsData[i].obsdata[j][C1C]);//只存入C1C和C2W的观测值
				CodeData.push_back(ObsData[i].obsdata[j][C2W]);
			}
			else
			{
				int k = 0;
				for (; k < Obs_Types.size(); k++)
				{
					if (Obs_Types[k].SatSystem == 'C')
						break;
				}
				int C2I = 0, C6I = 0;
				for (; C2I < Obs_Types[k].ObsTypes.size(); C2I++)
				{
					if (Obs_Types[k].ObsTypes[C2I] == "C2I")
						break;
				}
				for (; C6I < Obs_Types[k].ObsTypes.size(); C6I++)
				{
					if (Obs_Types[k].ObsTypes[C6I] == "C6I")
						break;
				}
				CodeData.push_back(ObsData[i].obsdata[j][C2I]);
				CodeData.push_back(ObsData[i].obsdata[j][C6I]);
			}
			if (abs(CodeData[0] - 0) > 1 && abs(CodeData[1] - 0) > 1)//统一处理，伪距值如果是零就舍弃
			{
				data.obsdata.push_back(CodeData);
				data.SatNum.push_back(ObsData[i].SatNum[j]);
				data.SatSystem.push_back(ObsData[i].SatSystem[j]);
				CodeData.clear();
			}
			else
			{
				CodeData.clear();
				continue;
			}
		}
		data.SatSystem.push_back('\0');//最后一个字符赋值为尾零
		CodeEpoch.push_back(data);
	}
	return;
}


//寻找共视卫星，将其数据保存再MatchCode中,共视卫星的排列顺序与CodeEpoch2相同
void MatchSat(vector<obs_epoch>& CodeEpoch1, vector<obs_epoch>& CodeEpoch2, vector<obs_epoch>& MatchCode)
{
	for (int i = 0; i < CodeEpoch2.size(); i++)
	{
		obs_epoch TimeData = CodeEpoch2[i];
		TimeData.obsdata.clear();
		TimeData.SatNum.clear();
		TimeData.SatSystem.clear();
		vector<double> CodeData;
		for (int j = 0; j < CodeEpoch2[i].obsdata.size(); j++)
		{
			for (int k = 0; k < CodeEpoch1[i].obsdata.size(); k++)
			{
				if (CodeEpoch2[i].SatSystem[j] == CodeEpoch1[i].SatSystem[k] && CodeEpoch2[i].SatNum[j] == CodeEpoch1[i].SatNum[k])
				{
					//这里存伪距的时候存四个值，分别为测站2、测站1的波段1和测站2、测站1的波段2
					CodeData.push_back(CodeEpoch2[i].obsdata[j][0]);
					CodeData.push_back(CodeEpoch1[i].obsdata[k][0]);
					CodeData.push_back(CodeEpoch2[i].obsdata[j][1]);
					CodeData.push_back(CodeEpoch1[i].obsdata[k][1]);
				}
			}
			if(CodeData.size()>1)
			{
				if (abs(CodeData[0] - 0) > 1 && abs(CodeData[1] - 0) > 1 && abs(CodeData[2] - 0) > 1 && abs(CodeData[3] - 0) > 1)//统一处理，伪距值如果是零就舍弃
				{
					TimeData.obsdata.push_back(CodeData);
					TimeData.SatNum.push_back(CodeEpoch2[i].SatNum[j]);
					TimeData.SatSystem.push_back(CodeEpoch2[i].SatSystem[j]);
					CodeData.clear();
				}
				else
				{
					CodeData.clear();
					continue;
				}
			}
		}
		TimeData.SatSystem.push_back('\0');//最后一个字符赋值为尾零
		TimeData.ObsSatNumber = TimeData.obsdata.size();
		MatchCode.push_back(TimeData);
	}
	cout << "共视卫星匹配完成" << endl;
	return;
}

//删除在两个测站上高度角不全大于15°的卫星
void RemoveAngle15(vector<obs_epoch>& MatchCode, vector<vector<SatCoor>>& Satsendcoor, vector<double>& ApproxXYZ1,
	vector<double>& ApproxXYZ2,vector<vector<vector<double>>>& SatAngle)
{
	for (int i = 0; i < MatchCode.size(); i++)
	{
		double a = 0, f = 0;
		vector<vector<double>> TimeAngle;
		for (int j = 0; j < Satsendcoor[i].size(); j++)
		{
			if (Satsendcoor[i][j].SatSys == 'G')
			{
				a = WGS84.a;
				f = WGS84.f;
			}
			else
			{
				a = BDCS.a;
				f = BDCS.f;
			}
			SatCoor singleSatdata = Satsendcoor[i][j];
			vector<double> SatXYZ = { singleSatdata.X,singleSatdata.Y ,singleSatdata.Z };
			vector<double> RecBLH1 = calculate_BLH(a, f, ApproxXYZ1);
			vector<double> SatNEU1 = calculateNEU(kongjian2chidao(SatXYZ, ApproxXYZ1), RecBLH1);
			vector<double> SatRAEL1 = calculateRAEL(SatNEU1);
			double Ei1 = SatRAEL1[2];
			vector<double> RecBLH2 = calculate_BLH(a, f, ApproxXYZ2);
			vector<double> SatNEU2 = calculateNEU(kongjian2chidao(SatXYZ, ApproxXYZ2), RecBLH2);
			vector<double> SatRAEL2 = calculateRAEL(SatNEU2);
			double Ei2 = SatRAEL2[2];
			if (Ei1 * 180 / M_PI > 15.0 && Ei2 * 180 / M_PI > 15.0)
			{
				TimeAngle.push_back({ Ei1,Ei2 });
				continue;
			}
			else
			{
				MatchCode[i].obsdata.erase(MatchCode[i].obsdata.begin() + j);
				MatchCode[i].SatNum.erase(MatchCode[i].SatNum.begin() + j);
				MatchCode[i].SatSystem.erase(MatchCode[i].SatSystem.begin() + j);
				MatchCode[i].ObsSatNumber--;
				Satsendcoor[i].erase(Satsendcoor[i].begin() + j);
				j--;
			}
		}
		SatAngle.push_back(TimeAngle);
	}
	return;
}


// 用于计算B、l、P等矩阵
void CalculateMatrix(SatCoor& Satsendcoor1, SatCoor& Satsendcoor2, vector<vector<double>>& B,
	vector<vector<double>>& Q, vector<vector<double>>& l, vector<double>& ApproxXYZ1, vector<double>& ApproxXYZ2,
	vector<double>& MatchCode1, vector<double>& MatchCode2, vector<double>Angle, vector<vector<double>>& C, bool system)
{
	double rou_11 = sqrt(pow(Satsendcoor1.X - ApproxXYZ1[0], 2) + pow(Satsendcoor1.Y - ApproxXYZ1[1], 2)
		+ pow(Satsendcoor1.Z - ApproxXYZ1[2], 2));//卫星1和测站1的距离
	double rou_21 = sqrt(pow(Satsendcoor1.X - ApproxXYZ2[0], 2) + pow(Satsendcoor1.Y - ApproxXYZ2[1], 2)
		+ pow(Satsendcoor1.Z - ApproxXYZ2[2], 2));//卫星1和测站2的距离
	double rou_12 = sqrt(pow(Satsendcoor2.X - ApproxXYZ1[0], 2) + pow(Satsendcoor2.Y - ApproxXYZ1[1], 2)
		+ pow(Satsendcoor2.Z - ApproxXYZ1[2], 2));
	double rou_22 = sqrt(pow(Satsendcoor2.X - ApproxXYZ2[0], 2) + pow(Satsendcoor2.Y - ApproxXYZ2[1], 2)
		+ pow(Satsendcoor2.Z - ApproxXYZ2[2], 2));
	//计算权阵
	double Ei = (Angle[0] + Angle[1]) / 2;
	
	vector<double> Ci;
	if (system == 0)
		Ci = { -1,0 };
	else
		Ci = { 0,-1 };
	for (int i = 0; i < Q.size() - 2; i++)
	{
		Ci.push_back(0);
		C[i].push_back(0);
	}
	Ci.push_back(1);
	C.push_back(Ci);

	vector<double> Qi;
	for (int i = 0; i < Q.size(); i++)
	{
		Q[i].push_back(0.0);
		Qi.push_back(0.0);
	}
	Qi.push_back(2);
	Q.push_back(Qi);

	//计算残差矩阵,MatchCode的排序为测站2、测站1
	double DD_Code = MatchCode2[0] - MatchCode1[0] - MatchCode2[1] + MatchCode1[1];
	double DD_Rou = (rou_22 - rou_21 - rou_12 + rou_11);
	vector<double> li = { DD_Code - DD_Rou };
	l.push_back(li);

	//计算设计矩阵B
	double a1 = (ApproxXYZ2[0] - Satsendcoor2.X) / rou_22 - (ApproxXYZ2[0] - Satsendcoor1.X) / rou_21;
	double a2 = (ApproxXYZ2[1] - Satsendcoor2.Y) / rou_22 - (ApproxXYZ2[1] - Satsendcoor1.Y) / rou_21;
	double a3 = (ApproxXYZ2[2] - Satsendcoor2.Z) / rou_22 - (ApproxXYZ2[2] - Satsendcoor1.Z) / rou_21;
	vector<double> Bi = { a1,a2,a3 };
	B.push_back(Bi);
	return;
}

//计算测站的坐标
void GetRecCoor(vector<obs_epoch>& MatchCode, vector<vector<SatCoor>>& Satsendcoor, vector<double>& ApproxXYZ1,
	vector<double>& ApproxXYZ2, vector<vector<double>>& RecCoor, vector<vector<vector<double>>>& SatAngle)
{
	for (int i = 0; i < MatchCode.size(); i++)
	{
		//第一颗卫星坐标
		SatCoor other_Satsendcoor1;
		int num_other = 0;
		SatCoor Satsendcoor1 = Satsendcoor[i][0];

		vector<double> trueXYZ2 = ApproxXYZ2;
		vector<vector<double>> delta = { {0},{0},{0} };
		int diedai = 0;

		do
		{
			diedai++;
			trueXYZ2[0] += delta[0][0], trueXYZ2[1] += delta[1][0];
			trueXYZ2[2] += delta[2][0];
			vector<vector<double>> B;
			vector<vector<double>> Qsd = { {2,0} ,{0,2} };
			vector<vector<double>> C;
			vector<vector<double>> l;
			int flag_othersys = 0;
			for (int j = 1; j < MatchCode[i].obsdata.size(); j++)
			{
				SatCoor Satsendcoor2 = Satsendcoor[i][j];
				if (Satsendcoor2.SatSys == Satsendcoor1.SatSys)
				{
					CalculateMatrix(Satsendcoor1, Satsendcoor2, B, Qsd, l, ApproxXYZ1, trueXYZ2, MatchCode[i].obsdata[0]
						, MatchCode[i].obsdata[j], SatAngle[i][j], C, 0);
				}
				else
				{
					if (flag_othersys == 0)//获取用于计算双差观测值的另一个参考卫星
					{
						flag_othersys = 1;
						other_Satsendcoor1 = Satsendcoor2;
						num_other = j;
					}
					else
					{
						CalculateMatrix(other_Satsendcoor1, Satsendcoor2, B, Qsd, l, ApproxXYZ1, trueXYZ2, MatchCode[i].obsdata[num_other]
							, MatchCode[i].obsdata[j], SatAngle[i][j], C, 1);
					}
				}
			}
			///*
			//* 这里将非差观测值视为等权，利用非差到双差的转换矩阵来定权
			//*/
			//{
			//	P.clear();
			//	for (int i = 0; i < l.size(); i++)
			//	{
			//		vector<double>Pi;
			//		for (int j = 0; j < l.size(); j++)
			//		{
			//			if (j == i)
			//			{
			//				Pi.push_back(l.size() / (2.0 * (l.size() + 1.0)));
			//			}
			//			else
			//				Pi.push_back(-1.0 / (2.0 * (l.size() + 1.0)));
			//		}
			//		P.push_back(Pi);
			//	}
			//}
			vector<vector<double>> CT = transposeMatrix(C);
			vector<vector<double>> Q = multiplyMatrices(C, multiplyMatrices(Qsd, CT));
			vector<vector<double>> P = inverseMatrix(Q);
			vector<vector<double>> BT = transposeMatrix(B);
			vector<vector<double>> BTP = multiplyMatrices(BT, P);
			vector<vector<double>> BTPB = multiplyMatrices(BTP, B);
			vector<vector<double>> Qx = inverseMatrix(BTPB);
			vector<vector<double>> QxBT = multiplyMatrices(Qx, BT);
			vector<vector<double>> QxBTP = multiplyMatrices(QxBT, P);
			delta = multiplyMatrices(QxBTP, l);
		} while ((abs(delta[0][0]) > 0.000001 || abs(delta[1][0]) > 0.000001 || abs(delta[2][0]) > 0.000001) && diedai < 10);
		RecCoor.push_back({ trueXYZ2[0] ,trueXYZ2[1] ,trueXYZ2[2] , double(diedai) });
	}
	cout << "测站2坐标计算完毕" << endl;
	return;
}

//保存误差数据，用于在Matlab里绘图
void MySave(vector<vector<double>> RecCoor)
{
	
	vector<double> TrueXYZ2 = { -.446707539935103E+07 ,0.268301185244666E+07 , -.366700685945280E+07 };//str2
	
	ofstream output_file("product/3RecCoorDataSTR2.txt");
	for (const auto& row : RecCoor)
	{
		for (int i = 0; i < row.size(); i++)
		{
			double data = row[i];
			if (i < 3)
			{
				data -= TrueXYZ2[i];
				//data = abs(data);
			}
			output_file << data << " ";
		}
		output_file << endl;
	}
	output_file.close();

	double a = BDCS.a;
	double f = BDCS.f;
	vector<double> RecBLH = calculate_BLH(a, f, TrueXYZ2);
	vector<vector<double>> NEUerror;
	for (int i = 0; i < RecCoor.size(); i++)
	{
		vector<double> timeXYZ = { RecCoor[i][0],RecCoor[i][1] ,RecCoor[i][2] };
		vector<double> data = calculateNEU(kongjian2chidao(timeXYZ, TrueXYZ2), RecBLH);
		NEUerror.push_back(data);
	}
	//ofstream output_file1("product/NEUerrorDataBDS2023.txt");
	ofstream output_file1("product/3NEUerrorDataSTR2.txt");
	for (const auto& row : NEUerror)
	{
		for (int i = 0; i < row.size(); i++)
		{
			double data = row[i];
			output_file1 << data << " ";
		}
		output_file1 << endl;
	}
	output_file1.close();
	cout << "误差文件保存完毕" << endl;
	return;
}

int main()
{
	vector<ionospheric_corr> ionos_corr;
	vector<time_system_corr> time_corr;
	vector<vector<GPSdata>> gpsdata;//每行都是同一个测站的数据
	vector<vector<SatCoor>> satcoor;//每行都是同一个时间的数据

	string NavPathname = "data/BRDM00DLR_S_20240830000_01D_MN.rnx";
	ReadNavFile(NavPathname, ionos_corr, time_corr, gpsdata, myLeaps);

	vector<obs_epoch> ObsData1, ObsData2;
	vector<double> ApproxXYZ1, ApproxXYZ2;
	vector<double> AntDelta1, AntDelta2;
	vector<obs_types> Obs_Types1, Obs_Types2;
	vector<obs_epoch> CodeEpoch1, CodeEpoch2;
	vector<vector<SatCoor>> Satsendcoor;//每一行是同一时间各卫星的坐标
	string ObsPathname1 = "data/STR100AUS_R_20240830000_01D_30S_MO.rnx";
	string ObsPathname2 = "data/STR200AUS_R_20240830000_01D_30S_MO.rnx";
	ReadObsFile(ObsPathname1, ObsData1, ApproxXYZ1, AntDelta1, Obs_Types1);
	GetCode(ObsData1, Obs_Types1, CodeEpoch1);
	ReadObsFile(ObsPathname2, ObsData2, ApproxXYZ2, AntDelta2, Obs_Types2);
	GetCode(ObsData2, Obs_Types2, CodeEpoch2);

	ApproxXYZ1 = { -.446710334569616E+07 ,0.268303947689997E+07 ,-.366694855820990E+07 };
	double Baseline = sqrt(pow(ApproxXYZ1[0] - ApproxXYZ2[0], 2) + pow(ApproxXYZ1[1] -
		ApproxXYZ2[1], 2) + pow(ApproxXYZ1[2] - ApproxXYZ2[2], 2));
	

	vector<obs_epoch> MatchCode;
	MatchSat(CodeEpoch1, CodeEpoch2, MatchCode);
	GetSatsendCoor(gpsdata, MatchCode, Satsendcoor);//计算所有观测到的卫星发射时刻坐标

	vector<vector<vector<double>>> SatAngle;
	RemoveAngle15(MatchCode, Satsendcoor, ApproxXYZ1, ApproxXYZ2, SatAngle);//计算每颗卫星的两个高度角并移除小于15的卫星
	vector<vector<double>> RecCoor;
	GetRecCoor(MatchCode, Satsendcoor, ApproxXYZ1, ApproxXYZ2, RecCoor, SatAngle);

	MySave(RecCoor);
	return 0;
}