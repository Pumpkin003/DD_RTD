#include "mytool.h"
using namespace std;

//将字符串转化为电离层误差数据
//手动转换数据，字符数组中都没有尾零
void string2iosc(const string& buf, vector<ionospheric_corr>& ionos_corr)
{
	ionospheric_corr data;
	char buffer[1024] = { 0 };
	char* start = buffer;
	strcpy(buffer, buf.c_str());
	{
		std::copy(start, start + 4, data.corr_type);
		start = start + 5;
		data.corr_type[4] = '\0';
	}// 读完4位空1格 A4，1X
	for (int i = 0; i < 4; i++)
	{
		char dnumber[12];
		std::copy(start, start+12 , dnumber);
		start = start + 12;
		data.parameters[i] = std::strtod(dnumber, nullptr);
	}//读12位，小数点后有4位  4D12.4
	{
		start++;
		data.time_mark=*start; 
	}// 1X, A1
	{
		start = start + 2;
		char inumber[2];
		std::copy(start, start + 2, inumber); 
		data.SV_ID = std::atoi(inumber);
	}//读2位 1X,I2
	ionos_corr.push_back(data);
	return;
}

//将字符串转化为时间系统误差数据
//手动转换数据，字符数组中都没有尾零
void string2timec(const string& buf, vector<time_system_corr>& time_corr)
{
	time_system_corr data;
	char buffer[1024] = { 0 };
	char* start = buffer;
	strcpy(buffer, buf.c_str());
	{
		std::copy(start, start + 4, data.type);
		start = start + 5;
		data.type[4] = '\0';
	}// 读完4位空1格 A4，1X
	{
		char dnumber[17];
		std::copy(start, start + 17, dnumber);
		start = start + 17;
		data.a0 = std::strtod(dnumber, nullptr);
	}//读17位，小数点后有10位  D17.10
	{
		char dnumber[16];
		std::copy(start, start + 16, dnumber);
		start = start + 16;
		data.a1 = std::strtod(dnumber, nullptr);
	}//读16位，小数点后有9位  D16.9
	{
		start++;
		char inumber[6];
		std::copy(start, start + 6, inumber);
		start = start + 6;
		data.t = std::atoi(inumber);
	}//参考时间，读6位  1XI6
	{
		start++;
		char inumber[4];
		std::copy(start, start + 4, inumber);
		start = start + 4;
		data.w = std::atoi(inumber);
	}//gps周，读4位  1XI4
	{
		start++;
		std::copy(start, start + 5, data.S);
		start = start + 6;
		data.S[5] = '\0';
	}//1X,A5,1X
	{
		char inumber[2];
		std::copy(start, start + 2, inumber);
		start = start + 2;
		data.U = std::atoi(inumber);
	}//读2位  I2,1X
	time_corr.push_back(data);
	return;
}

//判断数据块中的每一行有没有缺失（以number个空格为判断条件）
bool hasMoreSpaces(const std::string& str,int number) {
	int spaceCount = 0;
	for (char c : str) 
	{
		if (c == ' ') 
		{
			spaceCount++;
		}
	}
	
	if (spaceCount > number /*|| str.length() < 80*/) {
		return true;
	}
else
	return false;
}

//读取GPS和BDS的导航数据
//并将数据存储在gpsdata中
int ReadGps(ifstream& ifs,string &buf, vector<vector<GPSdata>>& gpsdata)
{
	GPSdata data;
	istringstream iss(buf);
	{
		if (hasMoreSpaces(buf,10))
			return 6;
		iss >> data.SatSys >> data.SatNumber;
		for (int i = 0; i < 6; ++i) {
			iss >> data.time[i];
		}
		iss >> data.a0 >> data.a1 >> data.a2;
	}
	{
		getline(ifs, buf);
		if (hasMoreSpaces(buf,9))
			return 5;
		iss.str(buf); // 
		iss.clear(); // 清除错误标志
		iss >> data.IODE >> data.Crs >> data.Delta_n >> data.M0;
	}
	{
		getline(ifs, buf);
		if (hasMoreSpaces(buf, 9))
			return 4;
		iss.str(buf); // 
		iss.clear(); // 清除错误标志
		iss >> data.Cuc >> data.e >> data.Cus >> data.sqrt_a;
	}
	{
		getline(ifs, buf);
		if (hasMoreSpaces(buf, 9))
			return 3;
		iss.str(buf); // 
		iss.clear(); // 清除错误标志
		iss >> data.toe >> data.Cic >> data.OMEGA0 >> data.Cis;
	}
	{
		getline(ifs, buf);
		if (hasMoreSpaces(buf, 9))
			return 2;
		iss.str(buf); // 
		iss.clear(); // 清除错误标志
		iss >> data.i0 >> data.Crc >> data.omega >> data.OMEGA_DOT;
	}
	{
		getline(ifs, buf);
		if (hasMoreSpaces(buf, 9))
			return 1;
		iss.str(buf); // 
		iss.clear(); // 清除错误标志
		iss >> data.I_DOT >> data.CodeL2 >> data.GPSWeek >> data.L2Flag;
	}
	{
		getline(ifs, buf);
		if (hasMoreSpaces(buf, 9))
			return 0;
		iss.str(buf); // 
		iss.clear(); // 清除错误标志
		iss >> data.SV_accuracy >> data.SV_health >> data.TGD >> data.IODC;
	}
	if (gpsdata.empty())
	{
		vector<GPSdata> DATA;
		DATA.push_back(data);
		gpsdata.push_back(DATA);
	}
	else
	{
		int size = gpsdata.size();
		if (gpsdata[size-1][0].SatSys == data.SatSys && gpsdata[size - 1][0].SatNumber == data.SatNumber)
		{
			gpsdata[size - 1].push_back(data);
		}
		else
		{
			vector<GPSdata> DATA;
			DATA.push_back(data);
			gpsdata.push_back(DATA);
		}
	}
	return 0;
}

//读导航文件的函数
void ReadNavFile(const string& pathname, vector<ionospheric_corr>& ionos_corr,
	vector<time_system_corr>& time_corr, vector<vector<GPSdata>>& gpsdata, leaps& myLeaps)
{
	ifstream ifs;
	ifs.open(pathname, ios::in);//以只读的方式打开文件

	if (!ifs.is_open())
	{
		cout << "文件打开失败了" << endl;
		return;
	}
	string buf;

	while (getline(ifs, buf))
	{
		//cout << ifs.eof() << endl;
		if (buf.find("IONOSPHERIC CORR") != string::npos)
		{
			string2iosc(buf, ionos_corr);
		}
		else if (buf.find("TIME SYSTEM CORR") != string::npos)
		{
			string2timec(buf, time_corr);
		}
		else if (buf.find("END OF HEADER") != string::npos)
		{
			cout << "头文件读取完毕" << endl;
			break;
		}
		else if (buf.find(" LEAP SECONDS") != string::npos)
		{
			sscanf(buf.c_str(), "%d %d %d %d",
				&myLeaps.leapnumber,
				&myLeaps.fleap,
				&myLeaps.Rweek,
				&myLeaps.Rday);
		}
		else
		{
			continue;
		}
	}//读头文件

	while (getline(ifs, buf))
	{
		if (buf.find("G") != string::npos|| buf.find("C") != string::npos)//GPS
		{
			int line=ReadGps(ifs, buf, gpsdata);
			for (int i = -1; i < line; i++)
				getline(ifs, buf);
			
		}
		else if (buf.find("S") != string::npos || buf.find("R") != string::npos)
		{
			for (int i = 0; i < 3; i++)
			{
				getline(ifs, buf);
			}
		}
		else if (buf.find("E") != string::npos || buf.find("I") != string::npos || buf.find("J") != string::npos)
		{
			for (int i = 0; i < 7; i++)
			{
				getline(ifs, buf);
			}
		}
	}
	ifs.close();
	cout << "导航数据读取完毕" << endl;
	return;
}

//读入并存储测站概略坐标及天线相关偏移量
void string2XYZorAnt(const string& buf, vector<double>& data)
{
	char buffer[1024] = { 0 };
	char* start = buffer;
	strcpy(buffer, buf.c_str());
	for (int i = 0; i < 3; i++)
	{
		char dnumber[14];
		std::copy(start, start + 14, dnumber);
		start = start + 14;
		data.push_back(std::strtod(dnumber, nullptr));
	}// 读14位数字,读三次
	return;
}

//读入并储存头文件中记录的各系统的观测值数目与类型
void string2obstype(ifstream& ifs, string& buf, vector<obs_types>& Obs_Types)
{
	obs_types data;
	char buffer[1024] = { 0 };
	char* start = buffer;
	strcpy(buffer, buf.c_str());
	{
		data.SatSystem = *start;
		start = start + 3;
	}//读入导航系统的名称,并空两格
	{
		char inumber[3];
		std::copy(start, start + 3, inumber);
		start = start + 3;
		data.ObsNum = std::atoi(inumber);
	}//读入该导航系统下观测值的数目
	if (data.ObsNum <= 13)
	{
		for (int i = 0; i < data.ObsNum; i++)
		{
			string type(3, ' ');
			start++;
			std::copy(start, start + 3, type.begin());
			start = start + 3;
			data.ObsTypes.push_back(type);
		}//读入一个观测类型并存储
	}
	else
	{
		for (int i = 0; i < 13; i++)
		{
			string type(3, ' ');
			start++;
			std::copy(start, start + 3, type.begin());
			start = start + 3;
			data.ObsTypes.push_back(type);
		}//读入一个观测类型并存储13次
		getline(ifs, buf);
		char buffer[1024] = { 0 };
		start = buffer + 6;
		strcpy(buffer, buf.c_str());
		for (int i = 0; i < data.ObsNum - 13; i++)
		{
			string type(3, ' ');
			start++;
			std::copy(start, start + 3, type.begin());
			start = start + 3;
			data.ObsTypes.push_back(type);
		}//读入剩余的观测类型并存储
	}
	Obs_Types.push_back(data);
	return;
}

//读入数据块的函数
void ReadObsData(ifstream& ifs, string& buf, vector<obs_epoch>& ObsData, const vector<obs_types>& Obs_Types)
{
	obs_epoch data;
	int epochnum = 0;
	while (getline(ifs, buf))
	{
		if (buf.find(">") != string::npos)
		{
			if (epochnum != 0)
			{
				data.SatSystem.push_back('\0');
				ObsData.push_back(data);
				data.SatSystem.clear();
				data.SatNum.clear();
				data.obsdata.clear();
			}
			char buffer[1024] = { 0 };
			char* start = buffer;
			strcpy(buffer, buf.c_str());
			start = start + 2;
			{
				char ynumber[4];
				std::copy(start, start + 4, ynumber);
				start += 5;
				data.y = std::stoi(ynumber);
			}//读入年份
			{
				char inumber[2];
				std::copy(start, start + 2, inumber);
				start += 3;
				data.m = std::stoi(inumber);
				std::copy(start, start + 2, inumber);
				start += 3;
				data.d = std::stoi(inumber);
				std::copy(start, start + 2, inumber);
				start += 3;
				data.h = std::stoi(inumber);
				std::copy(start, start + 2, inumber);
				start += 2;
				data.min = std::stoi(inumber);
			}//读入并存储月、日、时、分
			{
				char dnumber[11];
				std::copy(start, start + 11, dnumber);
				start += 13;
				data.sec = std::strtod(dnumber, nullptr);
			}//读入秒
			char state = *start;
			if (state != '0')
				continue;
			start++;
			data.DataState = state - '0';//读入数据状态
			epochnum++;
			{
				char inumber[3];
				std::copy(start, start + 3, inumber);
				data.ObsSatNumber = std::stoi(inumber);
			}
		}
		else
		{
			char buffer[1024] = { 0 };
			char* start = buffer;
			strcpy(buffer, buf.c_str());
			if (buf.find("G") != string::npos || buf.find("C") != string::npos)
			{
				char Sys = *start; 
				start++;
				int i = 0;
				for (; i < Obs_Types.size(); i++)
				{
					if (Sys == Obs_Types[i].SatSystem)
						break;
				}
				data.SatSystem.push_back(Sys);//读入导航系统的名称

				char inumber[2];
				std::copy(start, start + 2, inumber);
				start += 2;
				int num = std::atoi(inumber);
				data.SatNum.push_back(num);//读入卫星编号

				vector<double> Satdata;
				for (int j = 0; j < Obs_Types[i].ObsNum; j++)
				{
					double singledata = 0;
					char dnumber[14];
					std::copy(start, start + 14, dnumber);
					start += 16;
					singledata = std::strtod(dnumber, nullptr);
					Satdata.push_back(singledata);
				}
				data.obsdata.push_back(Satdata);
			}
			else
				continue;
		}
	}
	data.SatSystem.push_back('\0');
	ObsData.push_back(data);//把最后一个数据放进去
	return;
}

//读观测文件的函数
void ReadObsFile(const string& pathname, vector<obs_epoch>& ObsData, vector<double>& ApproxXYZ, vector<double>& AntDelta, vector<obs_types>& Obs_Types)
{
	ifstream ifs;
	ifs.open(pathname, ios::in);//以只读的方式打开文件

	if (!ifs.is_open())
	{
		cout << "文件打开失败了" << endl;
		return;
	}
	string buf;
	while (getline(ifs, buf))
	{
		//cout << ifs.eof() << endl;
		if (buf.find("APPROX POSITION XYZ") != string::npos)
		{
			string2XYZorAnt(buf, ApproxXYZ);
		}
		else if (buf.find("ANTENNA: DELTA H/E/N") != string::npos)
		{
			string2XYZorAnt(buf, AntDelta);
		}
		else if (buf.find("SYS / # / OBS TYPES") != string::npos)
		{
			string2obstype(ifs, buf, Obs_Types);
		}
		else if (buf.find("END OF HEADER") != string::npos)
		{
			cout << "头文件读取完毕" << endl;
			break;
		}
		else
			continue;
	}//读头文件
	ReadObsData(ifs, buf, ObsData, Obs_Types);//读入数据快
	cout << "观测文件读取完毕" << endl;
	ifs.close();
	return;
}

