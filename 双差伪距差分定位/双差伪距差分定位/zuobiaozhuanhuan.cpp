#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip> // ��������������ȵ�ͷ�ļ�
using namespace std;

#define M_PI 3.14159265358979323846

/*���ڽ����ǵĿռ�ֱ������תΪ��վ�ĳ��ֱ������
  ���룺���ǵĿռ�ֱ�����ꡢ��վ�Ŀռ�ֱ������
  ���أ����ǵĲ�վ���ֱ������*/
vector<double> kongjian2chidao(vector<double> array1, vector<double> array2)
{
	vector<double> array;
	size_t vectorSize = array1.size();
	for (int i = 0; i < vectorSize; i++)
	{
		array.push_back(array1[i] - array2[i]);
	}
	return array;
}
/*
	���ڼ���ĳ���BLH
	���룺�ռ�ֱ������
	���أ��������
*/
vector<double> calculate_BLH(double a, double f, vector<double> array1)
{
	double X = array1[0]; double Y = array1[1]; double Z = array1[2];
	double e2 = 2 * f - f * f;//����e��ƽ��
	double L = atan2(Y , X);//ʹ��atan2������ֵΪ-pi��pi
	double B0 = 0, N = 0;
	double B1 = atan2(Z ,sqrt(X * X + Y * Y));
	int counter = 0;
	do
	{
		B0 = B1;
		double W = sqrt(1 - e2 * sin(B0) * sin(B0));
		N = a / W;
		B1 = atan2((Z + N * e2 * sin(B0)) , sqrt(X * X + Y * Y));
		counter++;
	} while (abs(B1 - B0) < 0.00001 && counter < 100);
	double B = B1;
	double W = sqrt(1 - e2 * sin(B) * sin(B));
	N = a / W;
	double H = Z / sin(B) - N * (1 - e2);
	vector<double> array = {B, L, H};
	return array;
}

/*
	�������
	���룺��������
	���أ���������ĳ˻�
*/
vector<vector<double>> multiplyMatrices(const vector<vector<double>>& mat1, const vector<vector<double>>& mat2)
{
	int rows1 = mat1.size();
	int cols1 = mat1[0].size();
	int cols2 = mat2[0].size();

	// Check if the matrices can be multiplied
	if (cols1 != mat2.size()) {
		cerr << "Error: Matrix dimensions are not compatible for multiplication." << endl;
		exit(1); // You can handle this error condition as needed
	}

	vector<vector<double>> result(rows1, vector<double>(cols2, 0));

	for (int i = 0; i < rows1; ++i) {
		for (int j = 0; j < cols2; ++j) {
			for (int k = 0; k < cols1; ++k) {
				result[i][j] += mat1[i][k] * mat2[k][j];
			}
		}
	}
	return result;
}

// ��������ת��
vector<vector<double>> transposeMatrix(const vector<vector<double>>& mat) {
	int rows = mat.size();
	int cols = mat[0].size();

	vector<vector<double>> result(cols, vector<double>(rows, 0));

	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			result[j][i] = mat[i][j];
		}
	}

	return result;
}

// LU�ֽ⺯��
void luDecomposition(vector<vector<double>>& A, vector<vector<double>>& L, vector<vector<double>>& U)
{
	int n = A.size();
	L = vector<vector<double>>(n, vector<double>(n, 0));
	U = vector<vector<double>>(n, vector<double>(n, 0));

	for (int i = 0; i < n; ++i) {
		L[i][i] = 1;

		for (int k = i; k < n; ++k) {
			double sum = 0;
			for (int j = 0; j < i; ++j) {
				sum += L[i][j] * U[j][k];
			}
			U[i][k] = A[i][k] - sum;
		}

		for (int k = i + 1; k < n; ++k) {
			double sum = 0;
			for (int j = 0; j < i; ++j) {
				sum += L[k][j] * U[j][i];
			}
			L[k][i] = (A[k][i] - sum) / U[i][i];
		}
	}
	return;
}

// ��� Ly = b
vector<double> solveLowerTriangular(const vector<vector<double>>& L, const vector<double>& b)
{
	int n = L.size();
	vector<double> y(n, 0);

	for (int i = 0; i < n; ++i) {
		double sum = 0;
		for (int j = 0; j < i; ++j) {
			sum += L[i][j] * y[j];
		}
		y[i] = b[i] - sum;
	}

	return y;
}

// ��� Ux = y
vector<double> solveUpperTriangular(const vector<vector<double>>& U, const vector<double>& y) {
	int n = U.size();
	vector<double> x(n, 0);

	for (int i = n - 1; i >= 0; --i) {
		double sum = 0;
		for (int j = i + 1; j < n; ++j) {
			sum += U[i][j] * x[j];
		}
		x[i] = (y[i] - sum) / U[i][i];
	}

	return x;
}

// ���������
vector<vector<double>> inverseMatrix(vector<vector<double>>& A) {
	int n = A.size();
	vector<vector<double>> L, U, invA(n, vector<double>(n, 0));

	// LU�ֽ�
	luDecomposition(A, L, U);

	// �Ե�λ�����ÿһ�н������
	for (int i = 0; i < n; ++i) {
		vector<double> b(n, 0);
		b[i] = 1;

		// �� Ly = b
		vector<double> y = solveLowerTriangular(L, b);

		// �� Ux = y
		vector<double> x = solveUpperTriangular(U, y);

		// �������ֵ�������ĵ� i ��
		for (int j = 0; j < n; ++j) {
			invA[j][i] = x[j];
		}
	}
	return invA;
}


/*
  ����������վ�ĵ�ƽ����ϵ�е�����
  ���룺���ǵ�վ�ĳ������ϵ���꣬��վ�Ĵ������
  ��������ǵ�վ�ĵ�ƽ����ϵ����
*/ 
vector<double> calculateNEU(const vector<double>& XYZ_TES, const vector<double>& BLH)
{
	double B = BLH[0], L = BLH[1];
	vector<vector<double>> XYZ_tes;
	XYZ_tes.push_back(XYZ_TES);
	XYZ_tes = transposeMatrix(XYZ_tes);
	vector<vector<double>> R = {
		{ -sin(B) * cos(L),-sin(L),cos(B) * cos(L)},
		{-sin(B) * sin(L),cos(L),cos(B) * sin(L)},
		{cos(B),0,sin(B)}
	};
	vector<vector<double>> R_ni = inverseMatrix(R);
	vector<vector<double>> result = multiplyMatrices(R_ni, XYZ_tes);
	vector<double> NEU = { result[0][0],result[1][0],result[2][0] };
	return NEU;
}

/*
  ����������Բ�վ�ĸ߶Ƚ��뷽λ��
  ���룺���ǵ�վ�ĵ�ƽ����ϵ����
  �����������Բ�վ�ĸ߶ȽǼ���λ�ǵ�
*/
vector<double> calculateRAEL(const vector<double>& NEU)
{
	double N = NEU[0];
	double E = NEU[1];
	double U = NEU[2];
	double R = sqrt(N * N + E * E + U * U);
	double A = atan2(E, N);
	double EL = asin(U / R);
	return { R,A,EL };
}

