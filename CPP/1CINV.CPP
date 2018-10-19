//1CINV.CPP
//复矩阵求逆
#include "stdafx.h"
#include  <iostream>
#include  <cmath>
#include  <fstream>

using namespace std;

namespace XuAlgrithms {

	public ref class  cinv
	{
	private:
		int n;
		double  **ar, **ai;
	public:
		cinv(int nn)
		{
			int i;
			n = nn;
			ar = new double*[n];   //动态分配内存空间
			for (i = 0; i < n; i++) ar[i] = new double[n];
			ai = new double*[n];
			for (i = 0; i < n; i++) ai[i] = new double[n];
		}
		void input(double**, double **);    //从文件读入矩阵A的元素
		void inv();      //执行求逆运算
		void output(double**, double**);   //输出逆矩阵到文件
		~cinv()          //释放内存空间
		{
			int i;
			for (i = 0; i < n; i++) { delete[] ar[i]; }
			delete[] ar;
			for (i = 0; i < n; i++) { delete[] ai[i]; }
			delete[] ai;
		}
	};

	void cinv::input(double** arr, double** aii) {
		int i, j;
		for (i = 0; i < n; i++)                       //读入矩阵A
			for (j = 0; j < n; j++)  ar[i][j] = arr[i][j];
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)  ai[i][j] = aii[i][j];
	}

	void cinv::inv()          //执行求逆运算
	{
		int *is, *js, i, j, k;
		double p, q, s, t, d, b;
		is = new int[n];
		js = new int[n];
		for (k = 0; k <= n - 1; k++)
		{
			d = 0.0;
			for (i = k; i <= n - 1; i++)
				for (j = k; j <= n - 1; j++)
				{

					p = ar[i][j] * ar[i][j] + ai[i][j] * ai[i][j];
					if (p > d) { d = p; is[k] = i; js[k] = j; }
				}
			if (d + 1.0 == 1.0)
			{
				delete[] is, js;
				cout << "\n矩阵奇异！没有逆矩阵.  ";
				exit(1);
			}
			if (is[k] != k)
				for (j = 0; j <= n - 1; j++)
				{

					t = ar[k][j]; ar[k][j] = ar[is[k]][j]; ar[is[k]][j] = t;
					t = ai[k][j]; ai[k][j] = ai[is[k]][j]; ai[is[k]][j] = t;
				}
			if (js[k] != k)
				for (i = 0; i <= n - 1; i++)
				{
					t = ar[i][k]; ar[i][k] = ar[i][js[k]]; ar[i][js[k]] = t;
					t = ai[i][k]; ai[i][k] = ai[i][js[k]]; ai[i][js[k]] = t;
				}
			ar[k][k] = ar[k][k] / d; ai[k][k] = -ai[k][k] / d;
			for (j = 0; j <= n - 1; j++)
				if (j != k)
				{
					p = ar[k][j] * ar[k][k]; q = ai[k][j] * ai[k][k];
					s = (ar[k][j] + ai[k][j])*(ar[k][k] + ai[k][k]);
					ar[k][j] = p - q; ai[k][j] = s - p - q;
				}
			for (i = 0; i <= n - 1; i++)
				if (i != k)
				{
					for (j = 0; j <= n - 1; j++)
						if (j != k)
						{
							p = ar[k][j] * ar[i][k]; q = ai[k][j] * ai[i][k];
							s = (ar[k][j] + ai[k][j])*(ar[i][k] + ai[i][k]);
							t = p - q; b = s - p - q;
							ar[i][j] = ar[i][j] - t;
							ai[i][j] = ai[i][j] - b;
						}
				}
			for (i = 0; i <= n - 1; i++)
				if (i != k)
				{
					p = ar[i][k] * ar[k][k]; q = ai[i][k] * ai[k][k];
					s = (ar[i][k] + ai[i][k])*(ar[k][k] + ai[k][k]);
					ar[i][k] = q - p; ai[i][k] = p + q - s;
				}
		}
		for (k = n - 1; k >= 0; k--)
		{
			if (js[k] != k)
				for (j = 0; j <= n - 1; j++)
				{
					t = ar[k][j]; ar[k][j] = ar[js[k]][j]; ar[js[k]][j] = t;
					t = ai[k][j]; ai[k][j] = ai[js[k]][j]; ai[js[k]][j] = t;
				}
			if (is[k] != k)
				for (i = 0; i <= n - 1; i++)
				{
					t = ar[i][k]; ar[i][k] = ar[i][is[k]]; ar[i][is[k]] = t;
					t = ai[i][k]; ai[i][k] = ai[i][is[k]]; ai[i][is[k]] = t;
				}
		}
		delete[] is, js;
	}

	void cinv::output(double** arr, double** aii)       //输出逆矩阵到文件
	{
		int  i, j;
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				arr[i][j] = ar[i][j];
				aii[i][j] = ar[i][j];
			}
		}
	}
}