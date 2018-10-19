  //3LDLE.CPP
  //用分解法求解对称方程组
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  ldle
  {
  private: 
           int n, m;
		   double  **a, **c;
  public:
	       ldle (int nn, int mm)
		   {
			   int i;
			   n = nn;  m = mm;
	           a = new double*[n];   //动态分配内存空间
	           for (i=0; i<n; i++) a[i] = new double[n];
	           c = new double*[n];
	           for (i=0; i<n; i++) c[i] = new double[m];
		   }
	       void input ();    //从文件读入系数矩阵A以及m组常数向量C
		   void a_ldl ();    //执行分解法
           void output ();   //输出结果到文件并显示
		   ~ldle ()
		   {
			   int i;
			   for (i=0; i<n; i++) { delete [] a[i]; }
			   delete [] a;
			   for (i=0; i<n; i++) { delete [] c[i]; }
			   delete [] c;
		   }
  };

  void ldle::input ()      //从文件读入系数矩阵A以及m组常数向量C
  {
	  int  i, j;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
	  for (i=0; i<n; i++)                       //读入矩阵A
		  for (j=0; j<n; j++)  fin >>a[i][j];
	  for (i=0; i<m; i++)                       //读入m组常数向量C
		  for (j=0; j<n; j++)  fin >>c[j][i];
	  fin.close ();
  }

  void ldle::a_ldl ()          //执行分解法
  { 
	  int i,j,k,k1,k2,k3;
      double p;
      if (fabs(a[0][0])+1.0==1.0)
      { 
		  cout <<"\n系数矩阵奇异！无解. "  <<endl;
		  return;
	  }
      for (i=1; i<=n-1; i++)
          a[i][0]=a[i][0]/a[0][0];
      for (i=1; i<=n-2; i++)
      { 
          for (j=1; j<=i; j++)
          { 
              a[i][i]=a[i][i]-a[i][j-1]*a[i][j-1]*a[j-1][j-1];
          }
          p=a[i][i];
          if (fabs(p)+1.0==1.0)
		  { 
		      cout <<"\n系数矩阵奇异！无解. "  <<endl;
		      return;
		  }
          for (k=i+1; k<=n-1; k++)
          { 
              for (j=1; j<=i; j++)
              { 
		          a[k][i]=a[k][i]-a[k][j-1]*a[i][j-1]*a[j-1][j-1];
              }
              a[k][i]=a[k][i]/p;
          }
      }
      for (j=1; j<=n-1; j++)
      { 
          a[n-1][n-1]=a[n-1][n-1]-a[n-1][j-1]*a[n-1][j-1]*a[j-1][j-1];
      }
      p=a[n-1][n-1];
      if (fabs(p)+1.0==1.0)
      { 
		  cout <<"\n系数矩阵奇异！无解. "  <<endl;
		  return;
	  }
      for (j=0; j<=m-1; j++)
      for (i=1; i<=n-1; i++)
      { 
          for (k=1; k<=i; k++)
          { 
              c[i][j]=c[i][j]-a[i][k-1]*c[k-1][j];
          }
      }
      for (i=1; i<=n-1; i++)
      { 
          for (j=i; j<=n-1; j++)
          { 
              a[i-1][j]=a[i-1][i-1]*a[j][i-1];
          }
      }
      for (j=0; j<=m-1; j++)
      { 
          c[n-1][j]=c[n-1][j]/p;
          for (k=1; k<=n-1; k++)
          { 
			  k1=n-k; k3=k1-1;
              for (k2=k1; k2<=n-1; k2++)
              { 
                  c[k3][j]=c[k3][j]-a[k3][k2]*c[k2][j];
              }
              c[k3][j]=c[k3][j]/a[k3][k3];
          }
      }
  }

  void ldle::output ()       //输出结果到文件并显示
  {
	  int  i, j;
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  fout <<endl;  cout <<endl;
	  for (i=0; i<m; i++)
	  {
		  for (j=0; j<n; j++)
		  {
			  fout <<c[j][i] <<"   ";
		      cout <<c[j][i] <<"   ";
		  }
	      fout <<endl;  cout <<endl;
	  }
	  fout.close ();
  }

  void main ()      //主函数
  {
	  ldle  b(5, 2); 
	  b.input ();         //从文件读入系数矩阵A以及m组常数向量C
	  b.a_ldl ();         //执行分解
	  b.output ();        //输出结果到文件并显示
  }

