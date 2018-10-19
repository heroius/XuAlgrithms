  //3CHLK.CPP
  //用Cholesky分解法求解对称正定方程组
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  chlk
  {
  private: 
           int n, m;
		   double  **a, **d;
  public:
	       chlk (int nn, int mm)
		   {
			   int i;
			   n = nn;  m = mm;
	           a = new double*[n];  //动态分配内存空间
	           for (i=0; i<n; i++) a[i] = new double[n];
	           d = new double*[n];
	           for (i=0; i<n; i++) d[i] = new double[m];
		   }
	       void input ();    //从文件读入系数矩阵A以及m组常数向量D
		   void chol ();      //执行Cholesky分解
           void output ();   //输出结果到文件并显示
		   ~chlk ()
		   {
			   int i;
			   for (i=0; i<n; i++) { delete [] a[i]; }
			   delete [] a;
			   for (i=0; i<n; i++) { delete [] d[i]; }
			   delete [] d;
		   }
  };

  void chlk::input ()      //从文件读入系数矩阵A以及m组常数向量D
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
	  for (i=0; i<m; i++)                       //读入m组常数向量D
		  for (j=0; j<n; j++)  fin >>d[j][i];
	  fin.close ();
  }

  void chlk::chol ()       //执行Cholesky分解
  { 
	  int i,j,k;
      if ((a[0][0]+1.0==1.0)||(a[0][0]<0.0))
      { 
		  cout <<"\n程序工作失败！无解. " <<endl;
		  return;
	  }
      a[0][0]=sqrt(a[0][0]);
      for (j=1; j<=n-1; j++)  a[0][j]=a[0][j]/a[0][0];
      for (i=1; i<=n-1; i++)
      { 
          for (j=1; j<=i; j++)
          { 
              a[i][i]=a[i][i]-a[j-1][i]*a[j-1][i];
          }
          if ((a[i][i]+1.0==1.0)||(a[i][i]<0.0))
		  { 
		      cout <<"\n程序工作失败！无解. " <<endl;
		      return;
		  }
          a[i][i]=sqrt(a[i][i]);
          if (i!=(n-1))
          { 
			  for (j=i+1; j<=n-1; j++)
              { 
                  for (k=1; k<=i; k++)
                      a[i][j]=a[i][j]-a[k-1][i]*a[k-1][j];
                  a[i][j]=a[i][j]/a[i][i];
              }
          }
      }
      for (j=0; j<=m-1; j++)
      { 
		  d[0][j]=d[0][j]/a[0][0];
          for (i=1; i<=n-1; i++)
          { 
              for (k=1; k<=i; k++)
                  d[i][j]=d[i][j]-a[k-1][i]*d[k-1][j];
              d[i][j]=d[i][j]/a[i][i];
          }
      }
      for (j=0; j<=m-1; j++)
      { 
          d[n-1][j]=d[n-1][j]/a[n-1][n-1];
          for (k=n-1; k>=1; k--)
          { 
              for (i=k; i<=n-1; i++)
              { 
                  d[k-1][j]=d[k-1][j]-a[k-1][i]*d[i][j];
              }
              d[k-1][j]=d[k-1][j]/a[k-1][k-1];
          }
      }
  }

  void chlk::output ()       //输出结果到文件并显示
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
			  fout <<d[j][i] <<"   ";
		      cout <<d[j][i] <<"   ";
		  }
	      fout <<endl;  cout <<endl;
	  }
	  fout.close ();
  }

  void main ()      //主函数
  {
	  chlk  c(4, 2);
	  c.input ();         //从文件读入系数矩阵A以及m组常数向量D
	  c.chol ();          //执行Cholesky分解
	  c.output ();        //输出结果到文件并显示
  }

