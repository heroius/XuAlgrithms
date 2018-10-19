  //3BAND.CPP
  //求解一般带型方程组
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  band
  {
  private: 
           int n, h, ih, m;
		   double  **b, **d;
  public:
	       band (int nn, int hh, int mm)
		   {
			   int i;
			   n = nn;  h = hh;  m = mm;
	           ih = 2*h + 1;
	           b = new double*[n];    //动态分配内存空间
	           for (i=0; i<n; i++) b[i] = new double[ih];
	           d = new double*[n];
	           for (i=0; i<n; i++) d[i] = new double[m];
		   }
	       void input ();      //从文件读入带型矩阵数组B以及m组常数向量D
		   void band_gauss ();   //执行列选主元Gauss消去法
           void output ();     //输出结果到文件并显示
		   ~band ()
		   {
			   int i;
			   for (i=0; i<n; i++) { delete [] b[i]; }
			   delete [] b;
			   for (i=0; i<n; i++) { delete [] d[i]; }
			   delete [] d;
		   }
  };

  void band::input ()      //从文件读入带型矩阵数组B以及m组常数向量D
  {
	  int  i, j;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
	  for (i=0; i<n; i++)                       //读入带型矩阵A
		  for (j=0; j<ih; j++)  fin >>b[i][j];
	  for (i=0; i<m; i++)                       //读入m组常数向量D
		  for (j=0; j<n; j++)  fin >>d[j][i];
	  fin.close ();
  }

  void band::band_gauss ()     //执行列选主元Gauss消去法
  { 
	  int ls,k,i,j,is;
      double p,t;
      ls=h;
      for (k=0;k<=n-2;k++)
      { 
		  p=0.0;
          for (i=k;i<=ls;i++)
          { 
			  t=fabs(b[i][0]);
              if (t>p) {p=t; is=i;}
          }
          if (p+1.0==1.0)
          { 
			  cout <<"\n程序工作失败！无解. " <<endl;
			  return;
		  }
          for (j=0;j<=m-1;j++)
          { 
              t=d[k][j]; 
			  d[k][j]=d[is][j]; 
			  d[is][j]=t;
          }
          for (j=0;j<=ih-1;j++)
          { 
              t=b[k][j]; 
			  b[k][j]=b[is][j]; 
			  b[is][j]=t;
          }
          for (j=0;j<=m-1;j++)
              d[k][j]=d[k][j]/b[k][0];
          for (j=1;j<=ih-1;j++)
              b[k][j]=b[k][j]/b[k][0];
          for (i=k+1;i<=ls;i++)
          { 
			  t=b[i][0];
              for (j=0;j<=m-1;j++)
                  d[i][j]=d[i][j]-t*d[k][j];
              for (j=1;j<=ih-1;j++)
                  b[i][j-1]=b[i][j]-t*b[k][j];
              b[i][ih-1]=0.0;
          }
          if (ls!=(n-1)) ls=ls+1;
      }
      p=b[n-1][0];
      if (fabs(p)+1.0==1.0)
      { 
		  cout <<"\n程序工作失败！无解. " <<endl;
		  return;
	  }
      for (j=0;j<=m-1;j++)
          d[n-1][j]=d[n-1][j]/p;
      ls=1;
      for (i=n-2;i>=0;i--)
      { 
		  for (k=0;k<=m-1;k++)
          { 
              for (j=1;j<=ls;j++)
              { 
                  d[i][k]=d[i][k]-b[i][j]*d[i+j][k];
              }
          }
        if (ls!=(ih-1)) ls=ls+1;
      }
  }

  void band::output ()       //输出结果到文件并显示
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
		  fout <<endl;  cout <<endl;
	  }
	  fout.close ();
  }

  void main ()      //主函数
  {
	  band  c(8, 2, 3);
	  c.input ();         //从文件读入带型矩阵数组B以及m组常数向量D
	  c.band_gauss ();    //执行列选主元Gauss消去法
	  c.output ();        //输出结果到文件并显示
  }

