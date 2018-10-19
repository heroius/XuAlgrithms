  //3GAUS.CPP
  //全选主元Gauss消去法求解实系数方程组
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  gaus
  {
  private: 
           int n;
		   double  **a, *b;
  public:
	       gaus (int nn)
		   {
			   int i;
			   n = nn;
	           a = new double*[n];   //动态分配内存空间
	           for (i=0; i<n; i++) a[i] = new double[n];
	           b = new double[n];
		   }
	       void input ();    //从文件读入系数矩阵A以及常数向量B
		   void gauss ();    //执行全选主元Gauss消去法
           void output ();   //输出结果到文件并显示
		   ~gaus ()
		   {
			   int i;
			   for (i=0; i<n; i++) { delete [] a[i]; }
			   delete [] a;
			   delete [] b;
		   }
  };

  void gaus::input ()      //从文件读入系数矩阵A以及常数向量B
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
	  for (i=0; i<n; i++)  fin >>b[i];           //读入常数向量B
	  fin.close ();
  }

  void gaus::gauss ()         //执行全选主元Gauss消去法
  { 
	  int *js,l,k,i,j,is;
      double d,t;
      js = new int[n];
      l=1;
      for (k=0; k<=n-2; k++)
      { 
		  d=0.0;
          for (i=k;i<=n-1;i++)
          for (j=k;j<=n-1;j++)
          { 
			  t=fabs(a[i][j]);
              if (t>d) { d=t; js[k]=j; is=i;}
          }
          if (d+1.0==1.0) l=0;
          else
          { if (js[k]!=k)
              for (i=0;i<=n-1;i++)
              { 
                  t=a[i][k]; 
				  a[i][k]=a[i][js[k]]; 
				  a[i][js[k]]=t;
              }
              if (is!=k)
              { 
				  for (j=k;j<=n-1;j++)
                  { 
                      t=a[k][j]; 
					  a[k][j]=a[is][j]; 
					  a[is][j]=t;
                  }
                  t=b[k]; b[k]=b[is]; b[is]=t;
              }
          }
          if (l==0)
          { 
			  delete [] js;
			  cout <<"\n系数矩阵奇异！无解." <<endl;
              return;
          }
          d=a[k][k];
          for (j=k+1;j<=n-1;j++)
              a[k][j]=a[k][j]/d;
          b[k]=b[k]/d;
          for (i=k+1;i<=n-1;i++)
          { 
			  for (j=k+1;j<=n-1;j++)
                  a[i][j]=a[i][j]-a[i][k]*a[k][j];
              b[i]=b[i]-a[i][k]*b[k];
          }
      }
      d=a[n-1][n-1];
      if (fabs(d)+1.0==1.0)
      { 
		  delete [] js;
		  cout <<"\n系数矩阵奇异！无解." <<endl;
          return;
      }
      b[n-1]=b[n-1]/d;
      for (i=n-2;i>=0;i--)
      { 
		  t=0.0;
          for (j=i+1;j<=n-1;j++)
              t=t+a[i][j]*b[j];
          b[i]=b[i]-t;
      }
      js[n-1]=n-1;
      for (k=n-1;k>=0;k--)
        if (js[k]!=k)
        { 
			t=b[k]; b[k]=b[js[k]]; b[js[k]]=t;
		}
    delete [] js;
  }

  void gaus::output ()       //输出结果到文件并显示
  {
	  int  i;
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  fout <<endl;  cout <<endl;
	  for (i=0; i<n; i++)
	  {
		  fout <<b[i] <<"   ";
		  cout <<b[i] <<"   ";
	  }
	  fout <<endl;  cout <<endl;
	  fout.close ();
  }

  void main ()      //主函数
  {
	  gaus  c(4); 
	  c.input ();         //从文件读系数矩阵A以及常数向量B
	  c.gauss ();         //执行全选主元Gauss消去法
	  c.output ();        //输出结果到文件并显示
  }
