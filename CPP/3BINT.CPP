  //3BINT.CPP
  //求解病态方程组
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  bint
  {
  private: 
           int n;
		   double  **a, *b, *x, eps;
		   double  **p, *r, *e;
  public:
	       bint (int nn)
		   {
			   int i;
			   n = nn;
	           a = new double*[n];   //动态分配内存空间
	           for (i=0; i<n; i++) a[i] = new double[n];
               p = new double*[n];
	           for (i=0; i<n; i++) p[i] = new double[n];
               r = new double[n];
               e = new double[n];
	           b = new double[n];
	           x = new double[n];
		   }
	       void input ();  //从文件读入系数矩阵A以及常数向量B与eps
		   void a_bint ();    //执行求解病态方程组的步骤
		   void gauss ();     //全选主元Gauss消去法
           void output ();   //输出结果到文件并显示
		   ~bint ()
		   {
			   int i;
			   for (i=0; i<n; i++) { delete [] a[i]; }
			   delete [] a;
			   for (i=0; i<n; i++) { delete [] p[i]; }
			   delete [] p;
			   delete [] r, e, b, x;
		   }
  };

  void bint::input ()    //从文件读入系数矩阵A以及常数向量B与eps
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

  void bint::a_bint ()       //执行求解病态方程组的步骤
  { 
	  int i,j,k;
      double q,qq;
      i=100;
      for (k=0; k<=n-1; k++)
      for (j=0; j<=n-1; j++)  p[k][j]=a[k][j];
      for (k=0; k<=n-1; k++)  r[k]=b[k]; 
      gauss ();
      for (k=0; k<n; k++)  x[k] = r[k];
      q=1.0+eps;
      while (q>=eps)
      { 
		  if (i==0)  return;
          i=i-1;
		  for (k=0; k<n; k++)
		  {
			  e[k] = 0.0;
			  for (j=0; j<n; j++)  
				  e[k] = e[k] + a[k][j]*x[j];
		  }
          for ( k=0; k<=n-1; k++)  r[k]=b[k]-e[k];
          for ( k=0; k<=n-1; k++)
          for ( j=0; j<=n-1; j++)  p[k][j]=a[k][j];
          gauss ();
          q=0.0;
          for ( k=0; k<=n-1; k++)
          { 
			  qq=fabs(r[k])/(1.0+fabs(x[k]+r[k]));
              if (qq>q) q=qq;
          }
          for ( k=0; k<=n-1; k++)  x[k]=x[k]+r[k];
      }
  }

  void bint::gauss ()         //执行全选主元Gauss消去法
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
			  t=fabs(p[i][j]);
              if (t>d) { d=t; js[k]=j; is=i;}
          }
          if (d+1.0==1.0) l=0;
          else
          { if (js[k]!=k)
              for (i=0;i<=n-1;i++)
              { 
                  t=p[i][k]; 
				  p[i][k]=p[i][js[k]]; 
				  p[i][js[k]]=t;
              }
              if (is!=k)
              { 
				  for (j=k;j<=n-1;j++)
                  { 
                      t=p[k][j]; 
					  p[k][j]=p[is][j]; 
					  p[is][j]=t;
                  }
                  t=r[k]; r[k]=r[is]; r[is]=t;
              }
          }
          if (l==0)
          { 
			  delete [] js;
			  cout <<"\n系数矩阵奇异！无解." <<endl;
              return;
          }
          d=p[k][k];
          for (j=k+1;j<=n-1;j++)
              p[k][j]=p[k][j]/d;
          r[k]=r[k]/d;
          for (i=k+1;i<=n-1;i++)
          { 
			  for (j=k+1;j<=n-1;j++)
                  p[i][j]=p[i][j]-p[i][k]*p[k][j];
              r[i]=r[i]-p[i][k]*r[k];
          }
      }
      d=p[n-1][n-1];
      if (fabs(d)+1.0==1.0)
      { 
		  delete [] js;
		  cout <<"\n系数矩阵奇异！无解." <<endl;
          return;
      }
      r[n-1]=r[n-1]/d;
      for (i=n-2;i>=0;i--)
      { 
		  t=0.0;
          for (j=i+1;j<=n-1;j++)
              t=t+p[i][j]*r[j];
          r[i]=r[i]-t;
      }
      js[n-1]=n-1;
      for (k=n-1;k>=0;k--)
        if (js[k]!=k)
        { 
			t=r[k]; r[k]=r[js[k]]; r[js[k]]=t;
		}
    delete [] js;
  }

  void bint::output ()       //输出结果到文件并显示
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
		  fout <<x[i] <<"   ";
		  cout <<x[i] <<"   ";
	  }
	  fout <<endl;  cout <<endl;
	  fout.close ();
  }

  void main ()      //主函数
  {
	  bint  c(4); 
	  c.input ();         //从文件读入系数矩阵A以及常数向量B与eps
	  c.a_bint ();         //求解病态方程组
	  c.output ();        //输出结果到文件并显示
  }

