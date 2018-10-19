  //4NETN.CPP
  //拟牛顿法求解非线性方程组
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  netn
  {
  private: 
           int n, Max;
		   double  *x, *f, eps, t, h;
		   double  **a, *b;
  public:
	       netn (int nn)
		   {
			   int k;
			   n = nn;
               x = new double[n];  //动态分配内存
	           f = new double[n];
	           b = new double[n];
	           a = new double*[n];
	           for (k=0; k<n; k++) a[k] = new double[n];
		   }
	       void input ();    //由文件读入精度要求eps、控制变量t、
		                 //增量h初值、最大迭代次数Max以及解的一组初值
		   void netn_root ();      //执行拟牛顿法
		   void gauss ();          //全选主元Gauss消去法
           void output ();         //输出非线性方程组的一组解到文件并显示
		   void func (double *);   //计算方程组各左端函数值f
		   ~netn ()
		   {
			   int  k;
			   for (k=0; k<n; k++) { delete [] a[k]; }
			   delete [] a;
			   delete [] x, f, b; 
		   }
  };

  void netn::input ()      //由文件读入精度要求eps、控制变量t、
		              //增量h初值、最大迭代次数Max以及解的一组初值
  {
	  int  k;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
	  fin >>eps;  fin >>t; fin >>h;  fin >>Max;  //读入eps、t、h与Max
	  for (k=0; k<n; k++)  fin>>x[k];     //读入初值
	  fin.close ();
  }

  void netn::netn_root ()         //执行拟牛顿法
  { 
	  int i,j,l;
      double am,z,beta,d;
      l=0; am=1.0+eps;
      while (am>=eps)
      { 
          func (x);
          am=0.0;
          for (i=0; i<=n-1; i++)
          {
			  b[i] = f[i];
			  z=fabs(f[i]);
              if (z>am) am=z;
          }
          if (am>=eps)
          { 
			  l=l+1;
              if (l==Max)
              {
				  cout <<"\n可能不收敛！" <<endl;
                  return;
              }
              for (j=0; j<=n-1; j++)
              { 
				  z=x[j]; x[j]=x[j]+h;
                  func (x);
                  for (i=0; i<=n-1; i++) a[i][j] = f[i];
                  x[j]=z;
              }
              gauss();
              beta=1.0;
              for (i=0; i<=n-1; i++) beta=beta-b[i];
              if (fabs(beta)+1.0==1.0)
              {
				  cout <<"\nBeta=0. 程序工作失败！" <<endl;
                  return;
              }
              d=h/beta;
              for (i=0; i<=n-1; i++) x[i]=x[i]-d*b[i];
              h=t*h;
          }
      }
      cout <<"\n实际迭代次数为: " <<l <<endl;
  }

  void netn::gauss ()         //执行全选主元Gauss消去法
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

  void netn::output ()       //输出非线性方程组的一组解到文件并显示
  {
	  int k;
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  fout <<endl;  cout <<endl;
	  for (k=0; k<n; k++)
	  {
		  fout <<x[k] <<"   ";
		  cout <<x[k] <<"   ";
	  }
	  fout <<endl;  cout <<endl;
	  fout.close ();
  }

  void netn::func (double *x)//计算方程组各左端函数值f
  {
	  f[0] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] - 1.0;
	  f[1] = 2*x[0]*x[0] + x[1]*x[1] - 4*x[2];
	  f[2] = 3*x[0]*x[0] - 4*x[1] + x[2]*x[2];
  }

  void main ()      //主函数
  {
	  netn  root(3);
	  root.input ();      //由文件读入精度要求eps、控制变量t、
		                   //增量h初值、最大迭代次数Max以及解的一组初值
	  root.netn_root ();    //执行拟牛顿法
	  root.output ();      //输出非线性方程组的一组解到文件并显示
  }

//0.785197   0.496611   0.369923   
//实际迭代5次
