  //7RUNG_KUTA_1.CPP
  //变步长Runge-Kutta法
  #include  <iostream>
  #include  <iomanip>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  rung_kuta_1
  {
  private: 
	       int n, k;
		   double  t, h, eps, *y, **z;
		   double  *d, *g, *b, *c, *e;
  public:
	       rung_kuta_1 (int nn, int kk)
		   {
			   int i;
			   n = nn; k = kk;
               y = new double[n];      //动态分配内存
	           d = new double[n];
	           g = new double[n];
	           b = new double[n];
	           c = new double[n];
	           e = new double[n];
	           z = new double*[n];
	           for (i=0; i<n; i++) { z[i] = new double[k]; }
		   }
	       void  input ();         //由文件读入数据成员t,h,eps
		                      //以及n个未知函数在起始点t处的函数值
		   void solution ();    //执行变步长Runge-Kutta法
           void output ();        //输出k个积分点上
		                          //的未知函数值到文件并显示
		   void func (double,double*,double*);   
		                         //计算微分方程组中各方程右端函数值
		   ~rung_kuta_1 ()
		   { 
			   int i;
			   for (i=0; i<n; i++) { delete [] z[i]; }
			   delete [] z;
			   delete [] y, d, g, b, c, e;
		   }
  };

  void rung_kuta_1::input ()     //由文件读入数据成员t,h,eps
		                      //以及n个未知函数在起始点t处的函数值
  {
	  int  i;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
      fin >>t >>h >>eps;        //读入t,h,eps
      for (i=0; i<n; i++)  fin >>y[i];  
	                            //读入n个未知函数在起始点t处的函数值	  
	  fin.close ();
  }

  void rung_kuta_1::solution ()  //执行变步长Runge-Kutta法
  { 
	  int m,i,j,kk,u;
      double hh,p,dt,x,tt,q,a[4],t1;
	  t1 = t;
      for (u=0; u<k; u++)
	  {
		  t = t1 + u*h;
          hh=h; m=1; p=1.0+eps; x=t;
          for (i=0; i<=n-1; i++) c[i]=y[i];
          while (p>=eps)
		  { 
			  a[0]=hh/2.0; a[1]=a[0]; a[2]=hh; a[3]=hh;
              for (i=0; i<=n-1; i++)
			  { g[i]=y[i]; y[i]=c[i]; }
              dt=h/m; t=x;
              for (j=0; j<=m-1; j++)
			  { 
				  func (t,y,d);
                  for (i=0; i<=n-1; i++) 
				  { b[i]=y[i]; e[i]=y[i]; }
                  for (kk=0; kk<=2; kk++)
				  { 
					  for (i=0; i<=n-1; i++)
					  { 
						  y[i]=e[i]+a[kk]*d[i];
                          b[i]=b[i]+a[kk+1]*d[i]/3.0;
					  }
                     tt=t+a[kk];
                     func (tt,y,d);
				  }
                  for (i=0; i<=n-1; i++)  y[i]=b[i]+hh*d[i]/6.0;
                  t=t+dt;
			  }
              p=0.0;
             for (i=0; i<=n-1; i++)
			 { 
				 q=fabs(y[i]-g[i]);
                 if (q>p) p=q;
			 }
             hh=hh/2.0; m=m+m;
		  }
          for (i=0; i<=n-1; i++)  z[i][u] = y[i];
	  }
	  t = t1;
  }

  void rung_kuta_1::output ()   //输出k个积分点上
	                      //的未知函数值到文件并显示
  {
	  int i, j;
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
      for (i=0; i<k; i++)
	  {
		  cout <<"t = " <<t+(i+1)*h <<endl;
		  for (j=0; j<n; j++)
		  {
			  fout <<z[j][i] <<"   ";
			  cout <<"y(" <<j <<")=" <<setw(12) <<z[j][i] <<"  ";
		  }
		  fout <<endl;  cout <<endl;
	  }
	  fout.close ();
  }

  void rung_kuta_1::func (double t, double y[], double d[])
  {
	  d[0] = y[1];
	  d[1] = -y[0];
  }

  void main ()      //主函数
  {
	  rung_kuta_1  s(2, 10);           //创建对象
	  s.input ();        //由文件读入数据成员t,h,eps
		                 //以及n个未知函数在起始点t处的函数值
	  s.solution ();    //执行变步长Runge-Kutta法
	  s.output (); //输出k个积分点上的未知函数值到文件并显示
  }

