  //7GJFQ.CPP
  //双边法
  #include  <iostream>
  #include  <iomanip>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  gjfq
  {
  private: 
	       int n, k;
		   double  t, h, eps, *y, **z;
		   double  *d, *p, *u, *v, *w, *b, *c, *e;
  public:
	       gjfq (int nn, int kk)
		   {
			   int i;
			   n = nn; k = kk;
               y = new double[n];    //动态分配内存
	           d = new double[n];
	           p = new double[n];
	           u = new double[n];
	           v = new double[n];
	           w = new double[n];
	           b = new double[n];
	           c = new double[n];
	           e = new double[n];
	           z = new double*[n];
	           for (i=0; i<n; i++) { z[i] = new double[k]; }
		   }
	       void  input ();         //由文件读入数据成员t,h,eps
		                      //以及n个未知函数在起始点t处的函数值
		   void solution ();    //执行双边法
		   void rkt ();        //执行变步长Runge-Kutta法
           void output ();        //输出k个积分点上
		                          //的未知函数值到文件并显示
		   void func (double,double*,double*);   
		                         //计算微分方程组中各方程右端函数值
		   ~gjfq ()
		   { 
			   int i;
			   for (i=0; i<n; i++) { delete [] z[i]; }
			   delete [] z;
			   delete [] y, d, p, u, v, w, b, c, e;
		   }
  };

  void gjfq::input ()     //由文件读入数据成员t,h,eps
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

  void gjfq::solution ()   //执行双边法
  { 
      int i,j;
      double a,qq;
      for (i=0; i<=n-1; i++)   w[i] = y[i]; 
      a=t;
      func (t,y,d);
      for (j=0; j<=n-1; j++) u[j]=d[j];
      rkt ();
      t=a+h;
      func (t,y,d);
      for (j=0; j<=n-1; j++)
      { z[j][0]=y[j]; v[j]=d[j];}
      for (j=0; j<=n-1; j++)
      { 
		  p[j]=-4.0*z[j][0]+5.0*w[j]+2.0*h*(2.0*v[j]+u[j]);
          y[j]=p[j];
      }
      t=a+2.0*h;
      func (t,y,d);
      for (j=0; j<=n-1; j++)
      { 
		  qq=2.0*h*(d[j]-2.0*v[j]-2.0*u[j])/3.0;
          qq=qq+4.0*z[j][0]-3.0*w[j];
          z[j][1]=(p[j]+qq)/2.0;
          y[j]=z[j][1];
      }
      for (i=2; i<=k-1; i++)
      { 
		  t=a+i*h;
          func (t,y,d);
          for (j=0; j<=n-1; j++)
          { u[j]=v[j]; v[j]=d[j];}
          for (j=0; j<=n-1; j++)
          { 
			  qq=-4.0*z[j][i-1]+5.0*z[j][i-2];
              p[j]=qq+2.0*h*(2.0*v[j]+u[j]);
              y[j]=p[j];
          }
          t=t+h;
          func (t,y,d);
          for (j=0; j<=n-1; j++)
          { 
			  qq=2.0*h*(d[j]-2.0*v[j]-2.0*u[j])/3.0;
              qq=qq+4.0*z[j][i-1]-3.0*z[j][i-2];
              y[j]=(p[j]+qq)/2.0;
              z[j][i]=y[j];
          }
      }
	  t = a;
  }

  void gjfq::rkt ()  //执行变步长Runge-Kutta法
  { 
	  int m,i,j,k;
      double hh,pp,dt,x,tt,q,a[4];
      hh=h; m=1; pp=1.0+eps; x=t;
      for (i=0; i<=n-1; i++) c[i]=y[i];
      while (pp>=eps)
      { 
		  a[0]=hh/2.0; a[1]=a[0]; a[2]=hh; a[3]=hh;
          for (i=0; i<=n-1; i++)
          { p[i]=y[i]; y[i]=c[i];}
          dt=h/m; t=x;
          for (j=0; j<=m-1; j++)
          { 
			  func (t,y,d);
              for (i=0; i<=n-1; i++) 
              { b[i]=y[i]; e[i]=y[i];}
              for (k=0; k<=2; k++)
              { 
				  for (i=0; i<=n-1; i++)
                  { 
					  y[i]=e[i]+a[k]*d[i];
                      b[i]=b[i]+a[k+1]*d[i]/3.0;
                  }
                  tt=t+a[k];
                  func (tt,y,d);
              }
              for (i=0; i<=n-1; i++)
                  y[i]=b[i]+hh*d[i]/6.0;
              t=t+dt;
          }
          pp=0.0;
          for (i=0; i<=n-1; i++)
          { 
			  q=fabs(y[i]-p[i]);
              if (q>pp) pp=q;
          }
          hh=hh/2.0; m=m+m;
      }
  }

  void gjfq::output ()   //输出k个积分点上
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

  void gjfq::func (double t, double y[], double d[])
  {
	  d[0] = -y[1];
	  d[1] = y[0];
  }

  void main ()      //主函数
  {
	  gjfq  s(2, 10);           //创建对象
	  s.input ();        //由文件读入数据成员t,h,eps
		                 //以及n个未知函数在起始点t处的函数值
	  s.solution ();    //执行双边法
	  s.output (); //输出k个积分点上的未知函数值到文件并显示
  }

