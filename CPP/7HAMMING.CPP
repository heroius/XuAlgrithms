  //7HAMMING.CPP
  //Hamming方法
  #include  <iostream>
  #include  <iomanip>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  hamming
  {
  private: 
	       int n, k;
		   double  t, h, eps, *y, **z;
		   double  *d, *b, *u, *v, *w, *g, *bb;
  public:
	       hamming (int nn,  int kk)
		   {
			   int i;
			   n = nn; k = kk;
               y = new double[n];    //动态分配内存
	           d = new double[n];
	           b = new double[4*n];
	           u = new double[n];
	           v = new double[n];
	           w = new double[n];
	           g = new double[n];
	           bb = new double[n];
	           z = new double*[n];
	           for (i=0; i<n; i++) { z[i] = new double[k]; }
		   }
	       void  input ();         //由文件读入数据成员t,h,eps
		                      //以及n个未知函数在起始点t处的函数值
		   void solution ();    //执行Hamming方法
		   void rkt();          //执行四阶Runge-Kutta公式
           void output ();        //输出k个积分点上
		                          //的未知函数值到文件并显示
		   void func (double,double*,double*);   
		                         //计算微分方程组中各方程右端函数值
		   ~hamming ()
		   { 
			   int i;
			   for (i=0; i<n; i++) { delete [] z[i]; }
			   delete [] z;
			   delete [] y, d, b, u, v, w, g, bb;
		   }
  };

  void hamming::input ()     //由文件读入数据成员t,h,eps
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

  void hamming::solution ()  //执行Hamming方法
  { 
      int i,j,m;
      double a,q,p;
      a=t;
      for (i=0; i<=n-1; i++) w[i]=y[i];
      func (t,y,d);
      for (i=0; i<=n-1; i++) b[i]=d[i];
      for (i=0; i<3; i++)
          if (i<=k-1)
		  { 
			  t=a+i*h;
              rkt();
              for (m=0; m<=n-1; m++) z[m][i]=y[m];
			  t = a + (i+1)*h;
              func (t,y,d);
              for (m=0; m<=n-1; m++) b[(i+1)*n+m]=d[m];
		  }
      for (i=0; i<=n-1; i++) u[i]=0.0;
      for (i=3; i<k; i++)
      { 
		  for (j=0; j<=n-1; j++)
          { 
			  q=2.0*b[3*n+j]-b[n+n+j]+2.0*b[n+j];
			  p=z[j][i-4];
			  if (i==3) p = w[j];
              y[j]=p + 4.0*h*q/3.0;
          }
          for (j=0; j<=n-1; j++)
              y[j]=y[j]+112.0*u[j]/121.0;
          t=a+(i+1)*h;
          func (t,y,d);
          for (j=0; j<=n-1; j++)
          { 
			  q=9.0*z[j][i-1]-z[j][i-3];
              q=(q+3.0*h*(d[j]+2.0*b[3*n+j]-b[n+n+j]))/8.0;
              u[j]=q-y[j];
              z[j][i]=q-9.0*u[j]/121.0;
              y[j]=z[j][i];
              b[n+j]=b[n+n+j];
              b[n+n+j]=b[n+n+n+j];
          }
          func (t,y,d);
          for (m=0; m<=n-1; m++) b[3*n+m]=d[m];
      }
      t = a;
  }

  void hamming::rkt()  //执行四阶Runge-Kutta公式
  { 
	  int m,i,j,k;
      double hh,p,dt,x,tt,q,a[4];
      hh=h; m=1; p=1.0+eps; x=t;
      for (i=0; i<=n-1; i++) v[i]=y[i];
      while (p>=eps)
      { 
		  a[0]=hh/2.0; a[1]=a[0]; a[2]=hh; a[3]=hh;
          for (i=0; i<=n-1; i++)
          { g[i]=y[i]; y[i]=v[i]; }
          dt=h/m; t=x;
          for (j=0; j<=m-1; j++)
          { 
			  func (t,y,d);
              for (i=0; i<=n-1; i++) 
              { bb[i]=y[i]; u[i]=y[i];}
              for (k=0; k<=2; k++)
              { 
				  for (i=0; i<=n-1; i++)
                  { 
					  y[i]=u[i]+a[k]*d[i];
                      bb[i]=bb[i]+a[k+1]*d[i]/3.0;
                  }
                  tt=t+a[k];
                  func (tt,y,d);
              }
              for (i=0; i<=n-1; i++)
                  y[i]=bb[i]+hh*d[i]/6.0;
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
  }

  void hamming::output ()   //输出k个积分点上
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

  void hamming::func (double t, double y[], double d[])
  {
	  d[0] = y[1];
	  d[1] = -y[0];
	  d[2] = y[2];
  }

  void main ()      //主函数
  {
	  hamming  s(3, 10);           //创建对象
	  s.input ();        //由文件读入数据成员t,h,eps
		                 //以及n个未知函数在起始点t处的函数值
	  s.solution ();    //执行Hamming方法
	  s.output (); //输出k个积分点上的未知函数值到文件并显示
  }

