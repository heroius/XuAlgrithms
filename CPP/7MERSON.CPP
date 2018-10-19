  //7MERSON.CPP
  //变步长Merson方法
  #include  <iostream>
  #include  <iomanip>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  merson
  {
  private: 
	       int n, k;
		   double  t, h, eps, *y, **z;
		   double  *a, *b, *c, *d, *u, *v;
  public:
	       merson (int nn, int kk)
		   {
			   int i;
			   n = nn; k = kk;
               y = new double[n];    //动态分配内存
	           a = new double[n];
	           b = new double[n];
	           c = new double[n];
	           d = new double[n];
	           u = new double[n];
	           v = new double[n];
	           z = new double*[n];
	           for (i=0; i<n; i++) { z[i] = new double[k]; }
		   }
	       void  input ();         //由文件读入数据成员t,h,eps
		                      //以及n个未知函数在起始点t处的函数值
		   void solution ();    //执行变步长Merson方法
           void output ();        //输出k个积分点上
		                          //的未知函数值到文件并显示
		   void func (double,double*,double*);   
		                         //计算微分方程组中各方程右端函数值
		   ~merson ()
		   { 
			   int i;
			   for (i=0; i<n; i++) { delete [] z[i]; }
			   delete [] z;
			   delete [] y, a, b, c, d, u,v;
		   }
  };

  void merson::input ()     //由文件读入数据成员t,h,eps
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

  void merson::solution ()   //执行变步长Merson方法
  { 
	  int i,j,m,nn;
      double aa,bb,x,hh,p,dt,t0,qq;
      aa=t;
      for (i=0; i<k; i++)
      { 
		  x=aa+i*h; nn=1; hh=h;
          for (j=0; j<=n-1; j++) u[j]=y[j];
          p=1.0+eps;
          while (p>=eps)
          { 
			  for (j=0; j<=n-1; j++)
              { v[j]=y[j]; y[j]=u[j]; }
              dt=h/nn; t=x;
              for (m=0; m<=nn-1; m++)
              { 
				  func (t,y,d);
                  for (j=0; j<=n-1; j++)
                  { a[j]=d[j]; y[j]=y[j]+hh*d[j]/3.0; }
                  t0=t+hh/3.0;
                  func (t0,y,d);
                  for (j=0; j<=n-1; j++)
                  { b[j]=d[j]; y[j]=y[j]+hh*(d[j]-a[j])/6.0; }
                  func (t0,y,d);
                  for (j=0; j<=n-1; j++)
                  { 
					  b[j]=d[j];
                      bb=(d[j]-4.0*(b[j]+a[j]/4.0)/9.0)/8.0;
                      y[j]=y[j]+3.0*hh*bb;
                  }
                  t0=t+hh/2.0;
                  func (t0,y,d);
                  for (j=0; j<=n-1; j++)
                  { 
					  c[j]=d[j];
                      qq=d[j]-15.0*(b[j]-a[j]/5.0)/16.0;
                      y[j]=y[j]+2.0*hh*qq;
                  }
                  t0=t+hh;
                  func (t0,y,d);
                  for (j=0; j<=n-1; j++)
                  { 
					  qq=c[j]-9.0*(b[j]-2.0*a[j]/9.0)/8.0;
                      qq=d[j]-8.0*qq;
                      y[j]=y[j]+hh*qq/6.0;
                  }
                  t=t+dt;
              }
              p=0.0;
              for (j=0; j<=n-1; j++)
              { 
				  qq=fabs(y[j]-v[j]);
                  if (qq>p) p=qq;
              }
              hh=hh/2.0; nn=nn+nn;
          }
          for (j=0; j<=n-1; j++) z[j][i]=y[j];
      }
	  t = aa;
  }

  void merson::output ()   //输出k个积分点上
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

  void merson::func (double t, double y[], double d[])
  {
	  double p;
	  p = 60.0*(0.06 + t*(t-0.6));
	  d[0] = p*y[1];
	  d[1] = -p*y[0];
  }

  void main ()      //主函数
  {
	  merson  s(2, 10);           //创建对象
	  s.input ();        //由文件读入数据成员t,h,eps
		                 //以及n个未知函数在起始点t处的函数值
	  s.solution ();    //执行变步长Merson方法
	  s.output (); //输出k个积分点上的未知函数值到文件并显示
  }

