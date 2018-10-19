  //4CMTCL.CPP
  //Monte-Carlo法求方程一个复根
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  cmtcl
  {
  private: 
           int m;
		   double  b, eps, xr, xi;
		   double  r;
  public:
	       cmtcl(int mm,double bb, double es, double xr0, double xi0)
		   { m=mm; b=bb; eps=es;  xr=xr0; xi=xi0; } //提供b,m,eps以及根的初值xr,xi
		   void solution ();        //执行Monte-Carlo法
		   double rnd ();           //产生一个随机数
           void output ();         //输出根到文件并显示
		   double func (double, double);   //计算方程左端函数f(z)的模
  };

  void cmtcl::solution ()     //执行Monte-Carlo法
  { 
	  int k;
      double xx,yy,a,z,x1,y1,z1;
      a=b; k=1; r=1.0; xx=xr; yy=xi;
      z=func (xx,yy);
      while (a>=eps)
      { 
		  x1=-a+2.0*a*rnd(); x1=xx+x1; 
          y1=-a+2.0*a*rnd(); y1=yy+y1;
          z1=func (x1,y1);
          k=k+1;
          if (z1>=z)
          { if (k>m) { k=1; a=a/2.0; }}
          else
          { 
			  k=1; xx=x1; yy=y1;  z=z1;
              if (z<eps)
              { xr=xx; xi=yy; return; }
          }
      }
      xr=xx; xi=yy;
  }

  double cmtcl::rnd ()      //产生一个随机数
  { 
	  int m;
      double s,u,v,p;
      s=65536.0; u=2053.0; v=13849.0;
      m=(int)(r/s); r=r-m*s;
      r=u*r+v; m=(int)(r/s);
      r=r-m*s; p=r/s;
      return(p);
  }

  void cmtcl::output ()       //输出根到文件并显示
  {
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  fout <<endl <<xr <<"   " <<xi <<endl;
	  cout <<endl <<"z = " <<xr <<"  +j  " <<xi  <<endl;
	  fout.close ();
  }

//  double cmtcl::func (double x, double y)    //计算方程左端函数f(z)的模
//  {
//	  double z, u, v;
//	  u = x*x - y*y - 6*x + 13;
//	  v = 2*x*y - 6*y;
//	  z = sqrt(u*u + v*v);
//	  return z;
//  }
//z = 3 +j 2
  double cmtcl::func (double x, double y)    //计算方程左端函数f(z)的模
  {
	  double z, u, v;
	  u = x*x - y*y + x - y - 2;
	  v = 2*x*y + x + y + 2;
	  z = sqrt(u*u + v*v);
	  return z;
  }
//z = -2 +j 3.14415e-006
  void main ()      //主函数
  {
	  cmtcl  root(10, 1.0, 0.00001, 0.5, 0.5); //创建对象并提供m,b,eps,xr0,xi0
	  root.solution ();    //执行Monte-Carlo法
	  root.output ();     //输出根到文件并显示
  }

//
