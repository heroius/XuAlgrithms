  //4RMTCL.CPP
  //Monte-Carlo法求方程一个实根
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  rmtcl
  {
  private: 
           int m;
		   double  b, eps, x;
		   double  r;
  public:
	       rmtcl(int mm,double bb, double es, double x0)
		   { m=mm; b=bb; eps=es;  x=x0; }    //提供b,m,eps以及根的初值x
		   void solution ();        //执行Monte-Carlo法
		   double rnd ();           //产生一个随机数
           void output ();         //输出根到文件并显示
		   double func (double);   //计算方程左端函数f(x)值
  };

  void rmtcl::solution ()       //执行Monte-Carlo法
  { 
	  int k;
      double xx,a,y,x1,y1;
      a=b; k=1; r=1.0; xx=x; y=func (xx);
      while (a>=eps)
      { x1=rnd(); x1=-a+2.0*a*x1;
        x1=xx+x1; y1=func (x1);
        k=k+1;
        if (fabs(y1)>=fabs(y))
          { if (k>m) { k=1; a=a/2.0; }}
        else
          { k=1; xx=x1; y=y1;
            if (fabs(y)<eps)
              { x=xx; return; }
          }
      }
      x=xx; 
  }

  double rmtcl::rnd ()      //产生一个随机数
  { 
	  int m;
      double s,u,v,p;
      s=65536.0; u=2053.0; v=13849.0;
      m=(int)(r/s); r=r-m*s;
      r=u*r+v; m=(int)(r/s);
      r=r-m*s; p=r/s;
      return(p);
  }

  void rmtcl::output ()       //输出根到文件并显示
  {
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  fout <<endl <<x <<endl;
	  cout <<endl <<"x = " <<x <<endl;
	  fout.close ();
  }

  double rmtcl::func (double x)    //计算方程左端函数f(x)值
  {
	  double y;
	  y = exp(-x*x*x) - sin(x)/cos(x) + 800;
	  return y;
  }

  void main ()      //主函数
  {
	  rmtcl  root(10, 1.0, 0.00001, 0.5); //创建对象并提供m,b,eps,x0
	  root.solution ();    //执行Monte-Carlo法
	  root.output ();     //输出根到文件并显示
  }

//x = 1.56955
