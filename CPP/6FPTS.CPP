  //6FPTS.CPP
  //自适应梯形求积法
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  fpts 
  {
  private: 
		   double  a, b, eps, d, integ;
  public:
	       fpts (double aa,double bb,double es,double dd)
			                           //顺序提供a,b,eps,d值的构造函数
           { a = aa; b = bb; eps = es; d = dd; }
		   void integration ();    //执行自适应梯形求积法
		   void ppp(double,double,double,double,double,double,double);
           void output ();         //输出积分值到文件并显示
		   double func (double);   //计算被积函数值
  };

  void fpts::integration ()    //执行自适应梯形求积法
  { 
	  double h,f0,f1,t0;
      h=b-a; integ=0.0;
      f0=func (a); f1=func (b);
      t0=h*(f0+f1)/2.0;
      ppp(a,b,h,f0,f1,t0,eps);
  }

  void fpts::ppp(double x0,double x1,double h,double f0,
	             double f1,double t0,double eps)
  { 
	  double x,f,t1,t2,p,g,eps1;
      x=x0+h/2.0; f=func (x);
      t1=h*(f0+f)/4.0; t2=h*(f+f1)/4.0;
      p=fabs(t0-(t1+t2));
      if ((p<eps)||(h/2.0<d))
      { integ=integ+(t1+t2); return;}
      else
      { 
		  g=h/2.0; eps1=eps/1.4;
          ppp(x0,x,g,f0,f,t1,eps1);
          ppp(x,x1,g,f,f1,t2,eps1);
          return;
      }
  }

  void fpts::output ()       //输出积分值到文件并显示
  {
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  fout <<integ <<endl;  cout <<"积分值 = " <<integ <<endl;
	  fout.close ();
  }

  double fpts::func (double x)
  {  return 1.0/(1+25*x*x); }

  void main ()      //主函数
  {
	  fpts  solution(-1.0, 1.0, 0.000001,0.0001); 
	                         //创建对象并顺序提供a, b, eps,d值
	  solution.integration ();    //执行自适应梯形求积法
	  solution.output ();    //输出积分值到文件并显示
  }

