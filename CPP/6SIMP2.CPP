  //6SIMP2.CPP
  //变步长Simpson二重积分法
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  simp2
  {
  private: 
		   double  a, b, eps, integ;
  public:
	       simp2 (double aa,double bb,double es)  
			                        //顺序提供a,b,eps值的构造函数
		   { a = aa; b = bb; eps = es; }
		   void integration ();    //执行变步长Simpson二重积分法
		   double simp1(double);  //一维积分函数
           void output ();         //输出积分值到文件并显示
		   double func (double,double);   //计算被积函数值
		   double lower(double);          //计算内层积分下限值y0
		   double upper(double);          //计算内层积分上限值y1
  };

  void simp2::integration ()  //执行变步长Simpson二重积分法
  { 
      int n,j;
      double h,d,s1,s2,t1,x,t2,g,s,s0,ep;
      n=1; h=0.5*(b-a);
      d=fabs((b-a)*1.0e-06);
      s1=simp1(a); s2=simp1(b);
      t1=h*(s1+s2);
      s0=1.0e+35; ep=1.0+eps;
      while (((ep>=eps)&&(fabs(h)>d))||(n<16))
      { 
		  x=a-h; t2=0.5*t1;
          for (j=1;j<=n;j++)
          { 
			  x=x+2.0*h;
              g=simp1(x);
              t2=t2+h*g;
          }
          s=(4.0*t2-t1)/3.0;
          ep=fabs(s-s0)/(1.0+fabs(s));
          n=n+n; s0=s; t1=t2; h=h*0.5;
      }
      integ = s;
  }

  double simp2::simp1(double x)   //一维积分函数
  { 
	  int n,i;
      double y0,y1,h,d,t1,yy,t2,g,ep,g0;
      n=1;
	  y0 = lower(x);
	  y1 = upper(x);
      h=0.5*(y1-y0);
      d=fabs(h*2.0e-06);
      t1=h*(func (x,y0)+func (x,y1));
      ep=1.0+eps; g0=1.0e+35;
      while (((ep>=eps)&&(fabs(h)>d))||(n<16))
      { 
		  yy=y0-h;
          t2=0.5*t1;
          for (i=1;i<=n;i++)
          { 
			  yy=yy+2.0*h;
              t2=t2+h*func(x,yy);
          }
          g=(4.0*t2-t1)/3.0;
          ep=fabs(g-g0)/(1.0+fabs(g));
          n=n+n; g0=g; t1=t2; h=0.5*h;
      }
      return(g);
  }

  void simp2::output ()       //输出积分值到文件并显示
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

  double simp2::func (double x, double y)    //计算被积函数值
  {  return exp(x*x+y*y); }

  double simp2::lower(double x)   //计算内层积分下限值y0
  {  return (-sqrt(1-x*x));  }

  double simp2::upper(double x)   //计算内层积分上限值y1
  {  return (sqrt(1-x*x));  }

  void main ()      //主函数
  {
	  simp2  solution (0.0, 1.0, 0.0001); 
	                        //创建对象并顺序提供a,b,eps的值
	  solution.integration ();    //执行变步长Simpson二重积分法
	  solution.output ();    //输出积分值到文件并显示
  }

