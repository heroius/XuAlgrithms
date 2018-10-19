  //6FFTS.CPP
  //变步长梯形求积法
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  ffts 
  {
  private: 
		   double  a, b, eps, integ;
  public:
	       ffts (double aa, double bb, double es)//顺序提供a,b,eps值的构造函数
           { a = aa; b = bb; eps = es; }
		   void integration ();    //执行变步长梯形求积法
           void output ();         //输出积分值到文件并显示
		   double func (double);   //计算被积函数值
  };

  void ffts::integration ()       //执行变步长梯形求积法
  { 
	  int n,k;
      double fa,fb,h,t1,p,s,x,t;
      fa=func (a);  fb=func (b);
      n=1; h=b-a;
      t1=h*(fa+fb)/2.0;
      p=eps+1.0;
      while (p>=eps)
      { 
		  s=0.0;
          for (k=0;k<=n-1;k++)
          { 
			  x=a+(k+0.5)*h;
              s=s+func (x);
          }
          t=(t1+h*s)/2.0;
          p=fabs(t1-t);
          t1=t; n=n+n; h=h/2.0;
      }
      integ = t;
  }

  void ffts::output ()       //输出积分值到文件并显示
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

  double ffts::func (double x)
  {  return exp(-x*x); }

  void main ()      //主函数
  {
	  ffts  solution(0.0, 1.0, 0.000001); //创建对象并顺序提供a, b, eps值
	  solution.integration ();    //执行变步长梯形求积法
	  solution.output ();    //输出积分值到文件并显示
  }

