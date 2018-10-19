  //6SIMP.CPP
  //变步长Simpson求积法
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  simp 
  {
  private: 
		   double  a, b, eps, integ;
  public:
	       simp (double aa, double bb, double es)//顺序提供a,b,eps值的构造函数
           { a = aa; b = bb; eps = es; }
		   void integration ();    //执行变步长Simpson求积法
           void output ();         //输出积分值到文件并显示
		   double func (double);   //计算被积函数值
  };

  void simp::integration ()   //执行变步长Simpson求积法
  { 
	  int n,k;
      double h,t1,t2,s1,s2,ep,p,x;
      n=1; h=b-a;
      t1=h*(func (a) + func (b))/2.0;
      s1=t1;
      ep=eps+1.0;
      while (ep>=eps)
      { 
		  p=0.0;
          for (k=0;k<=n-1;k++)
          { 
			  x=a+(k+0.5)*h;
              p=p+func (x);
          }
          t2=(t1+h*p)/2.0;
          s2=(4.0*t2-t1)/3.0;
          ep=fabs(s2-s1);
          t1=t2; s1=s2; n=n+n; h=h/2.0;
      }
      integ = s2;
  }

  void simp::output ()       //输出积分值到文件并显示
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

  double simp::func (double x)
  {  double y;
     y = log(1.0+x)/(1.0+x*x);
	 return y; 
  }

  void main ()      //主函数
  {
	  simp  solution(0.0, 1.0, 0.000001); //创建对象并顺序提供a, b, eps值
	  solution.integration ();    //执行变步长Simpson求积法
	  solution.output ();    //输出积分值到文件并显示
  }

