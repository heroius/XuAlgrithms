  //8REGRES.CPP
  //一元线性回归分析程序
  #include  <iostream>
  #include  <iomanip>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  regres
  {
  private: 
	       int n;
		   double  *x, *y, a, b;
		   double  q, s, p, umax, umin, u;
  public:
	       regres (int nn)
		   {
			   n = nn;
               x = new double[n];   //动态分配内存
	           y = new double[n];
		   }
	       void  input ();      //由文件读入n个观测点(x, y)
		   void analysis ();    //执行一元线性回归分析
           void output ();  //输出回归系数a，b以及各统计量到文件并显示
		   ~regres ()
		   {   delete [] x, y; }
  };

  void regres::input ()     //由文件读入n以及n个观测点(x, y)
  {
	  int  i;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
      for (i=0; i<n; i++)  fin >>x[i];  
      for (i=0; i<n; i++)  fin >>y[i];  
	  fin.close ();
  }

  void regres::analysis()   //执行一元线性回归分析
  { 
	  int i;
      double xx,yy,e,f;
      xx=0.0; yy=0.0;
      for (i=0; i<=n-1; i++)
	  { xx=xx+x[i]/n; yy=yy+y[i]/n; }
      e=0.0; f=0.0;
      for (i=0; i<=n-1; i++)
      { 
		  q=x[i]-xx; e=e+q*q;
          f=f+q*(y[i]-yy);
      }
      a = f/e;  b = yy-a*xx;
      q=0.0; u=0.0; p=0.0;
      umax=0.0; umin=1.0e+30;
      for (i=0; i<=n-1; i++)
      { 
		  s=a*x[i]+b;
          q=q+(y[i]-s)*(y[i]-s);
          p=p+(s-yy)*(s-yy);
          e=fabs(y[i]-s);
          if (e>umax) umax=e;
          if (e<umin) umin=e;
          u=u+e/n;
      }
    s = sqrt(q/n);
  }

  void regres::output () //输出回归系数a，b以及各统计量到文件并显示
  {
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
      cout <<"a = " <<a <<"     " <<"b = " <<b <<endl;
	  fout <<a <<"     " <<b <<endl;
	  cout <<"偏差平方和 = " <<q <<endl;    fout <<q <<endl;
	  cout <<"平均标准偏差 = " <<s <<endl;  fout <<s <<endl;
	  cout <<"回归平方和 = " <<p <<endl;    fout <<p <<endl;
	  cout <<"最大偏差 = " <<umax <<endl;   fout <<umax <<endl;
	  cout <<"最小偏差 = " <<umin <<endl;   fout <<umin <<endl;
	  cout <<"偏差平均值 = " <<u <<endl;    fout <<u <<endl;
	  fout.close ();
  }

  void main ()      //主函数
  {
	  regres  s(11);           //创建对象
	  s.input ();        //由文件读入n个观测点(x, y)
	  s.analysis ();     //执行一元线性回归分析
	  s.output (); //输出回归系数a，b以及各统计量到文件并显示
  }

