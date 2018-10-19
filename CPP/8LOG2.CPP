  //8LOG2.CPP
  //对数数据相关程序
  #include  <iostream>
  #include  <iomanip>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  log2
  {
  private: 
	       int n;
		   double  *x, *y;
		   double  a, b, q, s, umax, umin, u;
  public:
	       log2 (int nn)
		   {
			   n = nn;
               x = new double[n];    //动态分配内存
	           y = new double[n];
		   }
	       void  input ();      //由文件读入n个数据点(x, y)
		   void solution ();    //执行对数数据相关
           void output ();  //输出拟合函数参数a,b以及各统计量到文件并显示
		   ~log2 ()
		   {   delete [] x, y; }
  };

  void log2::input ()     //由文件读入n个数据点(x, y)
  {
	  int  i;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
      for (i=0; i<n; i++)  fin >>x[i];    //读入n个自变量值
      for (i=0; i<n; i++)  fin >>y[i];    //读入n个因变量值
	  fin.close ();
  }

  void log2::solution ()    //执行对数数据相关
  { 
	  int i;
      double xx,yy,dx,dxy;
      xx=0.0; yy=0.0;
      for (i=0; i<=n-1; i++)
      { 
		  xx=xx+log(x[i])/n; 
          yy=yy+log(y[i])/n;
      }
      dx=0.0; dxy=0.0;
      for (i=0; i<=n-1; i++)
      { 
		  q=log(x[i])-xx; dx=dx+q*q;
          dxy=dxy+q*(log(y[i])-yy);
      }
      a=dxy/dx; b=yy-a*xx;
      b=exp(b);
      q=0.0; u=0.0; umax=0.0; umin=1.0e+30;
      for (i=0; i<=n-1; i++)
      { 
		  s=a*log(x[i]); 
		  s=b*exp(s);
          q=q+(y[i]-s)*(y[i]-s);
          dx=fabs(y[i]-s);
          if (dx>umax) umax=dx;
          if (dx<umin) umin=dx;
          u=u+dx/n;
      }
      s=sqrt(q/n);
  }

  void log2::output ()     //输出拟合函数参数a,b以及各统计量到文件并显示
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
	  cout <<"最大偏差 = " <<umax <<endl;   fout <<umax <<endl;
	  cout <<"最小偏差 = " <<umin <<endl;   fout <<umin <<endl;
	  cout <<"偏差平均值 = " <<u <<endl;    fout <<u <<endl;
	  fout.close ();
  }

  void main ()      //主函数
  {
	  log2  s(10);           //创建对象
	  s.input ();        //由文件读入n个数据点(x, y)
	  s.solution ();     //执行对数数据相关
	  s.output (); //输出拟合函数参数a,b以及各统计量到文件并显示
  }

