  //5LAGR.CPP
  //Lagrange插值
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  lagr
  {
  private: 
           int n;
		   double  *x, *y, t, z;
  public:
	       lagr (int nn)
		   {
			   n = nn;
               x = new double[n];   //动态分配内存
	           y = new double[n];
		   }
	       void input ();          //由文件读入n个数据点(x, y)
		   double interp (double);   //执行Lagrange插值
           void output ();         //输出插值点t处的近似值z到文件
		   ~lagr ()
		   {  delete [] x, y;  }
  };

  void lagr::input ()      //由文件读入n个数据点(x, y)
  {
	  int  k;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
      for (k=0; k<n; k++)              //读入n个数据点
	  { fin >>x[k];  fin >>y[k]; }
	  fin.close ();
  }

  double lagr::interp (double tt)         //执行Lagrange插值
  { 
	  int i,j,k,m;
      double s;
	  t = tt;
      z=0.0;
      if (n<1)  { z = 0.0;  return z;  }
      if (n==1) { z = y[0];  return z;}
      if (n==2)
      { 
		  z=(y[0]*(t-x[1])-y[1]*(t-x[0]))/(x[0]-x[1]);
          return z;
      }
      i=0;
      while ((x[i]<t)&&(i<n)) i=i+1;
      k=i-4;
      if (k<0) k=0;
      m=i+3;
      if (m>n-1) m=n-1;
      for (i=k;i<=m;i++)
      { 
		  s=1.0;
          for (j=k;j<=m;j++)
            if (j!=i) s=s*(t-x[j])/(x[i]-x[j]);
          z=z+s*y[i];
      }
      return z;
  }

  void lagr::output ()       //输出插值点t处的近似值z到文件
  {
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2,ios::app);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  fout <<endl <<t <<"   " <<z <<endl;
	  fout.close ();
  }

  void main ()      //主函数
  {
	  lagr  solution(10); 
	  solution.input ();          //由文件读入n个数据点(x, y)
	  cout <<"z(0.63) = "<<solution.interp(0.63)<<endl; //执行Lagrange插值
	  solution.output ();       //输出插值点t处的近似值z到文件
  }

//z(0.63) = 0.532591
