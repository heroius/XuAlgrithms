  //5HMT.CPP
  //Hermite插值
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  hmt
  {
  private: 
           int n;
		   double  *x, *y, *dy, t, z;
  public:
	       hmt (int nn)
		   {
			   n = nn;
               x = new double[n];    //动态分配内存
	           y = new double[n];
	           dy = new double[n];
		   }
	       void input ();          //由文件读入n个数据点(x, y)
		   double interp (double);    //执行Hermite插值
           void output ();         //输出插值点t处的近似值z到文件
		   ~hmt ()
		   {  delete [] x, y, dy;  }
  };

  void hmt::input ()      //由文件读入n个数据点(x, y)
  {
	  int  k;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
      for (k=0; k<n; k++)              //读入n个数据点
	  { fin >>x[k];  fin >>y[k];  fin >>dy[k]; }
	  fin.close ();
  }

  double hmt::interp (double tt)      //执行Hermite插值
  { 
	  int i,j;
      double p,q,s;
	  t = tt;
      z=0.0;
      for (i=1;i<=n;i++)
      { 
		  s=1.0;
          for (j=1;j<=n;j++)
              if (j!=i) s=s*(t-x[j-1])/(x[i-1]-x[j-1]);
          s=s*s;
          p=0.0;
          for (j=1;j<=n;j++)
              if (j!=i) p=p+1.0/(x[i-1]-x[j-1]);
          q=y[i-1]+(t-x[i-1])*(dy[i-1]-2.0*y[i-1]*p);
          z=z+q*s;
      }
      return z;
  }

  void hmt::output ()       //输出插值点t处的近似值z到文件
  {
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2, ios::app);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  fout <<endl <<t <<"   " <<z <<endl;
	  fout.close ();
  }

  void main ()      //主函数
  {
	  hmt  solution(10); 
	  solution.input ();          //由文件读入n个数据点(x, y)
	  cout <<"z(0.356) = "<<solution.interp (0.356)<<endl;  //执行插值
	  solution.output ();          //输出插值点t处的近似值z到文件
  }

//z(0.356) = 0.70048
