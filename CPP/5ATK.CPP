  //5ATK.CPP
  //Aitken插值
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  atk
  {
  private: 
           int n;
		   double  *x, *y, t, z;
  public:
	       atk (int nn)
		   {
			   n = nn;
               x = new double[n];   //动态分配内存
	           y = new double[n];
		   }
	       void input ();          //由文件读入n个数据点(x, y)
		   double interp (double, double);    //执行Aitken插值
           void output ();         //输出插值点t处的近似值z到文件
		   ~atk ()
		   {  delete [] x, y;  }
  };

  void atk::input ()      //由文件读入n个数据点(x, y)
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

  double atk::interp (double tt,double eps)   //执行Aitken插值
  { 
	  int i,j,k,m,l;
      double xx[10],yy[10];
	  t = tt;
      z=0.0;
      if (n<1) return z ;
      if (n==1) { z=y[0]; return z ;}
      m=10;
      if (m>n) m=n;
      if (t<=x[0]) k=1;
      else if (t>=x[n-1]) k=n;
      else
      { 
		  k=1; j=n;
          while ((k-j!=1)&&(k-j!=-1))
          { 
			  l=(k+j)/2;
              if (t<x[l-1]) j=l;
              else k=l;
          }
          if (fabs(t-x[l-1])>fabs(t-x[j-1])) k=j;
      }
      j=1; l=0;
      for (i=1;i<=m;i++)
      { 
		  k=k+j*l;
          if ((k<1)||(k>n))
          { l=l+1; j=-j; k=k+j*l;}
          xx[i-1]=x[k-1]; yy[i-1]=y[k-1];
          l=l+1; j=-j;
      }
      i=0;
      do
      { 
		  i=i+1; z=yy[i];
          for (j=0;j<=i-1;j++)
              z=yy[j]+(t-xx[j])*(yy[j]-z)/(xx[j]-xx[i]);
          yy[i]=z;
      }
      while ((i!=m-1)&&(fabs(yy[i]-yy[i-1])>eps));
      return z;
  }

  void atk::output ()       //输出插值点t处的近似值z到文件
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
	  atk  solution(11); 
	  solution.input ();          //由文件读入n个数据点(x, y)
	  cout <<"z(-0.75) = "<<solution.interp (-0.75,0.000001)<<endl;  //执行Aitken插值
	  solution.output ();          //输出插值点t处的近似值z到文件
	  cout <<"z( 0.05) = "<<solution.interp (0.05,0.000001)<<endl;  //执行Aitken插值
	  solution.output ();          //输出插值点t处的近似值z到文件
  }

//z(-0.75) = -0.00308891
//z( 0.05) = 0.959859
