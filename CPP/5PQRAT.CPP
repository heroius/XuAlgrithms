  //5PQRAT.CPP
  //连分式插值
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  pqrat
  {
  private: 
           int n;
		   double  *x, *y, t, z;
  public:
	       pqrat (int nn)
		   {
			   n = nn;
               x = new double[n];    //动态分配内存
	           y = new double[n];
		   }
	       void input ();          //由文件读入n个数据点(x, y)
		   double interp (double);   //执行连分式插值
		   void output ();         //输出插值点t处的近似值z到文件
		   ~pqrat ()
		   {  delete [] x, y;  }
  };

  void pqrat::input ()      //由文件读入n个数据点(x, y)
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

  double pqrat::interp (double tt)     //执行连分式插值
  { 
	  int i,j,k,m,l;
      double h,b[8];
	  t = tt;
      z=0.0;
      if (n<1) { z = 0.0; return z; }
      if (n==1) { z=y[0]; return z; }
      if (n<=8) { k=0; m=n;}
      else if (t<x[4]) { k=0; m=8;}
      else if (t>x[n-5]) { k=n-8; m=8; }
      else
      { 
		  k=1; j=n;
          while (j-k!=1)
          { 
			  i=(k+j)/2;
              if (t<x[i-1]) j=i;
              else k=i;
          }
          k=k-4; m=8;
      }
      b[0]=y[k];
      for (i=2;i<=m;i++)
      { 
		  h=y[i+k-1]; l=0; j=1;
          while ((l==0)&&(j<=i-1))
          { 
			  if (fabs(h-b[j-1])+1.0==1.0) l=1;
              else h=(x[i+k-1]-x[j+k-1])/(h-b[j-1]);
              j=j+1;
          }
          b[i-1]=h;
          if (l!=0) b[i-1]=1.0e+35;
      }
      z=b[m-1];
      for (i=m-1;i>=1;i--) z=b[i-1]+(t-x[i+k-1])/z;
	  return z;
  }

  void pqrat::output ()       //输出插值点t处的近似值z到文件
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
	  pqrat  solution(10); 
	  solution.input ();          //由文件读入n个数据点(x, y)
	  cout <<"z(-0.85) = "<<solution.interp(-0.85)<<endl; //执行连分式插值
	  solution.output ();     //输出插值点t处的近似值z到文件
	  cout <<"z( 0.25) = "<<solution.interp(0.25)<<endl;  //执行连分式插值
	  solution.output ();     //输出插值点t处的近似值z到文件
  }

//z(-0.85) = 0.0524591
//z( 0.25) = 0.390244
