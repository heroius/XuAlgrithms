  //7EULER1.CPP
  //变步长Euler方法
  #include  <iostream>
  #include  <iomanip>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  euler1
  {
  private: 
	       int n, k;
		   double  t, h, eps, *y, **z;
		   double  *d, *a, *b, *c;
  public:
	       euler1 (int nn, int kk)
		   {
			   int i;
			   n = nn; k = kk;
               y = new double[n];     //动态分配内存
	           d = new double[n];
	           a = new double[n];
	           b = new double[n];
	           c = new double[n];
	           z = new double*[n];
	           for (i=0; i<n; i++) { z[i] = new double[k]; }
		   }
	       void  input ();         //由文件读入数据成员t,h,eps
		                      //以及n个未知函数在起始点t处的函数值
		   void solution ();    //执行变步长Euler方法
           void output ();   //输出k个积分点上
		                     //的未知函数值到文件并显示
		   void func (double,double*,double*);   
		                         //计算微分方程组中各方程右端函数值
		   ~euler1 ()
		   { 
			   int i;
			   for (i=0; i<n; i++) { delete [] z[i]; }
			   delete [] z;
			   delete [] y, d, a, b, c;
		   }
  };

  void euler1::input ()     //由文件读入数据成员t,h,eps
		                      //以及n个未知函数在起始点t处的函数值
  {
	  int  i;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
      fin >>t >>h >>eps;        //读入t,h,eps
      for (i=0; i<n; i++)  fin >>y[i];  
	                            //读入n个未知函数在起始点t处的函数值	  
	  fin.close ();
  }

  void euler1::solution ()  //执行变步长Euler方法
  { 
	  int i,j,m, kk;
      double hh,p,x,q, tt;
      tt = t;
	  for (kk=0; kk<k; kk++)
	  {
		  hh=h; m=1; p=1.0+eps;
          for (i=0; i<=n-1; i++) a[i]=y[i];
          while (p>=eps)
		  { 
		      for (i=0; i<=n-1; i++)
			  { b[i]=y[i]; y[i]=a[i]; }
              for (j=0; j<=m-1; j++)
			  { 
			      for (i=0; i<=n-1; i++) c[i]=y[i];
                  x=t+j*hh;
                  func (x,y,d);
                  for (i=0; i<=n-1; i++)
                      y[i]=c[i]+hh*d[i];
                  x=t+(j+1)*hh;
                  func (x,y,d);
                  for (i=0; i<=n-1; i++)
                      d[i]=c[i]+hh*d[i];
                  for (i=0; i<=n-1; i++)
                      y[i]=(y[i]+d[i])/2.0;
			  }
              p=0.0;
              for (i=0; i<=n-1; i++)
			  { 
			      q=fabs(y[i]-b[i]);
                  if (q>p) p=q;
			  }
              hh=hh/2.0; m=m+m;
		  }
          for (i=0; i<n; i++)  z[i][kk] = y[i];
	      t = t + h;
	  }
	  t = tt;
  }

  void euler1::output () //输出k个积分点上的未知函数值到文件并显示
  {
	  int i, j;
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
      for (i=0; i<k; i++)
	  {
		  cout <<"t = " <<t+(i+1)*h <<endl;
		  for (j=0; j<n; j++)
		  {
			  fout <<z[j][i] <<"   ";
			  cout <<"y(" <<j <<")=" <<setw(12) <<z[j][i] <<"  ";
		  }
		  fout <<endl;  cout <<endl;
	  }
	  fout.close ();
  }

  void euler1::func (double t, double y[], double d[])
  {
	  d[0] = y[1];
	  d[1] = -y[0];
	  d[2] = -y[2];
  }

  void main ()      //主函数
  {
	  euler1  s(3, 10);           //创建对象
	  s.input ();        //由文件读入数据成员t,h,eps
		                 //以及n个未知函数在起始点t处的函数值
	  s.solution ();    //执行变步长Euler方法
	  s.output ();    //输出k个积分点上的未知函数值到文件并显示
  }

