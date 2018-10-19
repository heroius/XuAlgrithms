  //7RAT_PBS.CPP
  //连分式法
  #include  <iostream>
  #include  <iomanip>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  rat_pbs
  {
  private: 
	       int n, k;
		   double  t, h, eps, *y, **z;
		   double  *b, *d, *e, *u, *v, *w;
  public:
	       rat_pbs (int nn, int kk)
		   {
			   int i;
			   n = nn; k = kk;
               y = new double[n];      //动态分配内存
	           b = new double[10*n];
	           d = new double[n];
	           e = new double[n];
	           u = new double[n];
	           v = new double[n];
               w = new double[n];
	           z = new double*[n];
	           for (i=0; i<n; i++) { z[i] = new double[k]; }
		   }
	       void  input ();         //由文件读入数据成员t,h,eps
		                      //以及n个未知函数在起始点t处的函数值
		   void solution ();    //执行连分式法
		   void rkt(double,double);  //执行Runge-Kutta法
           void output ();        //输出k个积分点上
		                          //的未知函数值到文件并显示
		   void func (double,double*,double*);   
		                         //计算微分方程组中各方程右端函数值
		   ~rat_pbs ()
		   { 
			   int i;
			   for (i=0; i<n; i++) { delete [] z[i]; }
			   delete [] z;
			   delete [] y, b, d, e, u, v, w;
		   }
  };

  void rat_pbs::input ()     //由文件读入数据成员t,h,eps
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

  void rat_pbs::solution ()  //执行连分式法
  { 
      int i,j,jj,kk,m,nn,it;
      double x,hh,dd,q,p,g[10],tt;
      tt = t;
	  for (jj=0; jj<k; jj++)
	  {
		  t = tt + jj*h;
          for (j=0; j<=n-1; j++) v[j]=y[j];
          x=t; nn=1; hh=h; g[0]=hh;
          rkt(x,hh);
          for (j=0; j<=n-1; j++)
		  { b[j]=y[j]; u[j]=y[j]; }
          kk=1; it=1;
          while (it==1)
		  { 
			  nn=nn+nn; hh=hh/2.0; it=0;
              g[kk]=hh;
              for (j=0; j<=n-1; j++) y[j]=v[j];
              t=x;
              for (j=0; j<=nn-1; j++)
			  { 
				  rkt(t,hh);
                  t=t+hh;
			  }
              for (j=0; j<=n-1; j++)
			  { 
				  dd=y[j]; m=0;
                  for (i=0; i<=kk-1; i++)
                  if (m==0)
				  { 
					  q=dd-b[i*n+j];
                      if (fabs(q)+1.0==1.0) m=1;
                      else dd=(g[kk]-g[i])/q;
				  }
                  b[kk*n+j]=dd;
                  if (m!=0) b[kk*n+j]=1.0e+35;
			  }
              for (j=0; j<=n-1; j++)
			  { 
				  dd=0.0;
                  for (i=kk-1; i>=0; i--)
                      dd=-g[i]/(b[(i+1)*n+j]+dd);
                  y[j]=dd+b[j];
			  }
              p=0.0;
              for (j=0; j<=n-1; j++)
			  { 
				  q=fabs(y[j]-u[j]);
                  if (q>p) p=q;
			  }
              if ((p>=eps)&&(kk<7))
			  { 
				  for (j=0; j<=n-1; j++) u[j]=y[j];
                  kk=kk+1; it=1;
			  }
		  }
	      for (j=0; j<n; j++)  z[j][jj] = y[j];
	  }
      t = tt;
  }

  void rat_pbs::rkt(double t, double h)
  { 
	  int i,kk;
      double a[4],tt;
      a[0]=h/2.0; a[1]=a[0]; a[2]=h; a[3]=h;
      func (t,y,d);
      for (i=0; i<=n-1; i++) { w[i]=y[i]; e[i]=y[i];}
      for (kk=0; kk<=2; kk++)
      { 
		  for (i=0; i<=n-1; i++)
          { 
			  y[i]=e[i]+a[kk]*d[i];
              w[i]=w[i]+a[kk+1]*d[i]/3.0;
          }
          tt=t+a[kk];
          func (tt,y,d);
      }
      for (i=0; i<=n-1; i++)  y[i]=w[i]+h*d[i]/6.0;
  }

  void rat_pbs::output ()   //输出k个积分点上
	                      //的未知函数值到文件并显示
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

  void rat_pbs::func (double t, double y[], double d[])
  {
	  d[0] = -y[1];
	  d[1] = y[0];
  }

  void main ()      //主函数
  {
	  rat_pbs  s(2, 10);           //创建对象
	  s.input ();        //由文件读入数据成员t,h,eps
		                 //以及n个未知函数在起始点t处的函数值
	  s.solution ();    //执行连分式法
	  s.output (); //输出k个积分点上的未知函数值到文件并显示
  }

