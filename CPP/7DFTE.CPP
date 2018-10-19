  //7DFTE.CPP
  //求解二阶微分方程边值问题
  #include  <iostream>
  #include  <iomanip>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  dfte
  {
  private: 
	       int n;
		   double  a, b, ya, yb, *y, u, v, w, f;
		   double  *g, *d;
  public:
	       dfte (int nn)
		   {
			   n = nn;
               y = new double[n];     //动态分配内存
	           d = new double[2*n];
	           g = new double[6*n];	
		   }
	       void  input ();       //由文件读入数据成员a,b,ya,yb
		   void solution ();    //执行有限差分法
		   void trde (int,double*,double*); //执行追赶法
           void output ();        //输出n个等距离散点上
		                          //的未知函数值到文件并显示
		   void func (double);   
		                  //计算二阶微分方程中的函数u(x),v(x),w(x),f(x)
		   ~dfte ()
		   { delete [] y, g, d;  }
  };

  void dfte::input ()     //由文件读入数据成员a,b,ya,yb
  {
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
      fin >>a >>b >>ya >>yb;        //读入a,b,ya,yb  
	  fin.close ();
  }

  void dfte::solution ()   //执行有限差分法
  { 
	  int i,k,nn;
      double h,x;
      h=(b-a)/(n-1.0); nn=2*n-1;
      g[0]=1.0; g[1]=0.0;
      y[0]=ya; y[n-1]=yb;
      g[3*n-3]=1.0; g[3*n-4]=0.0;
      for (i=2; i<=n-1; i++)
      { 
		  x=a+(i-1)*h;
          func (x);
          k=3*(i-1)-1;
          g[k]=u-h*v/2.0;
          g[k+1]=h*h*w-2.0*u;
          g[k+2]=u+h*v/2.0;
          y[i-1]=h*h*f;
      }
      trde(n,g,y);
      h=h/2.0;
      g[0]=1.0; g[1]=0.0;
      d[0]=ya; d[nn-1]=yb;
      g[3*nn-3]=1.0; g[3*nn-4]=0.0;
      for (i=2; i<=nn-1; i++)
      { 
		  x=a+(i-1)*h;
          func (x);
          k=3*(i-1)-1;
          g[k]=u-h*v/2.0;
          g[k+1]=h*h*w-2.0*u;
          g[k+2]=u+h*v/2.0;
          d[i-1]=h*h*f;
      }
      trde(nn,g,d);
      for (i=2; i<=n-1; i++)
      { 
		  k=i+i-1;
          y[i-1]=(4.0*d[k-1]-y[i-1])/3.0;
      }
  }

  void dfte::trde (int n, double b[], double d[])  //执行追赶法
  { 
	  int k,j,m;
      double s;
      m= 3*n-2;
      for (k=0;k<=n-2;k++)
      { 
		  j=3*k; s=b[j];
          if (fabs(s)+1.0==1.0)
          { 
			  cout <<"\n程序工作失败！无解. " <<endl;
			  return;
		  }
          b[j+1]=b[j+1]/s;
          d[k]=d[k]/s;
          b[j+3]=b[j+3]-b[j+2]*b[j+1];
          d[k+1]=d[k+1]-b[j+2]*d[k];
      }
      s=b[3*n-3];
      if (fabs(s)+1.0==1.0)
      { 
		  cout <<"\n程序工作失败！无解. " <<endl;
		  return;
	  }
      d[n-1]=d[n-1]/s;
      for (k=n-2;k>=0;k--)
          d[k]=d[k]-b[3*k+1]*d[k+1];
  }


  void dfte::output ()    //输出n个等距离散点上
		                          //的未知函数值到文件并显示
  {
	  int i;
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
      for (i=0; i<n; i++)
	  {
		  cout <<"x = " <<a+i*(b-a)/(n-1) <<"     " <<"y = " <<y[i] <<endl;
		  fout <<a+i*(b-a)/(n-1) <<"     " <<y[i] <<endl;
	  }
	  fout.close ();
  }

  void dfte::func (double x)
  {
	  u = -1.0;
	  v = 0.0;
	  w = 2.0/(x*x);
	  f = 1.0/x;
  }

  void main ()      //主函数
  {
	  dfte  s(11);           //创建对象
	  s.input ();        //由文件读入数据成员a,b,ya,yb
	  s.solution ();    //执行有限差分法
	  s.output ();      //输出n个等距离散点上
		                          //的未知函数值到文件并显示
  }

