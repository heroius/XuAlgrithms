  //5AKIMA.CPP
  //Akima方法
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  akima
  {
  private: 
           int n, k;
		   double  *x, *y, t, z, s[4];
  public:
	       akima (int nn)
		   {
			   n = nn;
               x = new double[n];    //动态分配内存
	           y = new double[n];
		   }
	       void input ();          //由文件读入n个数据点(x, y)
		   void interp (double);  //计算插值点t所在子区间上的三次多项式
		                            //并计算t点的近似值z
           void output (); //输出三次多项式系数以及t点的近似值z到文件并显示
		   ~akima ()
		   {  delete [] x, y;  }
  };

  void akima::input ()      //由文件读入n个数据点(x, y)
  {
	  int  j;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
      for (j=0; j<n; j++)              //读入n个数据点
	  { fin >>x[j];  fin >>y[j]; }
	  fin.close ();
  }

  void akima::interp(double tt)    ////计算插值点t所在子区间上的三次多项式
		                            //并计算t点的近似值z
  { 
	  int m,kk;
      double u[5],p,q;
	  t = tt;
      z=0.0; s[0]=0.0; s[1]=0.0; s[2]=0.0; s[3]=0.0;
      if (n<1) { k = 0; return; }
      if (n==1) { k = 0; s[0]=y[0]; z=y[0]; return;}
      if (n==2)
      { 
		  k = 0;
		  s[0]=y[0]; s[1]=(y[1]-y[0])/(x[1]-x[0]);
          z=(y[0]*(t-x[1])-y[1]*(t-x[0]))/(x[0]-x[1]);
          return;
      }
      if (t<=x[1]) k=0;
      else if (t>=x[n-1]) k=n-2;
      else
      { 
		  k=1; m=n;
          while (((k-m)!=1)&&((k-m)!=-1))
          { 
			  kk=(k+m)/2;
              if (t<x[kk-1]) m=kk;
              else k=kk;
          }
          k=k-1;
      }
      u[2]=(y[k+1]-y[k])/(x[k+1]-x[k]);
      if (n==3)
      { 
		  if (k==0)
          { 
			  u[3]=(y[2]-y[1])/(x[2]-x[1]);
              u[4]=2.0*u[3]-u[2];
              u[1]=2.0*u[2]-u[3];
              u[0]=2.0*u[1]-u[2];
          }
          else
          { 
			  u[1]=(y[1]-y[0])/(x[1]-x[0]);
              u[0]=2.0*u[1]-u[2];
              u[3]=2.0*u[2]-u[1];
              u[4]=2.0*u[3]-u[2];
          }
      }
      else
      { 
		  if (k<=1)
          { 
			  u[3]=(y[k+2]-y[k+1])/(x[k+2]-x[k+1]);
              if (k==1)
              { 
				  u[1]=(y[1]-y[0])/(x[1]-x[0]);
                  u[0]=2.0*u[1]-u[2];
                  if (n==4) u[4]=2.0*u[3]-u[2];
                  else u[4]=(y[4]-y[3])/(x[4]-x[3]);
              }
              else
              { 
				  u[1]=2.0*u[2]-u[3];
                  u[0]=2.0*u[1]-u[2];
                  u[4]=(y[3]-y[2])/(x[3]-x[2]);
              }
          }
          else if (k>=(n-3))
          { 
			  u[1]=(y[k]-y[k-1])/(x[k]-x[k-1]);
              if (k==(n-3))
              { 
				  u[3]=(y[n-1]-y[n-2])/(x[n-1]-x[n-2]);
                  u[4]=2.0*u[3]-u[2];
                  if (n==4) u[0]=2.0*u[1]-u[2];
                  else u[0]=(y[k-1]-y[k-2])/(x[k-1]-x[k-2]);
              }
              else
              { 
				  u[3]=2.0*u[2]-u[1];
                  u[4]=2.0*u[3]-u[2];
                  u[0]=(y[k-1]-y[k-2])/(x[k-1]-x[k-2]);
              }
          }
          else
          { 
			  u[1]=(y[k]-y[k-1])/(x[k]-x[k-1]);
              u[0]=(y[k-1]-y[k-2])/(x[k-1]-x[k-2]);
              u[3]=(y[k+2]-y[k+1])/(x[k+2]-x[k+1]);
              u[4]=(y[k+3]-y[k+2])/(x[k+3]-x[k+2]);
          }
      }
      s[0]=fabs(u[3]-u[2]);
      s[1]=fabs(u[0]-u[1]);
      if ((s[0]+1.0==1.0)&&(s[1]+1.0==1.0))
           p=(u[1]+u[2])/2.0;
      else p=(s[0]*u[1]+s[1]*u[2])/(s[0]+s[1]);
      s[0]=fabs(u[3]-u[4]);
      s[1]=fabs(u[2]-u[1]);
      if ((s[0]+1.0==1.0)&&(s[1]+1.0==1.0))
         q=(u[2]+u[3])/2.0;
      else q=(s[0]*u[2]+s[1]*u[3])/(s[0]+s[1]);
      s[0]=y[k];
      s[1]=p;
      s[3]=x[k+1]-x[k];
      s[2]=(3.0*u[2]-2.0*p-q)/s[3];
      s[3]=(q+p-2.0*u[2])/(s[3]*s[3]);
      p=t-x[k];
      z=s[0]+s[1]*p+s[2]*p*p+s[3]*p*p*p;
  }

  void akima::output ()//输出三次多项式系数以及t点的近似值z到文件并显示
  {
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2, ios::app);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
      fout <<endl;  cout <<endl;
	  fout <<k <<":" <<endl;
	  fout <<s[0] <<"   " <<s[1] <<"   " <<s[2] <<"   " <<s[3] <<endl;
	  cout <<k <<":" <<endl;
	  cout <<s[0] <<"   " <<s[1] <<"   " <<s[2] <<"   " <<s[3] <<endl;
	  fout <<endl <<t <<"   " <<z <<endl;
	  cout <<endl <<t <<"   " <<z <<endl;
	  fout.close ();
  }

  void main ()      //主函数
  {
	  akima  solution(11); 
	  solution.input ();          //由文件读入n个数据点(x, y)
	  solution.interp (-0.85);         //执行Akima方法
	  solution.output ();//输出三次多项式系数以及t点的近似值z到文件并显示
	  solution.interp (0.15);         //执行Akima方法
	  solution.output ();//输出三次多项式系数以及t点的近似值z到文件并显示	  
  }

