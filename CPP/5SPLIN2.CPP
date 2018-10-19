  //5SPLIN2.CPP
  //第二种边界条件的三次样条函数插值、微商与积分
  #include  <iostream>
  #include  <iomanip>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class splin2 
  {
  private: 
           int n, m;
		   double  *x, *y, *dy, *ddy, *t, *z, *dz, *ddz, integ;
		   double  *s;
  public:
	       splin2 (int nn, int mm)
		   {
			   n = nn; m = mm;
               x = new double[n];  //动态分配内存
	           y = new double[n];
	           dy = new double[n];
	           ddy = new double[n];
	           s = new double[n];
	           t = new double[m]; 
	           z = new double[m];
	           dz = new double[m];
	           ddz = new double[m];
		   }
	       void input ();          //由文件读入n个数据点(x, y)
		                           //以及两端点上的二阶导数值
		                           //以及m个指定插值点t
		   void interp ();    //执行三次样条函数插值、微商与积分
           void output ();         //输出n个数据点及其一阶与二阶导数值
		        //m个插值点t处的函数值及其一阶与二阶导数值与积分值到文件并显示
		   ~splin2 ()
		   {  delete [] x, y, dy, ddy, t, z, dz, ddz, s;  }
  };

  void splin2::input ()      //由文件读入n个数据点(x, y)
		                     //以及两端点上的二阶导数值
		                     //以及m个指定插值点t
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
	  fin >>ddy[0] >>ddy[n-1];        //读入两端点上的二阶导数值
      for (k=0; k<m; k++) fin >>t[k];    //读入m个插值点
	  fin.close ();
  }

  void splin2::interp ()       //执行三次样条函数插值、微商与积分
  { 
	  int i,j;
      double h0,h1,alpha,beta;
      dy[0]=-0.5;
      h0=x[1]-x[0];
      s[0]=3.0*(y[1]-y[0])/(2.0*h0)-ddy[0]*h0/4.0;
      for (j=1;j<=n-2;j++)
      { 
		  h1=x[j+1]-x[j];
          alpha=h0/(h0+h1);
          beta=(1.0-alpha)*(y[j]-y[j-1])/h0;
          beta=3.0*(beta+alpha*(y[j+1]-y[j])/h1);
          dy[j]=-alpha/(2.0+(1.0-alpha)*dy[j-1]);
          s[j]=(beta-(1.0-alpha)*s[j-1]);
          s[j]=s[j]/(2.0+(1.0-alpha)*dy[j-1]);
          h0=h1;
      }
      dy[n-1]=(3.0*(y[n-1]-y[n-2])/h1+ddy[n-1]*h1/
            2.0-s[n-2])/(2.0+dy[n-2]);
      for (j=n-2;j>=0;j--)
          dy[j]=dy[j]*dy[j+1]+s[j];
      for (j=0;j<=n-2;j++) s[j]=x[j+1]-x[j];
      for (j=0;j<=n-2;j++)
      { 
		  h1=s[j]*s[j];
          ddy[j]=6.0*(y[j+1]-y[j])/h1-2.0*(2.0*dy[j]+dy[j+1])/s[j];
      }
      h1=s[n-2]*s[n-2];
      ddy[n-1]=6.*(y[n-2]-y[n-1])/h1+2.*(2.*dy[n-1]+dy[n-2])/s[n-2];
      integ=0.0;
      for (i=0;i<=n-2;i++)
      { 
		  h1=0.5*s[i]*(y[i]+y[i+1]);
          h1=h1-s[i]*s[i]*s[i]*(ddy[i]+ddy[i+1])/24.0;
          integ=integ+h1;
      }
      for (j=0;j<=m-1;j++)
      { 
		  if (t[j]>=x[n-1]) i=n-2;
          else
          { 
			  i=0;
              while (t[j]>x[i+1]) i=i+1;
          }
          h1=(x[i+1]-t[j])/s[i];
          h0=h1*h1;
          z[j]=(3.0*h0-2.0*h0*h1)*y[i];
          z[j]=z[j]+s[i]*(h0-h0*h1)*dy[i];
          dz[j]=6.0*(h0-h1)*y[i]/s[i];
          dz[j]=dz[j]+(3.0*h0-2.0*h1)*dy[i];
          ddz[j]=(6.0-12.0*h1)*y[i]/(s[i]*s[i]);
          ddz[j]=ddz[j]+(2.0-6.0*h1)*dy[i]/s[i];
          h1=(t[j]-x[i])/s[i];
          h0=h1*h1;
          z[j]=z[j]+(3.0*h0-2.0*h0*h1)*y[i+1];
          z[j]=z[j]-s[i]*(h0-h0*h1)*dy[i+1];
          dz[j]=dz[j]-6.0*(h0-h1)*y[i+1]/s[i];
          dz[j]=dz[j]+(3.0*h0-2.0*h1)*dy[i+1];
          ddz[j]=ddz[j]+(6.0-12.0*h1)*y[i+1]/(s[i]*s[i]);
          ddz[j]=ddz[j]-(2.0-6.0*h1)*dy[i+1]/s[i];
      }
  }

  void splin2::output ()       //输出n个数据点及其一阶与二阶导数值
		 //m个插值点t处的函数值及其一阶与二阶导数值与积分值到文件并显示
  {
	  int i;
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  cout <<setw(16)<<"x[i]" <<setw(16)<<"y[i]" 
		   <<setw(16)<<"dy[i]" <<setw(16)<<"ddy[i]" <<endl;
	  for (i=0; i<n; i++)
	  {
		  fout <<setw(16)<<x[i]<<setw(16)<<y[i]
			   <<setw(16)<<dy[i]<<setw(16)<<ddy[i]<<endl;
		  cout <<setw(16)<<x[i]<<setw(16)<<y[i]
			   <<setw(16)<<dy[i]<<setw(16)<<ddy[i]<<endl;
	  }
	  fout <<endl;  cout <<endl;
	  cout <<setw(16)<<"t[i]" <<setw(16)<<"z[i]" 
		   <<setw(16)<<"dz[i]" <<setw(16)<<"ddz[i]" <<endl;
	  for (i=0; i<m; i++)
	  {
		  fout <<setw(16)<<t[i]<<setw(16)<<z[i]
			   <<setw(16)<<dz[i]<<setw(16)<<ddz[i]<<endl;
		  cout <<setw(16)<<t[i]<<setw(16)<<z[i]
			   <<setw(16)<<dz[i]<<setw(16)<<ddz[i]<<endl;
	  }
	  fout <<endl;  cout <<endl;
	  fout <<integ <<endl;  cout <<integ <<endl;
	  fout.close ();
  }

  void main ()      //主函数
  {
	  splin2  solution(12, 8); 
	  solution.input ();     //由文件读入n个数据点(x, y)
		                           //以及两端点上的二阶导数值
		                           //以及m个指定插值点t
	  solution.interp ();    //执行三次样条函数插值、微商与积分
	  solution.output ();    //输出n个数据点及其一阶与二阶导数值
		//m个插值点t处的函数值及其一阶与二阶导数值与积分值到文件并显示
  }

