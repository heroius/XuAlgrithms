  //5LAGR2.CPP
  //二元Lagrange插值
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  lagr2
  {
  private: 
           int n, m;
		   double  *x, *y, **z, u, v, f;
  public:
	       lagr2 (int nn, int mm)
		   {
			   int k;
			   n = nn;  m = mm;
               x = new double[n]; //动态分配内存
	           y = new double[m];
	           z = new double*[n];
	           for (k=0; k<n; k++) { z[k] = new double[m]; }
		   }
	       void input ();          //由文件读入n*m个数据点(x, y, z)
		   double interp (double, double);    //执行二元Lagrange插值
           void output ();         //输出插值点u,v处的近似值f到文件
		   ~lagr2 ()
		   {  
			   int i;
			   for (i=0; i<n; i++) { delete [] z[i]; }
			   delete [] z;
			   delete [] x, y;  
		   }
  };

  void lagr2::input ()      //由文件读入n*m个数据点(x, y, z)
  {
	  int  k, j;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
      for (k=0; k<n; k++) { fin >>x[k]; }
	  for (k=0; k<m; k++) { fin >>y[k]; }
	  for (k=0; k<n; k++)                   //读入n*m个数据点
	  for (j=0; j<m; j++)  { fin >>z[k][j]; }
	  fin.close ();
  }

  double lagr2::interp(double uu,double vv)
  { 
	  int ip,ipp,i,j,kk,iq,iqq,k;
      double h,b[10];
	  u = uu;  v = vv;
      if (u<=x[0]) { ip=1; ipp=4;}
      else if (u>=x[n-1]) { ip=n-3; ipp=n;}
      else
      { 
		  i=1; j=n;
          while (((i-j)!=1)&&((i-j)!=-1))
          { 
			  kk=(i+j)/2;
              if (u<x[kk-1]) j=kk;
              else i=kk;
          }
          ip=i-3; ipp=i+4;
      }
      if (ip<1) ip=1;
      if (ipp>n) ipp=n;
      if (v<=y[0]) { iq=1; iqq=4;}
      else if (v>=y[m-1]) { iq=m-3; iqq=m;}
      else
      { 
		  i=1; j=m;
          while (((i-j)!=1)&&((i-j)!=-1))
          { 
			  kk=(i+j)/2;
              if (v<y[kk-1]) j=kk;
              else i=kk;
          }
          iq=i-3; iqq=i+4;
      }
      if (iq<1) iq=1;
      if (iqq>m) iqq=m;
      for (i=ip-1;i<=ipp-1;i++)
      { 
		  b[i-ip+1]=0.0;
          for (j=iq-1;j<=iqq-1;j++)
          { 
			  h=z[i][j];
              for (k=iq-1;k<=iqq-1;k++)
                  if (k!=j) h=h*(v-y[k])/(y[j]-y[k]);
              b[i-ip+1]=b[i-ip+1]+h;
          }
      }
      f=0.0;
      for (i=ip-1;i<=ipp-1;i++)
      { 
		  h=b[i-ip+1];
          for (j=ip-1;j<=ipp-1;j++)
              if (j!=i) h=h*(u-x[j])/(x[i]-x[j]);
          f=f+h;
      }
      return f;
  }

  void lagr2::output ()       //输出插值点u,v处的近似值f到文件
  {
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2, ios::app);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  fout <<endl <<u <<"  " <<v <<"     " <<f <<endl;
	  fout.close ();
  }

  void main ()      //主函数
  {
	  int i, j;
	  lagr2  solution(11, 11); 
      ofstream ft ("lagr2.txt");
	  if (!ft)
	  { cout <<"\n不能打开文件lagr2.txt " <<endl; exit(1); }
      for (i=0; i<11; i++)  ft <<0.1*i <<"  ";
	  ft <<endl;
	  for (i=0; i<11; i++)  ft <<0.1*i <<"  ";
	  ft <<endl;
	  for (i=0; i<11; i++)
	  {
		  for (j=0; j<11; j++)
			  ft <<exp(-(0.1*i-0.1*j)) <<" ";
		  ft <<endl;
	  }
	  ft <<endl;
	  ft.close ();
	  solution.input ();    //由文件读入n*m个数据点(x, y, z)
	  cout <<"x = 0.350,  y = 0.650,  z(x, y) = ";
	  cout <<solution.interp (0.350, 0.650)<<endl;  //执行二元Lagrange插值
	  solution.output ();          //输出插值点u,v处的近似值f到文件
	  cout <<"x = 0.450,  y = 0.550,  z(x, y) = ";
	  cout <<solution.interp (0.450, 0.550)<<endl;  //执行二元Lagrange插值
	  solution.output ();          //输出插值点u,v处的近似值f到文件
  }

