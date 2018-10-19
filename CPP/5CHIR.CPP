  //5CHIR.CPP
  //Chebyshev曲线拟合
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  chir
  {
  private: 
           int n, m;
		   double  *x, *y, *a;
  public:
	       chir (int nn, int mm)
		   {
			   n = nn; m = mm;
               x = new double[n]; //动态分配内存
	           y = new double[n];
	           a = new double[m+2];
		   }
	       void input ();          //由文件读入n个数据点(x, y)
		   void fit ();         //执行Chebyshev曲线拟合
           void output ();         //输出m次拟合多项式的m+1个系数
		                           //以及偏差最大值到文件并显示
		   ~chir ()
		   {  delete [] x, y, a;  }
  };

  void chir::input ()      //由文件读入n个数据点(x, y)
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

  void chir::fit ()       //执行Chebyshev曲线拟合
  { 
	  int m1,i,j,l,ii,k,im,ix[21];
      double h[21],ha,hh,y1,y2,h1,h2,d,hm;
      for (i=0; i<=m+1; i++) a[i]=0.0;
      if (m>=n) m=n-2;
      if (m>=20) m=18;
      m1=m+2;
      ha=0.0;
      ix[0]=0; ix[m+1]=n-1;
      l=(n-1)/(m+1); j=l;
      for (i=1; i<=m; i++)
      { ix[i]=j; j=j+l;}
      while (1==1)
      { 
		  hh=1.0;
          for (i=0; i<=m+1; i++)
          { a[i]=y[ix[i]]; h[i]=-hh; hh=-hh;}
          for (j=1; j<=m+1; j++)
          { 
			  ii=m1; y2=a[ii-1]; h2=h[ii-1];
              for (i=j; i<=m+1; i++)
              { 
				  d=x[ix[ii-1]]-x[ix[m1-i-1]];
                  y1=a[m-i+j];
                  h1=h[m-i+j];
                  a[ii-1]=(y2-y1)/d;
                  h[ii-1]=(h2-h1)/d;
                  ii=m-i+j+1; y2=y1; h2=h1;
              }
          }
          hh=-a[m+1]/h[m+1];
          for (i=0; i<=m+1; i++)
              a[i]=a[i]+h[i]*hh;
          for (j=1; j<=m; j++)
          { 
			  ii=m-j+1; d=x[ix[ii-1]];
              y2=a[ii-1];
              for (k=m1-j; k<=m+1; k++)
              { 
				  y1=a[k-1]; a[ii-1]=y2-d*y1;
                  y2=y1; ii=k;
              }
          }
          hm=fabs(hh);
          if (hm<=ha) { a[m+1]=-hm; return;}
          a[m+1]=hm; ha=hm; im=ix[0]; h1=hh;
          j=0;
          for (i=0; i<=n-1; i++)
          { 
			  if (i==ix[j])
              { if (j<m+1) j=j+1; }
              else
              { 
				  h2=a[m];
                  for (k=m-1; k>=0; k--)
                      h2=h2*x[i]+a[k];
                  h2=h2-y[i];
                  if (fabs(h2)>hm)
                  { hm=fabs(h2); h1=h2; im=i;}
              }
          }
          if (im==ix[0]) return;
          i=0;l=1;
          while (l==1)
          { 
			  l=0;
              if (im>=ix[i])
              { 
				  i=i+1;
                  if (i<=m+1) l=1;
              }
          }
          if (i>m+1) i=m+1;
          if (i==(i/2)*2) h2=-hh;
          else h2=hh;
          if (h1*h2>=0.0) ix[i]=im;
          else
          { 
			  if (im<ix[0])
              { 
				  for (j=m; j>=0; j--)
                      ix[j+1]=ix[j];
                  ix[0]=im;
              }
              else
              { 
				  if (im>ix[m+1])
                  { 
					  for (j=1; j<=m+1; j++)
                          ix[j-1]=ix[j];
                      ix[m+1]=im;
                  }
                  else ix[i-1]=im;
              }
          }
      }
  }

  void chir::output ()       //输出m次拟合多项式的m+1个系数
		                           //以及偏差最大值到文件并显示
  {
	  int i;
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
      for (i=0; i<=m; i++)
	  {
		  fout <<a[i] <<endl;
		  cout <<a[i] <<endl;
	  }
	  fout <<endl <<a[m+1] <<endl;
	  cout <<endl <<a[m+1] <<endl;	  
	  fout.close ();
  }

  void main ()      //主函数
  {
	  int k;
	  double x, y;
	  chir  solution(101, 5); 
	  ofstream ft ("chir.txt");
	  if (!ft)
	  { cout <<"\n不能打开文件 chir.txt" <<endl; exit(1); }
	  for (k=0; k<=100; k++)
	  {
		  x = -1.0+0.02*k;  y = atan(x);
	      ft <<x <<"   " <<y <<endl;
	  }
	  ft.close ();
	  solution.input ();          //由文件读入n个数据点(x, y)
	  solution.fit ();            //执行最小二乘曲线拟合
	  solution.output ();          //输出m次拟合多项式的m+1个系数
		                           //以及偏差最大值到文件并显示
  }

