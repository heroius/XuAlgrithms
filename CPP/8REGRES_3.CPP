  //8REGRES_3.CPP
  //逐步回归分析程序
  #include  <iostream>
  #include  <iomanip>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  regres_3
  {
  private: 
	       int n, k;
		   double  **x, f1, f2, eps;
		   double  *xx, *b, *v, *s, C, F, *ye, *yr, **r;
  public:
	       regres_3 (int nn, int kk)
		   {
			   int i;
			   n = nn; k = kk;
               x = new double*[k];         //动态分配内存
	           for (i=0; i<k; i++) { x[i] = new double[n+1]; }
	           r = new double*[n+1];
	           for (i=0; i<n+1; i++) { r[i] = new double[n+1]; }
	           xx = new double[n+1];
	           b = new double[n+1];
	           v = new double[n+1];
	           s = new double[n+1];
	           ye = new double[k];
	           yr = new double[k];
		   }
	       void  input ();   //由文件读入自变量X与因变量y的k次观测值
		   void analysis (double,double,double); //执行逐步回归分析
           void output ();  //输出回归系数与各统计量到文件并显示
		   ~regres_3 ()
		   {   
			   int i;
			   for (i=0; i<k; i++) { delete [] x[i]; }
			   delete [] x;
			   for (i=0; i<n+1; i++) { delete [] r[i]; }
			   delete [] r;
			   delete [] xx, b, v, s, ye, yr; 
		   }
  };

  void regres_3::input ()  //由文件读入自变量X与因变量y的k次观测值
  {
	  int  i, j;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
      for (i=0; i<k; i++)   //读入自变量X与因变量y的k次观测值
		  for (j=0; j<=n; j++)  fin >>x[i][j];
	  fin.close ();
  }

  void regres_3::analysis (double ff1, double ff2, double es)//执行逐步回归分析
  { 
	  int i,j,ii,m,imi,imx,l,it;
      double z,phi,sd,vmi,vmx,q,fmi,fmx;
      f1 = ff1;  f2 = ff2;  eps = es;
      m=n+1; q=0.0;
      for (j=0; j<=n; j++)
      { 
		  z=0.0;
          for (i=0; i<=k-1; i++) z=z+x[i][j]/k;
          xx[j]=z;
      }
      for (i=0; i<=n; i++)
      for (j=0; j<=i; j++)
      { 
		  z=0.0;
          for (ii=0; ii<=k-1; ii++)
              z=z+(x[ii][i]-xx[i])*(x[ii][j]-xx[j]);
          r[i][j]=z;
      }
      for (i=0; i<=n; i++)  ye[i]=sqrt(r[i][i]);
      for (i=0; i<=n; i++)
      for (j=0; j<=i; j++)
      { 
		  r[i][j]=r[i][j]/(ye[i]*ye[j]);
          r[j][i]=r[i][j];
      }
      phi=k-1.0;
      sd=ye[n]/sqrt(k-1.0);
      it=1;
      while (it==1)
      { 
		  it=0;
          vmi=1.0e+35; vmx=0.0;
          imi=-1; imx=-1;
          for (i=0; i<=n; i++)
          { v[i]=0.0; b[i]=0.0; s[i]=0.0; }
          for (i=0; i<=n-1; i++)
              if (r[i][i]>=eps)
			  { 
				  v[i]=r[i][n]*r[n][i]/r[i][i];
                  if (v[i]>=0.0)
				  { 
					  if (v[i]>vmx) { vmx=v[i]; imx=i;}
				  }
                  else
				  { 
					  b[i]=r[i][n]*ye[n]/ye[i];
                      s[i]=sqrt(r[i][i])*sd/ye[i];
                      if (fabs(v[i])<vmi)
					  { vmi=fabs(v[i]); imi=i; }
				  }
			  }
          if (phi!=n-1.0)
          { 
			  z=0.0;
              for (i=0; i<=n-1; i++)  z=z+b[i]*xx[i];
              b[n]=xx[n]-z; s[n]=sd; v[n]=q;
          }
          else
          { b[n]=xx[n]; s[n]=sd; }
          fmi=vmi*phi/r[n][n];
          fmx=(phi-1.0)*vmx/(r[n][n]-vmx);
          if ((fmi<f2)||(fmx>=f1))
          { 
			  if (fmi<f2)  { phi=phi+1.0; l=imi; }
              else  { phi=phi-1.0; l=imx; }
              for (i=0; i<=n; i++)
                  if (i!=l)
                      for (j=0; j<=n; j++)
                          if (j!=l)
                              r[i][j]=r[i][j]-
                                       (r[l][j]/r[l][l])*r[i][l];
              for (j=0; j<=n; j++)
                  if (j!=l) r[l][j]=r[l][j]/r[l][l];
              for (i=0; i<=n; i++)
                  if (i!=l)  r[i][l]=-r[i][l]/r[l][l];
              r[l][l]=1.0/r[l][l];
              q=r[n][n]*ye[n]*ye[n];
              sd=sqrt(r[n][n]/phi)*ye[n];
              C=sqrt(1.0-r[n][n]);
              F=(phi*(1.0-r[n][n]))/((k-phi-1.0)*r[n][n]);
              it=1;
          }
      }
      for (i=0; i<=k-1; i++)
      { 
		  z=0.0;
          for (j=0; j<=n-1; j++) z=z+b[j]*x[i][j];
          ye[i]=b[n]+z; yr[i]=x[i][n]-ye[i];
      }
  }

  void regres_3::output () //输出回归系数与各统计量到文件并显示
  {
	  int i, j;
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2,ios::app);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  cout <<endl;
      cout <<"f1 = " <<f1 <<"    " <<"f2 = " <<f2 <<endl;
	  cout <<"eps = " <<eps <<endl <<endl;
	  fout <<endl;
	  fout <<f1 <<"    " <<f2 <<endl;
	  fout <<eps <<endl <<endl;
	  cout <<"自变量算术平均值与因变量平均值:" <<endl;
      for (i=0; i<n; i++)
	  {
		  cout <<"x(" <<i <<")=" <<xx[i] <<"    ";
		  fout <<xx[i] <<"    ";
	  }
	  cout <<"y=" <<xx[n] <<endl <<endl;
	  fout <<xx[n] <<endl <<endl;
	  cout <<"回归系数:" <<endl;
      for (i=0; i<n+1; i++)
	  {
		  cout <<"b(" <<i <<")=" <<b[i] <<"    ";
		  fout <<b[i] <<"    ";
	  }
      cout <<endl <<endl;
	  fout <<endl <<endl;
	  cout <<"偏回归平方和:" <<endl;
      for (i=0; i<n; i++)
	  {
		  cout <<"v(" <<i <<")=" <<v[i] <<"    ";
		  fout <<v[i] <<"    ";
	  }	  
	  cout <<endl <<"残差平方和q = " <<v[n] <<endl <<endl;
	  fout <<endl <<v[n] <<endl <<endl;
	  cout <<"各因子回归系数的标准偏差:" <<endl;
      for (i=0; i<n; i++)
	  {
		  cout <<"s(" <<i <<")=" <<s[i] <<"    ";
		  fout <<s[i] <<"    ";
	  }	  
	  cout <<endl <<"估计的标准偏差s = " <<s[n] <<endl <<endl;
	  fout <<endl <<s[n] <<endl <<endl;
      cout <<"复相关系数C = " <<C <<endl <<endl;
	  fout <<C <<endl <<endl;
	  cout <<"F-检验值F = " <<F <<endl <<endl;
	  fout <<F <<endl <<endl;
	  cout <<"因变量条件期望值的k个估计值以及k个观测值的残差:" <<endl;
      for (i=0; i<k; i++)
	  {
		  cout <<"ye(" <<i <<")=" <<ye[i] <<"    ";
          cout <<"yr(" <<i <<")=" <<yr[i] <<endl;
		  fout <<ye[i] <<"     " <<yr[i] <<endl;
	  }
	  cout <<endl;   fout <<endl;
	  cout <<"规格化的系数相关矩阵:" <<endl;
	  for (i=0; i<=n; i++)
	  {
		  for (j=0; j<=n; j++)
		  {
			  cout <<r[i][j] <<"    ";
			  fout <<r[i][j] <<"    ";
		  }
		  cout <<endl;  fout <<endl;
	  }

	  fout.close ();
  }

  void main ()      //主函数
  {
	  regres_3  s(4, 13);           //创建对象
	  s.input ();        //由文件读入自变量X与因变量y的k次观测值
	  s.analysis (1.46, 1.45, 1.0e-30); //执行逐步回归分析
	  s.output ();       //输出回归系数与各统计量到文件并显示

	  s.analysis (4.75, 4.67, 1.0e-30); //执行逐步回归分析
	  s.output ();       //输出回归系数与各统计量到文件并显示

	  s.analysis (9.33, 9.07, 1.0e-30); //执行逐步回归分析
	  s.output ();       //输出回归系数与各统计量到文件并显示
  }

