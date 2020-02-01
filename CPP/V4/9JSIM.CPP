  //9JSIM.CPP
  //单形调优法求n维极值
  #include  <iostream>
  #include  <iomanip>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  jsim
  {
  private: 
	       int n, k, flag;
		   double  d, u, v, eps, *x, **xx, *f;
		   double  *xt, *xe, *xf;
  public:
	       jsim (int nn)
		   {
			   int i;
			   n = nn;
               x = new double[n+1];      //动态分配内存
	           f = new double[n+1];
	           xt = new double[n];
	           xe = new double[n];
	           xf = new double[n];
	           xx = new double*[n];
	           for (i=0; i<n; i++) { xx[i] = new double[n+1]; }
		   }
	       void  input ();      //由文件读入k,d,u,v,eps
		   void solution ();    //执行单形调优法
           void output ();     //输出实际迭代次数以及极小值点
		                       //与极小值到文件并显示
		   double func (double*); //计算目标函数值
		   ~jsim ()
		   {
			   int i;
			   for (i=0; i<n; i++) { delete [] xx[i]; }
			   delete [] xx;
			   delete [] x, f, xt, xe, xf; 
		   }
  };

  void jsim::input ()     //由文件读入k,d,u,v,eps
  {
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
      fin >>k >>d >>u >>v >>eps;     //读入k,d,u,v,eps
	  fin.close ();
  }

  void jsim::solution ()    //执行单形调优法
  { 
	  int r,g,i,j,l,kk;
      double nn,fe,fr,fl,fg,ft,ff;
      kk=0; nn=1.0*n;
      fr=sqrt(nn+1.0);
      fl=d*(fr-1.0)/(1.414*nn);
      fg=d*(fr+nn-1.0)/(1.414*nn);
      for(i=0; i<=n-1; i++)
      for(j=0; j<=n; j++)  xx[i][j]=0.0;
      for(i=1; i<=n; i++)
      for(j=0; j<=n-1; j++) xx[j][i]=fl;
      for(i=1; i<=n; i++)  xx[i-1][i]=fg;
      for(i=0; i<=n; i++)
      { 
		  for(j=0; j<=n-1; j++)   xt[j]=xx[j][i]; 
          f[i]=func (xt);
      }
      ft=1.0+eps;
      while ((kk<k)&&(ft>eps))
      { 
		  kk=kk+1;
          fr=f[0]; fl=f[0]; r=0; l=0;
          for (i=1; i<=n; i++)
          { 
			  if (f[i]>fr) { r=i; fr=f[i];}
              if (f[i]<fl) { l=i; fl=f[i];}
          }
          g=0; fg=f[0];
          j=0;
          if (r==0) { g=1; fg=f[1]; j=1;}
          for (i=j+1; i<=n; i++)
              if ((i!=r)&&(f[i]>fg)) { g=i; fg=f[i];}
          for (j=0; j<=n-1; j++)
          { 
			  xf[j]=0.0;
              for (i=0; i<=n; i++)
                  if (i!=r) xf[j]=xf[j]+xx[j][i]/nn;
              xt[j]=2.0*xf[j]-xx[j][r];
          }
          ft=func (xt);
          if (ft<f[l])
          { 
			  for (j=0; j<=n-1; j++)  xf[j]=(1.0+u)*xt[j]-u*xf[j];
              ff=func (xf);
              if (ff<f[l])
              { 
				  for (j=0; j<=n-1; j++) xx[j][r]=xf[j];
                  f[r]=ff;
              }
              else
              { 
				  for (j=0; j<=n-1; j++)  xx[j][r]=xt[j];
                  f[r]=ft;
              }
          }
          else
          { 
			  if (ft<=f[g])
              { 
				  for (j=0; j<=n-1; j++) xx[j][r]=xt[j];
                  f[r]=ft;
              }
              else 
              { 
				  if (ft<=f[r])
                  { 
					  for (j=0; j<=n-1; j++)  xx[j][r]=xt[j];
                      f[r]=ft;
                  }
                  for (j=0; j<=n-1; j++)
                      xf[j]=v*xx[j][r]+(1.0-v)*xf[j];
                  ff=func (xf);
                  if (ff>f[r])
                      for (i=0; i<=n; i++)
					  { 
						  for (j=0; j<=n-1; j++)
						  { 
							  xx[j][i]=(xx[j][i]+xx[j][l])/2.0;
                              x[j]=xx[j][i]; xe[j]=x[j];
						  }
                          fe=func (xe);
						  f[i]=fe;
					  }
                  else
                  { 
					  for (j=0; j<=n-1; j++)  xx[j][r]=xf[j];
                      f[r]=ff;
                  }
              }
          }
          ff=0.0; ft=0.0;
          for (i=0; i<=n; i++)
          { 
			  ff=ff+f[i]/(1.0+nn);
              ft=ft+f[i]*f[i];
          }
          ft=(ft-(1.0+n)*ff*ff)/nn;
      }
      for (j=0; j<=n-1; j++)
      { 
		  x[j]=0.0;
          for (i=0; i<=n; i++) x[j]=x[j]+xx[j][i]/(1.0+nn);
          xe[j]=x[j];
      }
      fe=func (xe); x[n]=fe;
      flag = kk;
	  return;
  }

  double jsim::func (double *x)    //计算目标函数值
  {
	  double y;
	  y = x[1] - x[0]*x[0];
	  y = 100*y*y;
	  y = y + (1 - x[0])*(1 - x[0]);
	  return y;
  }

  void jsim::output ()     //输出实际迭代次数以及极小值点
		                       //与极小值到文件并显示
  {
	  int i;
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
      cout <<"flag = " <<flag <<endl;
	  fout <<flag <<endl;
	  for (i=0; i<n; i++)
	  {
		  cout <<"x(" <<i <<")=" <<x[i] <<endl;
		  fout <<x[i] <<endl;
	  }
	  cout <<"f = " <<x[n] <<endl;
	  fout <<x[n] <<endl;
	  fout.close ();
  }

  void main ()      //主函数
  {
	  jsim  s(2);           //创建对象
	  s.input ();        //由文件读入k,d,u,v,eps
	  s.solution ();     //执行单形调优法
	  s.output ();        //输出实际迭代次数以及极小值点
		                       //与极小值到文件并显示
  }

