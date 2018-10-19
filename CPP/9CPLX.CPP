  //9CPLX.CPP
  //复形调优法求约束条件下的n维极值
  #include  <iostream>
  #include  <iomanip>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  cplx
  {
  private: 
	       int m, n, k, flag;
		   double  alpha, eps, *a, *b, *x, **xx, r;
		   double  *c, *d, *w, *xt, *xf;
  public:
	       cplx (int mm, int nn)
		   {
			   int i;
			   m = mm; n = nn;
               a = new double[n];   //动态分配内存
	           b = new double[n];
	           x = new double[n+1];
	           c = new double[m];
	           d = new double[m];
	           w = new double[m];
	           xt = new double[n];
	           xf = new double[n];
	           xx = new double*[n+1];
	           for (i=0; i<n+1; i++) { xx[i] = new double[2*n]; }
		   }
	       void  input ();      //由文件读入k,alpha,eps
		                        //以及常量约束中的下限a与上限b
		                        //初始复形的第一个顶点坐标值
		   void solution ();    //执行复形调优法
		   double rnd ();      //产生随机数
           double func (double *);   //计算目标函数值 
           void lower_upper(double *);    //计算函数约束条件中的
		                                 //下限、上限与条件值
           void output ();     //输出实际迭代次数
	                   //以及极小值点与极小值到文件并显示
		   ~cplx ()
		   {
			   int i;
			   for (i=0; i<n+1; i++) { delete [] xx[i]; }
			   delete [] xx;
			   delete [] a, b, x, c, d, w, xt, xf; 
		   }
  };

  void cplx::input ()     //由文件读入k,alpha,eps
		                        //以及常量约束中的下限a与上限b
								//初始复形的第一个顶点坐标值
  {
	  int  i;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
      fin >>k >>alpha >>eps;   //读入k,alpha,eps
      for (i=0; i<n; i++)  fin >>a[i];     //读入向量a
      for (i=0; i<n; i++)  fin >>b[i];   //读入向量b
	  for (i=0; i<n; i++)  fin >>x[i];   
	                            //读入初始复形的第一个顶点坐标值
	  fin.close ();
  }

  void cplx::solution ()   //执行复形调优法
  { 
      int g,i,j,it,kt,jt,kk,rr;
      double fj,fr,fg,z;
      r=0.0;
      for (i=0; i<=n-1; i++) xx[i][0]=x[i];
      xx[n][0]=func (x);
      for (j=1; j<=2*n-1; j++)
      { 
		  for (i=0; i<=n-1; i++)
		  { 
			  xx[i][j]=a[i]+(b[i]-a[i])*rnd();
              x[i]=xx[i][j];
          }
          it=1;
          while (it==1)
          { 
			  it=0; rr=0; g=0;
              while ((rr<n)&&(g==0))
              { 
				  if ((a[rr]<=x[rr])&&(b[rr]>=x[rr])) rr=rr+1;
                  else g=1;
              }
              if (g==0)
              { 
				  lower_upper (x);
                  rr=0;
                  while ((rr<m)&&(g==0))
                  { 
					  if ((c[rr]<=w[rr])&&(d[rr]>=w[rr])) rr=rr+1;
                      else g=1;
                  }
              }
              if (g!=0)
              { 
				  for (rr=0; rr<=n-1; rr++)
                  { 
					  z=0.0;
                      for (g=0; g<=j-1; g++)  z=z+xx[rr][g]/(1.0*j);
                      xx[rr][j]=(xx[rr][j]+z)/2.0;
                      x[rr]=xx[rr][j];
                  }
                  it=1;
              }
              else xx[n][j]=func (x);
          }
      }
      kk=1; it=1;
      while (it==1)
      { 
		  it=0;
          fr=xx[n][0]; rr=0;
          for (i=1; i<=2*n-1; i++)
              if (xx[n][i]>fr)  { rr=i; fr=xx[n][i];}
          g=0; j=0; fg=xx[n][0];
          if (rr==0) { g=1; j=1; fg=xx[n][1];}
          for (i=j+1; i<=2*n-1; i++)
              if (i!=rr)
                  if (xx[n][i]>fg) { g=i; fg=xx[n][i];}
          for (i=0; i<=n-1; i++)
          { 
			  xf[i]=0.0;
              for (j=0; j<=2*n-1; j++)
                  if (j!=rr)  xf[i]=xf[i]+xx[i][j]/(2.0*n-1.0);
              xt[i]=(1.0+alpha)*xf[i]-alpha*xx[i][rr];
          }
          jt=1;
          while (jt==1)
          { 
			  jt=0;
              z=func (xt);
              while (z>fg)
              { 
				  for (i=0; i<=n-1; i++)  xt[i]=(xt[i]+xf[i])/2.0;
                  z=func (xt);
              }
              j=0;
              for (i=0; i<=n-1; i++)
              { 
				  if (a[i]>xt[i])  { xt[i]=xt[i]+0.000001; j=1;}
                  if (b[i]<xt[i])  { xt[i]=xt[i]-0.000001; j=1;}
              }
              if (j!=0) jt=1;
              else
              { 
				  lower_upper (xt);
                  j=0; kt=1;
                  while ((kt==1)&&(j<m))
                  { 
					  if ((c[j]<=w[j])&&(d[j]>=w[j])) j=j+1;
                      else kt=0;
                  }
                  if (j<m)
                  { 
					  for (i=0; i<=n-1; i++) xt[i]=(xt[i]+xf[i])/2.0;
                      jt=1;
                  }
              }
          }
          for (i=0; i<=n-1; i++) xx[i][rr]=xt[i];
          xx[n][rr]=z;
          fr=0.0; fg=0.0;
          for (j=0; j<=2*n-1; j++)
          { 
			  fj=xx[n][j];
              fr=fr+fj/(2.0*n);
              fg=fg+fj*fj;
          }
          fr=(fg-2.0*n*fr*fr)/(2.0*n-1.0);
          if (fr>=eps)
          { 
			  kk=kk+1;
              if (kk<k) it=1;
          }
      }
      for (i=0; i<=n-1; i++)
      { 
		  x[i]=0.0;
          for (j=0; j<=2*n-1; j++) x[i]=x[i]+xx[i][j]/(2.0*n);
      }
      x[n]=func (x); 
      flag = kk;
	  return;
  }

  double cplx::rnd ()      //产生一个随机数
  { 
	  int m;
      double s,u,v,p;
      s=65536.0; u=2053.0; v=13849.0;
      m=(int)(r/s); r=r-m*s;
      r=u*r+v; m=(int)(r/s);
      r=r-m*s; p=r/s;
      return(p);
  }

  void cplx::output () //输出实际迭代次数
	                   //以及极小值点与极小值到文件并显示
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

  double cplx::func (double *x)    //计算目标函数值
  {
	  double y;
	  y = -(9.0 - (x[0] - 3.0)*(x[0] - 3.0));
	  y = y*x[1]*x[1]*x[1]/(27.0*sqrt(3.0));
	  return y;
  }

  void cplx::lower_upper(double *x)     //计算函数约束条件中的
		                               //下限、上限与条件值
  {
	  c[0] = 0.0;   c[1] = 0.0;
	  d[0] = x[0]/sqrt(3.0);   d[1] = 6.0;
	  w[0] = x[1];  w[1] = x[0] + x[1]*sqrt(3.0);
  }

  void main ()      //主函数
  {
	  cplx  s(2, 2);           //创建对象
	  s.input ();        //由文件读入k,alpha,eps
		                        //以及常量约束中的下限a与上限b
	                            //初始复形的第一个顶点坐标值
	  s.solution ();     //执行复形调优法
	  s.output ();       //输出实际迭代次数
	                   //以及极小值点与极小值到文件并显示
  }

