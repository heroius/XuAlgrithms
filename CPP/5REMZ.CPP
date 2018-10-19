  //5REMZ.CPP
  //最佳一致逼近的Remez方法
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  remz
  {
  private: 
           int n;
		   double  a, b, eps, *p;
  public:
	       remz (double aa,double bb,double es,int nn)
			                     //顺序提供a,b,eps,n初值的构造函数
		   { 
			   a = aa; b = bb; eps = es; n = nn; 
			   p = new double[n+2];
		   }
		   void remez ();        //执行最佳一致逼近的Remez方法
		   double func (double);    //计算函数f(x)值
           void output ();         //输出最佳一致逼近多项式的n+1个系数
		                           //以及偏差绝对值到文件并显示
		   ~remz ()
		   {  delete [] p;  }
  };

  void remz::remez ()    //执行最佳一致逼近的Remez方法
  { 
	  int i,j,k,m,nn;
      double x[21],g[21],d,t,u,s,xx,x0,h,yy;
      if (n>19) n=19;
	  nn = n+1;
      m=nn+1; d=1.0e+35;
      for (k=0; k<=nn; k++)
      { 
		  t=cos((nn-k)*3.1415926/(1.0*nn));
          x[k]=(b+a+(b-a)*t)/2.0;
      }
      while (1==1)
      { 
		  u=1.0;
          for (i=0; i<=m-1; i++)
          { 
			  p[i]=func (x[i]);
              g[i]=-u; u=-u;
          }
          for (j=0; j<=nn-1; j++)
          { 
			  k=m; s=p[k-1]; xx=g[k-1];
              for (i=j; i<=nn-1; i++)
              { 
				  t=p[nn-i+j-1]; x0=g[nn-i+j-1];
                  p[k-1]=(s-t)/(x[k-1]-x[m-i-2]);
                  g[k-1]=(xx-x0)/(x[k-1]-x[m-i-2]);
                  k=nn-i+j; s=t; xx=x0;
              }
          }
          u=-p[m-1]/g[m-1];
          for (i=0; i<=m-1; i++)
              p[i]=p[i]+g[i]*u;
          for (j=1; j<=nn-1; j++)
          { 
			  k=nn-j; h=x[k-1]; s=p[k-1];
              for (i=m-j; i<=nn; i++)
              { 
				  t=p[i-1]; p[k-1]=s-h*t;
                  s=t; k=i;
              }
          }
          p[m-1]=fabs(u); u=p[m-1];
          if (fabs(u-d)<=eps) return;
          d=u; h=0.1*(b-a)/(1.0*nn);
          xx=a; x0=a;
          while (x0<=b)
          { 
			  s=func (x0); t=p[nn-1];
              for (i=nn-2; i>=0; i--)
                  t=t*x0+p[i];
              s=fabs(s-t);
              if (s>u) { u=s; xx=x0;}
              x0=x0+h;
          }
          s=func (xx); t=p[nn-1];
          for (i=nn-2; i>=0; i--)   t=t*xx+p[i];
          yy=s-t; i=1; j=nn+1;
          while ((j-i)!=1)
          { 
			  k=(i+j)/2;
              if (xx<x[k-1]) j=k;
              else i=k;
          }
          if (xx<x[0])
          { 
			  s=func (x[0]); t=p[nn-1];
              for (k=nn-2; k>=0; k--)  t=t*x[0]+p[k];
              s=s-t;
              if (s*yy>0.0) x[0]=xx;
              else
              { 
				  for (k=nn-1; k>=0; k--)  x[k+1]=x[k];
                  x[0]=xx;
              }
          }
          else
          { 
			  if (xx>x[nn])
              { 
				  s=func (x[nn]); t=p[nn-1];
                  for (k=nn-2; k>=0; k--)  t=t*x[nn]+p[k];
                  s=s-t;
                  if (s*yy>0.0) x[nn]=xx;
                  else
                  { 
					  for (k=0; k<=nn-1; k++)  x[k]=x[k+1];
                      x[nn]=xx;
                  }
              }
              else
              { 
				  i=i-1; j=j-1;
                  s=func (x[i]); t=p[nn-1];
                  for (k=nn-2; k>=0; k--) t=t*x[i]+p[k];
                  s=s-t;
                  if (s*yy>0.0) x[i]=xx;
                  else x[j]=xx;
              }
          }
      }
  }

  void remz::output ()       //输出最佳一致逼近多项式的n+1个系数
		                           //以及偏差绝对值到文件并显示
  {
	  int k;
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  for (k=0; k<=n; k++)
	  {
		  fout <<p[k] <<endl;
		  cout <<p[k] <<endl;
	  }
	  fout <<endl <<p[n+1] <<endl;
	  cout <<endl <<p[n+1] <<endl;
	  fout.close ();
  }

  double remz::func (double x)    //计算函数f(x)值
  {
      return exp(x);
  }

  void main ()      //主函数
  {
	  remz solution(-1.0, 1.0, 1.0e-10, 3);//创建对象
	                                   //并顺序提供a,b,eps,n的初值
	  solution.remez ();          //执行最佳一致逼近的Remez方法
	  solution.output ();          //输出最佳一致逼近多项式的n+1个系数
		                           //以及偏差绝对值到文件并显示
  }

