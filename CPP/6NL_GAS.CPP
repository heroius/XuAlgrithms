  //6NL_GAS.CPP
  //Legebdre-Gauss方法计算多重积分
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  nl_gas
  {
  private: 
	       int  n, m;
		   double  lower, upper, integ;
		   double  *x, *a, *b;
		   int  *is;
  public:
	       nl_gas (int nn, int mm)  
			   //顺序提供积分重数n与积分区间等分数m的构造函数
		   { 
			   n = nn; m = mm;
			   x = new double[n];
			   a = new double[2*n+2];
			   b = new double[n+1];
			   is = new int[2*n+2];
		   }
		   void integration ();    //执行Legebdre-Gauss方法
           void output ();         //输出积分值到文件并显示
		   double func (double *);   //计算被积函数值
		   void upper_lower(int,double *);   //计算各层上下限
		   ~nl_gas ()
		   { delete [] x, a, b, is; }
  };

  void nl_gas::integration ()   //执行Legebdre-Gauss方法
  { 
	  int mm,j,k,q,l;
      double p;
      static double t[5]={-0.9061798459,-0.5384693101,0.0,
                                 0.5384693101,0.9061798459};
      static double c[5]={0.2369268851,0.4786286705,0.5688888889,
                                         0.4786286705,0.2369268851};
      mm=1; l=1;
      a[n]=1.0; a[2*n+1]=1.0;
      while (l==1)
      { 
		  for (j=mm;j<=n;j++)
          { 
			  upper_lower (j-1, x);
              a[j-1]=0.5*(upper-lower)/m;
              b[j-1]=a[j-1]+lower;
              x[j-1]=a[j-1]*t[0]+b[j-1];
              a[n+j]=0.0;
              is[j-1]=1; is[n+j]=1;
          }
          j=n; q=1;
          while (q==1)
          { 
			  k=is[j-1];
              if (j==n) p=func(x);
              else p=1.0;
              a[n+j]=a[n+j+1]*a[j]*p*c[k-1]+a[n+j];
              is[j-1]=is[j-1]+1;
              if (is[j-1]>5)
                  if (is[n+j]>=m)
				  { 
					  j=j-1; q=1;
                      if (j==0)
					  { 
						  integ = a[n+1]*a[0];
				          return;
					  }
				  }
                  else
				  { 
					  is[n+j]=is[n+j]+1;
                      b[j-1]=b[j-1]+a[j-1]*2.0;
                      is[j-1]=1; k=is[j-1];
                      x[j-1]=a[j-1]*t[k-1]+b[j-1];
                      if (j==n) q=1;
                      else q=0;
				  }
              else
              { 
				  k=is[j-1];
                  x[j-1]=a[j-1]*t[k-1]+b[j-1];
                  if (j==n) q=1;
                  else q=0;
              }
          }
          mm=j+1;
      }
  }

  void nl_gas::output ()       //输出积分值到文件并显示
  {
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  fout <<integ <<endl;  cout <<"积分值 = " <<integ <<endl;
	  fout.close ();
  }

  double nl_gas::func (double *x)
  {  return (x[2]*x[2]); }

  void nl_gas::upper_lower (int k, double *x)
  {
	  switch (k)
	  {
	  case 0:  lower = 0.0;  upper = 1.0;  break; 
	  case 1:  lower = 0.0;
		       upper = sqrt(1.0 - x[0]*x[0]); break;
	  case 2:  lower = x[0]*x[0] + x[1]*x[1];
		       upper = sqrt(2.0 - lower);
			   lower = sqrt(lower); break;
	  default:  ;
	  }
  }

  void main ()      //主函数
  {
	  nl_gas  solution (3,4); 
	                //创建对象并顺序提供积分重数n与积分区间等分数m
	  solution.integration ();    //执行Legebdre-Gauss求积法
	  solution.output ();    //输出积分值到文件并显示
  }

