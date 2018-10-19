  //10SP_FUNC.CPP
  //特殊函数的计算
  #include  <iostream>
  #include  <iomanip>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  sp_func
  {
  private: 	       
  public:
	       double gamma(double);               //Gamma函数
		   double incom_gamma(double,double);  //不完全Gamma函数
		   double errf(double);                //误差函数
		   double bessel1(int,double);   //第一类整数阶Bessel函数
		   double bessel2(int,double);   //第二类整数阶Bessel函数
		   double bes1(int,double);   //变型第一类整数阶Bessel函数
		   double bes2(int,double);   //变型第二类整数阶Bessel函数
           double beta(double,double,double);  //不完全Beta函数
           double bt(double,double,double);
           double gauss(double,double,double); //正态分布函数
		   double student(double,int);         //Student分布函数
		   double chi_square(double,int);      //Chi-square分布函数
		   double F(double, int, int);         //F-分布函数
		   double sinn(double);                //正弦积分
		   double coss(double);                //余弦积分
		   double expp(double);                //指数积分
		   double elp1(double,double);         //第一类椭圆积分
		   double fk(double,double);
		   double elp2(double,double);         //第二类椭圆积分
		   double ek(double,double);
  };

  double sp_func::gamma(double x)              //Gamma函数
  {
      int i;
      double y,t,s,u;
      static double a[11]={ 0.0000677106,-0.0003442342,
               0.0015397681,-0.0024467480,0.0109736958,
               -0.0002109075,0.0742379071,0.0815782188,
                         0.4118402518,0.4227843370,1.0};
      if (x<=0.0)
      { 
		  cout <<"err**x<=0!" <<endl; 
		  return -1.0;
	  }
      y=x;
      if (y<=1.0) { t=1.0/(y*(y+1.0)); y=y+2.0; }
      else if (y<=2.0) { t=1.0/y; y=y+1.0; }
      else if (y<=3.0) t=1.0;
      else
      { 
		  t=1.0;
          while (y>3.0)  { y=y-1.0; t=t*y;}
      }
      s=a[0]; u=y-2.0;
      for (i=1; i<=10; i++)  s=s*u+a[i];
      s=s*t;
      return s;
  }

  double sp_func::incom_gamma(double a, double x) //不完全Gamma函数
  {
	  int n;
      double p,q,d,s,s1,p0,q0,p1,q1,qq;
      if ((a<=0.0)||(x<0.0))
      { 
		  if (a<=0.0) cout <<"err**a<=0!" <<endl;
          if (x<0.0)  cout <<"err**x<0!" <<endl;
          return -1.0;
      }
      if (x+1.0==1.0)  return 0.0;
      if (x>1.0e+35)  return 1.0;
      q=log(x); q=a*q; qq=exp(q);
      if (x<1.0+a)
      { 
		  p=a; d=1.0/a; s=d;
          for (n=1; n<=100; n++)
          { 
			  p=1.0+p; d=d*x/p; s=s+d;
	          if (fabs(d)<fabs(s)*1.0e-07)
              { 
				  s=s*exp(-x)*qq/gamma(a);
                  return s;
              }
          }
      }
      else
      { 
		  s=1.0/x; p0=0.0; p1=1.0; q0=1.0; q1=x;
          for (n=1; n<=100; n++)
          { 
			  p0=p1+(n-a)*p0; q0=q1+(n-a)*q0;
              p=x*p0+n*p1; q=x*q0+n*q1;
              if (fabs(q)+1.0!=1.0)
              { 
				  s1=p/q; p1=p; q1=q;
                  if (fabs((s1-s)/s1)<1.0e-07)
                  { 
					  s=s1*exp(-x)*qq/gamma(a);
                      return 1.0-s;
                  }
                  s=s1;
              }
              p1=p; q1=q;
          }
      }
      cout <<"a too large !" <<endl;
      s=1.0-s*exp(-x)*qq/gamma(a);
      return s;
  }

  double sp_func::errf(double x)               //误差函数
  {
      double y;
      if (x>=0.0) y=incom_gamma(0.5,x*x);
      else  y=-incom_gamma(0.5,x*x);
	  return y;
  }

  double sp_func::bessel1(int n, double x)  //第一类整数阶Bessel函数
  {
	  int i,m;
      double t,y,z,p,q,s,b0,b1;
      static double a[6]={ 57568490574.0,-13362590354.0,
                   651619640.7,-11214424.18,77392.33017,
                   -184.9052456};
      static double b[6]={ 57568490411.0,1029532985.0,
                9494680.718,59272.64853,267.8532712,1.0};
      static double c[6]={ 72362614232.0,-7895059235.0,
                   242396853.1,-2972611.439,15704.4826,
                  -30.16036606};
      static double d[6]={ 144725228443.0,2300535178.0,
               18583304.74,99447.43394,376.9991397,1.0};
      static double e[5]={ 1.0,-0.1098628627e-02,
               0.2734510407e-04,-0.2073370639e-05,
               0.2093887211e-06};
      static double f[5]={ -0.1562499995e-01,
             0.1430488765e-03,-0.6911147651e-05,
             0.7621095161e-06,-0.934935152e-07};
      static double g[5]={ 1.0,0.183105e-02,
              -0.3516396496e-04,0.2457520174e-05,
              -0.240337019e-06};
      static double h[5]={ 0.4687499995e-01,
              -0.2002690873e-03,0.8449199096e-05,
              -0.88228987e-06,0.105787412e-06};
      t=fabs(x);
      if (n<0) n=-n;
      if (n!=1)
      { 
		  if (t<8.0)
		  { 
			  y=t*t; p=a[5]; q=b[5];
	          for (i=4; i>=0; i--) { p=p*y+a[i]; q=q*y+b[i];}
              p=p/q;
          }
          else
          { 
			  z=8.0/t; y=z*z;  p=e[4]; q=f[4];
              for (i=3; i>=0; i--) { p=p*y+e[i]; q=q*y+f[i];}
              s=t-0.785398164;
              p=p*cos(s)-z*q*sin(s);
              p=p*sqrt(0.636619772/t);
          }
      }
      if (n==0) return p;
      b0=p;
      if (t<8.0)
      { 
		  y=t*t; p=c[5]; q=d[5];
          for (i=4; i>=0; i--) { p=p*y+c[i]; q=q*y+d[i];}
          p=x*p/q;
      }
      else
      { 
		  z=8.0/t; y=z*z; p=g[4]; q=h[4];
          for (i=3; i>=0; i--) { p=p*y+g[i]; q=q*y+h[i];}
          s=t-2.356194491;
          p=p*cos(s)-z*q*sin(s);
          p=p*x*sqrt(0.636619772/t)/t;
      }
      if (n==1) return p;
      b1=p;
      if (x==0.0) return 0.0;
      s=2.0/t;
      if (t>1.0*n)
      { 
		  if (x<0.0) b1=-b1;
          for (i=1; i<=n-1; i++)  { p=s*i*b1-b0; b0=b1; b1=p;}
      }
      else
      { 
		  m=(n+(int)sqrt(40.0*n))/2;
          m=2*m;
          p=0.0; q=0.0; b0=1.0; b1=0.0;
          for (i=m-1; i>=0; i--)
          { 
			  t=s*(i+1)*b0-b1;
              b1=b0; b0=t;
              if (fabs(b0)>1.0e+10)
              { 
				  b0=b0*1.0e-10; b1=b1*1.0e-10;
                  p=p*1.0e-10; q=q*1.0e-10;
              }
              if ((i+2)%2==0) q=q+b0;
              if ((i+1)==n) p=b1;
          }
          q=2.0*q-b0; p=p/q;
      }
      if ((x<0.0)&&(n%2==1)) p=-p;
      return p;
  }

  double sp_func::bessel2(int n, double x)  //第二类整数阶Bessel函数
  {
	  int i;
      double y,z,p,q,s,b0,b1;
      static double a[6]={ -2.957821389e+9,7.062834065e+9,
             -5.123598036e+8,1.087988129e+7,-8.632792757e+4,
             2.284622733e+2};
      static double b[6]={ 4.0076544269e+10,7.452499648e+8,
           7.189466438e+6,4.74472647e+4,2.261030244e+2,1.0};
      static double c[6]={ -4.900604943e+12,1.27527439e+12,
            -5.153438139e+10,7.349264551e+8,-4.237922726e+6,
             8.511937935e+3};
      static double d[7]={ 2.49958057e+13,4.244419664e+11,
            3.733650367e+9,2.245904002e+7,1.02042605e+5,
            3.549632885e+2,1.0};
      static double e[5]={ 1.0,-0.1098628627e-02,
             0.2734510407e-04,-0.2073370639e-05,
             0.2093887211e-06};
      static double f[5]={ -0.1562499995e-01,
             0.1430488765e-03,-0.6911147651e-05,
             0.7621095161e-06,-0.934935152e-07};
      static double g[5]={ 1.0,0.183105e-02,
             -0.3516396496e-04,0.2457520174e-05,
             -0.240337019e-06};
      static double h[5]={ 0.4687499995e-01,
             -0.2002690873e-03,0.8449199096e-05,
             -0.88228987e-06,0.105787412e-06};
      if (n<0) n=-n;
      if (x<0.0) x=-x;
      if (x==0.0) return(-1.0e+70);
      if (n!=1)
      { 
		  if (x<8.0)
          { 
			  y=x*x; p=a[5]; q=b[5];
	          for (i=4; i>=0; i--)  { p=p*y+a[i]; q=q*y+b[i];}
              p=p/q+0.636619772*bessel1(0,x)*log(x);
          }
          else
          { 
			  z=8.0/x; y=z*z;
              p=e[4]; q=f[4];
              for (i=3; i>=0; i--) { p=p*y+e[i]; q=q*y+f[i];}
              s=x-0.785398164;
              p=p*sin(s)+z*q*cos(s);
              p=p*sqrt(0.636619772/x);
          }
      }
      if (n==0) return p;
      b0=p;
      if (x<8.0)
      { 
		  y=x*x; p=c[5]; q=d[6];
          for (i=4; i>=0; i--)
          { p=p*y+c[i]; q=q*y+d[i+1];}
          q=q*y+d[0];
          p=x*p/q+0.636619772*(bessel1(1,x)*log(x)-1.0/x);;
      }
      else
      { 
		  z=8.0/x; y=z*z;
          p=g[4]; q=h[4];
          for (i=3; i>=0; i--) { p=p*y+g[i]; q=q*y+h[i];}
          s=x-2.356194491;
          p=p*sin(s)+z*q*cos(s);
          p=p*sqrt(0.636619772/x);
      }
      if (n==1) return(p);
      b1=p;
      s=2.0/x;
      for (i=1; i<=n-1; i++) { p=s*i*b1-b0; b0=b1; b1=p;}
      return p;
  }

  double sp_func::bes1(int n, double x)  //变型第一类整数阶Bessel函数
  {
	  int i,m;
      double t,y,p,b0,b1,q;
      static double a[7]={ 1.0,3.5156229,3.0899424,1.2067492,
                         0.2659732,0.0360768,0.0045813};
      static double b[7]={ 0.5,0.87890594,0.51498869,
              0.15084934,0.02658773,0.00301532,0.00032411};
      static double c[9]={ 0.39894228,0.01328592,0.00225319,
                        -0.00157565,0.00916281,-0.02057706,
                         0.02635537,-0.01647633,0.00392377};
      static double d[9]={ 0.39894228,-0.03988024,-0.00362018,
                        0.00163801,-0.01031555,0.02282967,
                        -0.02895312,0.01787654,-0.00420059};
      if (n<0) n=-n;
      t=fabs(x);
      if (n!=1)
      { 
		  if (t<3.75)
          { 
			  y=(x/3.75)*(x/3.75); p=a[6];
              for (i=5; i>=0; i--)  p=p*y+a[i];
          }
          else
          { 
			  y=3.75/t; p=c[8];
              for (i=7; i>=0; i--)  p=p*y+c[i];
              p=p*exp(t)/sqrt(t);
          }
      }
      if (n==0) return p;
      q=p;
      if (t<3.75)
      { 
		  y=(x/3.75)*(x/3.75); p=b[6];
          for (i=5; i>=0; i--) p=p*y+b[i];
          p=p*t;
      }
      else
      { 
		  y=3.75/t; p=d[8];
          for (i=7; i>=0; i--) p=p*y+d[i];
          p=p*exp(t)/sqrt(t);
      }
      if (x<0.0) p=-p;
      if (n==1) return p;
      if (x==0.0) return(0.0);
      y=2.0/t; t=0.0; b1=1.0; b0=0.0;
      m=n+(int)sqrt(40.0*n);
      m=2*m;
      for (i=m; i>0; i--)
      { 
		  p=b0+i*y*b1; b0=b1; b1=p;
          if (fabs(b1)>1.0e+10)
          { 
			  t=t*1.0e-10; b0=b0*1.0e-10;  b1=b1*1.0e-10;
          }
          if (i==n) t=b0;
      }
      p=t*q/b1;
      if ((x<0.0)&&(n%2==1)) p=-p;
      return(p);
  }

  double sp_func::bes2(int n, double x)  //变型第二类整数阶Bessel函数
  {
	  int i;
      double y,p,b0,b1;
      static double a[7]={ -0.57721566,0.4227842,0.23069756,
           0.0348859,0.00262698,0.0001075,0.0000074};
      static double b[7]={ 1.0,0.15443144,-0.67278579,
           -0.18156897,-0.01919402,-0.00110404,-0.00004686};
      static double c[7]={ 1.25331414,-0.07832358,0.02189568,
           -0.01062446,0.00587872,-0.0025154,0.00053208};
      static double d[7]={ 1.25331414,0.23498619,-0.0365562,
           0.01504268,-0.00780353,0.00325614,-0.00068245};
      if (n<0) n=-n;
      if (x<0.0) x=-x;
      if (x==0.0) return(1.0e+70);
      if (n!=1)
      { 
		  if (x<=2.0)
          { 
			  y=x*x/4.0; p=a[6];
              for (i=5; i>=0; i--) p=p*y+a[i];
              p=p-bes1(0,x)*log(x/2.0);
          }
          else
          { 
			  y=2.0/x; p=c[6];
              for (i=5; i>=0; i--) p=p*y+c[i];
              p=p*exp(-x)/sqrt(x);
          }
      }
      if (n==0) return p;
      b0=p;
      if (x<=2.0)
      { 
		  y=x*x/4.0; p=b[6];
          for (i=5; i>=0; i--) p=p*y+b[i];
          p=p/x+bes1(1,x)*log(x/2.0);
      }
      else
      {   y=2.0/x; p=d[6];
          for (i=5; i>=0; i--) p=p*y+d[i];
          p=p*exp(-x)/sqrt(x);
      }
      if (n==1) return p;
      b1=p;
      y=2.0/x;
      for (i=1; i<n; i++) { p=b0+i*y*b1; b0=b1; b1=p;}
      return p;
  }

  double sp_func::beta(double a, double b, double x) //不完全Beta函数
  {
	  double y;
      if (a<=0.0)
      { printf("err**a<=0!"); return(-1.0);}
      if (b<=0.0)
      { printf("err**b<=0!"); return(-1.0);}
      if ((x<0.0)||(x>1.0))
      { printf("err**x<0 or x>1 !");
        return(1.0e+70);
      }
      if ((x==0.0)||(x==1.0)) y=0.0;
      else
      { 
		  y=a*log(x)+b*log(1.0-x);
          y=exp(y);
          y=y*gamma(a+b)/(gamma(a)*gamma(b));
      }
      if (x<(a+1.0)/(a+b+2.0)) y=y*bt(a,b,x)/a;
      else y=1.0-y*bt(b,a,1.0-x)/b;
      return y;
  }

  double sp_func::bt(double a, double b, double x)
  { 
	  int k;
      double d,p0,q0,p1,q1,s0,s1;
      p0=0.0; q0=1.0; p1=1.0; q1=1.0;
      for (k=1; k<=100; k++)
      { 
		  d=(a+k)*(a+b+k)*x;
          d=-d/((a+k+k)*(a+k+k+1.0));
          p0=p1+d*p0; q0=q1+d*q0; s0=p0/q0;
          d=k*(b-k)*x;
          d=d/((a+k+k-1.0)*(a+k+k));
          p1=p0+d*p1; q1=q0+d*q1; s1=p1/q1;
          if (fabs(s1-s0)<fabs(s1)*1.0e-07)
          return s1 ;
      }
      cout <<"a or b too big !" <<endl;
      return s1 ;
  }

  double sp_func::gauss(double a, double d, double x)//正态分布函数
  {
	  double y;
      if (d<=0.0) d=1.0e-10;
      y=0.5+0.5*errf((x-a)/(sqrt(2.0)*d));
	  return y;
  }

  double sp_func::student(double t, int n)        //Student分布函数
  {
	  double y;
      if (t<0.0) t=-t;
      y=1.0-beta(n/2.0,0.5,n/(n+t*t));
	  return y;
  }

  double sp_func::chi_square(double x, int n)    //Chi-square分布函数
  {
	  double y;
      if (x<0.0) x=-x;
      y=incom_gamma(n/2.0,x/2.0);
	  return y;
  }

  double sp_func::F(double f, int n1, int n2)        //F-分布函数
  {
	  double y;
      if (f<0.0) f=-f;
      y=beta(n2/2.0,n1/2.0,n2/(n2+n1*f));
	  return y;
  }

  double sp_func::sinn(double x)               //正弦积分
  {
	  int m,i,j;
      double s,p,ep,h,aa,bb,w,xx,g;
      static double t[5]={-0.9061798459,-0.5384693101,0.0,
                         0.5384693101,0.9061798459};
      static double c[5]={0.2369268851,0.4786286705,0.5688888889,
                        0.4786286705,0.2369268851};
      m=1;
      if (x==0) return(0.0);
      h=fabs(x);  s=fabs(0.0001*h);
      p=1.0e+35; ep=0.000001; g=0.0;
      while ((ep>=0.0000001)&&(fabs(h)>s))
      { 
		  g=0.0;
          for (i=1;i<=m;i++)
          { 
			  aa=(i-1.0)*h; bb=i*h;
              w=0.0;
              for (j=0;j<=4;j++)
              { 
				  xx=((bb-aa)*t[j]+(bb+aa))/2.0;
                  w=w+sin(xx)/xx*c[j];
              }
              g=g+w;
          }
          g=g*h/2.0;
          ep=fabs(g-p)/(1.0+fabs(g));
          p=g; m=m+1; h=fabs(x)/m;
      }
      return g;
  }

  double sp_func::coss(double x)               //余弦积分
  {
	  int m,i,j;
      double s,p,ep,h,aa,bb,w,xx,g,r,q;
      static double t[5]={-0.9061798459,-0.5384693101,0.0,
                         0.5384693101,0.9061798459};
      static double c[5]={0.2369268851,0.4786286705,0.5688888889,
                        0.4786286705,0.2369268851};
      m=1;
      if (x==0) x=1.0e-35;
      if (x<0.0) x=-x;
      r=0.57721566490153286060651;
      q=r+log(x);
      h=x; s=fabs(0.0001*h);
      p=1.0e+35; ep=0.000001; g=0.0;
      while ((ep>=0.0000001)&&(fabs(h)>s))
      { 
		  g=0.0;
          for (i=1;i<=m;i++)
          { 
			  aa=(i-1.0)*h; bb=i*h;
              w=0.0;
              for (j=0;j<=4;j++)
              { 
				  xx=((bb-aa)*t[j]+(bb+aa))/2.0;
                  w=w+(1.0-cos(xx))/xx*c[j];
              }
              g=g+w;
          }
          g=g*h/2.0;
          ep=fabs(g-p)/(1.0+fabs(g));
          p=g; m=m+1; h=x/m;
      }
      g=q-g;
      return g;
  }

  double sp_func::expp(double x)               //指数积分
  {
	  int m,i,j;
      double s,p,ep,h,aa,bb,w,xx,g,r,q;
      static double t[5]={-0.9061798459,-0.5384693101,0.0,
                         0.5384693101,0.9061798459};
      static double c[5]={0.2369268851,0.4786286705,0.5688888889,
                        0.4786286705,0.2369268851};
      m=1;
      if (x==0) x=1.0e-10;
      if (x<0.0) x=-x;
      r=0.57721566490153286060651;
      q=r+log(x);
      h=x; s=fabs(0.0001*h);
      p=1.0e+35; ep=0.000001; g=0.0;
      while ((ep>=0.0000001)&&(fabs(h)>s))
      { 
		  g=0.0;
          for (i=1;i<=m;i++)
          { 
			  aa=(i-1.0)*h; bb=i*h;
              w=0.0;
              for (j=0;j<=4;j++)
              { 
				  xx=((bb-aa)*t[j]+(bb+aa))/2.0;
                  w=w+(exp(-xx)-1.0)/xx*c[j];
              }
              g=g+w;
          }
          g=g*h/2.0;
          ep=fabs(g-p)/(1.0+fabs(g));
          p=g; m=m+1; h=x/m;
      }
      g=q+g;
      return g;
  }

  double sp_func::elp1(double k, double f)        //第一类椭圆积分
  {
	  int n;
      double pi,y,e,ff;
      if (k<0.0) k=-k;
      if (k>1.0) k=1.0/k;
      pi=3.1415926; y=fabs(f);
      n=0;
      while (y>=pi) { n=n+1; y=y-pi;}
      e=1.0;
      if (y>=pi/2.0) { n=n+1; e=-e; y=pi-y;}
      if (n==0) ff=fk(k,y);
      else
      { 
		  ff=fk(k,pi/2.0);
          ff=2.0*n*ff+e*fk(k,y);
      }
      if (f<0.0) ff=-ff;
      return ff;
  }

  double sp_func::fk(double k, double f)
  { 
	  int m,i,j;
      double s,p,ep,h,aa,bb,w,xx,g,q;
      static double t[5]={-0.9061798459,-0.5384693101,0.0,
                         0.5384693101,0.9061798459};
      static double c[5]={0.2369268851,0.4786286705,0.5688888889,
                        0.4786286705,0.2369268851};
      m=1; g=0.0;
      h=fabs(f); s=fabs(0.0001*h);
      p=1.0e+35; ep=0.000001;
      while ((ep>=0.0000001)&&(fabs(h)>s))
      { 
		  g=0.0;
          for (i=1;i<=m;i++)
          { 
			  aa=(i-1.0)*h; bb=i*h;
              w=0.0;
              for (j=0;j<=4;j++)
              { 
				  xx=((bb-aa)*t[j]+(bb+aa))/2.0;
                  q=sqrt(1.0-k*k*sin(xx)*sin(xx));
                  w=w+c[j]/q;
              }
              g=g+w;
          }
          g=g*h/2.0;
          ep=fabs(g-p)/(1.0+fabs(g));
          p=g; m=m+m; h=0.5*h;
      }
      return g;
  }

  double sp_func::elp2(double k, double f)       //第二类椭圆积分
  {
	  int n;
      double pi,y,e,ff;
      if (k<0.0) k=-k;
      if (k>1.0) k=1.0/k;
      pi=3.1415926; y=fabs(f);
      n=0;
      while (y>=pi) { n=n+1; y=y-pi;}
      e=1.0;
      if (y>=pi/2.0) { n=n+1; e=-e; y=pi-y;}
      if (n==0) ff=ek(k,y);
      else
      { 
		  ff=ek(k,pi/2.0); ff=2.0*n*ff+e*ek(k,y);
      }
      if (f<0.0) ff=-ff;
      return ff;
  }

  double sp_func::ek(double k, double f)
  { 
	  int m,i,j;
      double s,p,ep,h,aa,bb,w,xx,g,q;
      static double t[5]={-0.9061798459,-0.5384693101,0.0,
                         0.5384693101,0.9061798459};
      static double c[5]={0.2369268851,0.4786286705,0.5688888889,
                        0.4786286705,0.2369268851};
      m=1; g=0.0;
      h=fabs(f); s=fabs(0.0001*h);
      p=1.0e+35; ep=0.000001;
      while ((ep>=0.0000001)&&(fabs(h)>s))
      { 
		  g=0.0;
          for (i=1;i<=m;i++)
          { 
			  aa=(i-1.0)*h; bb=i*h;
              w=0.0;
              for (j=0;j<=4;j++)
              { 
				  xx=((bb-aa)*t[j]+(bb+aa))/2.0;
                  q=sqrt(1.0-k*k*sin(xx)*sin(xx));
                  w=w+q*c[j];
              }
              g=g+w;
          }
          g=g*h/2.0;
          ep=fabs(g-p)/(1.0+fabs(g));
          p=g; m=m+m; h=0.5*h;
      }
      return g;
  }

  void main ()      //主函数
  {
	  sp_func  fun;           //创建对象
      int i, j;
      double x, y;
//Gamma函数例
      for (i=1; i<=10; i++)
      { 
		  x=0.5*i; y=fun.gamma(x);
		  cout <<"x = " <<setw(3) <<x <<setw(20) <<"gamma(x) = " <<y <<endl;
      }
      y=fun.gamma(1.5)*fun.gamma(2.5)/fun.gamma(4.0);
      cout <<"B(1.5,2.5) = " <<y <<endl <<endl;
//不完全Gamma函数例
      static double aa[3]={0.5,5.0,50.0};
      static double t[3]={0.1,1.0,10.0};
      for (i=0; i<=2; i++)
      for (j=0; j<=2; j++)
      { 
          y=fun.incom_gamma(aa[i], t[j]);
          cout <<"gamma(" <<aa[i] <<", " <<t[j] <<") = " <<y <<endl;
      }
      cout <<endl;
//误差函数例
      for (i=0; i<=4; i++)
      { 
		  for (j=0; j<=7; j++)
          { 
			  x=0.05*(8.0*i+j); 
              cout <<setw(10) <<fun.errf(x);
          }
          cout <<endl;
      }
      x=2.0; ;
      cout <<setw(10) <<fun.errf(x) <<endl <<endl;
//第一类整数阶Bessel函数例
      for (j=0; j<=5; j++)
      { 
		  x=0.05;
          for (i=1; i<=4; i++)
          { 
			  y=fun.bessel1(j,x);
			  cout <<"n = " <<j <<"  x = " <<setw(4) <<x 
				   <<"  J(n, x) = " <<y <<endl;
              x=x*10.0;
          }
      }
      cout <<endl;
//第二类整数阶Bessel函数例
      for (j=0; j<=5; j++)
      { 
		  x=0.05;
          for (i=1; i<=4; i++)
          { 
			  y=fun.bessel2(j,x);
			  cout <<"n = " <<j <<"  x = " <<setw(5) <<x 
				   <<"  Y(n, x) = " <<y <<endl;
              x=x*10.0;
          }
      }
      cout <<endl;
//变型第一类整数阶Bessel函数例
      for (j=0; j<=5; j++)
      { 
		  x=0.05;
          for (i=1; i<=4; i++)
          { 
			  y=fun.bes1(j,x);
			  cout <<"n = " <<j <<"  x = " <<setw(5) <<x 
				   <<"  I(n, x) = " <<y <<endl;
              x=x*10.0;
          }
      }
      cout <<endl;
//变型第二类整数阶Bessel函数例
      for (j=0; j<=5; j++)
      { 
		  x=0.05;
          for (i=1; i<=4; i++)
          { 
			  y=fun.bes2(j,x);
			  cout <<"n = " <<j <<"  x = " <<setw(5) <<x 
				   <<"  K(n, x) = " <<y <<endl;
              x=x*10.0;
          }
      }
      cout <<endl;
//不完全Beta函数例
      static double a[5]={ 0.5,0.5,1.0,5.0,8.0};
      static double b[5]={ 0.5,5.0,3.0,0.5,10.0};
      x=0.0;
      for (j=0; j<=5; j++)
      { 
		  cout <<"x = " <<x <<endl;
          for (i=0; i<=4; i++)
          { 
              y=fun.beta(a[i],b[i],x);
              cout <<"B(" <<setw(3) <<a[i] <<", " 
				   <<setw(3) <<b[i] <<") = " <<y <<endl;
          }
          x=x+0.2;
      }
      cout <<endl;
//正态分布函数例
      static double aaa[2]={ -1.0,3.0};
      static double d[2]={ 0.5,15.0};
      for (i=0; i<=1; i++)
      { 
		  x=-10.0;
          for (j=0; j<=4; j++)
          { 
			  y=fun.gauss(aaa[i],d[i],x);
              cout <<"P(" <<setw(3) <<aaa[i] <<", " <<setw(4) <<d[i] 
				   <<", " <<setw(4) <<x <<") = " <<y <<endl;
              x=x+5.0;
          }
      }
      cout <<endl;
//Student分布函数例
      for (j=1; j<=5; j++)
      { 
		  x=0.5; y=fun.student(x,j);
		  cout <<"P(" <<setw(3) <<x <<", " <<setw(3) 
			   <<j <<") = " <<y <<endl;
          x=5.0; y=fun.student(x,j);
		  cout <<"P(" <<setw(3) <<x <<", " <<setw(3) 
			   <<j <<") = " <<y <<endl;
      }
      cout <<endl;
//Chi-square分布函数例
      for (j=1; j<=5; j++)
      { 
		  x=0.5; y=fun.chi_square(x,j);
		  cout <<"P(" <<setw(3) <<x <<", " <<setw(3) 
			   <<j <<") = " <<y <<endl;
          x=5.0; y=fun.chi_square(x,j);
		  cout <<"P(" <<setw(3) <<x <<", " <<setw(3) 
			   <<j <<") = " <<y <<endl;
      }
      cout <<endl;
//F-分布函数例
      static int n[2]={ 2,5};
      static int m[2]={ 3,10};
      for (i=0; i<=1; i++)
      { 
		  x=3.5;  y=fun.F(x,n[i],m[i]);
		  cout <<"P(" <<setw(4) <<x <<", " <<setw(4) <<n[i] 
				   <<", " <<setw(4) <<m[i] <<") = " <<y <<endl;
          x=9.0;  y=fun.F(x,n[i],m[i]);
		  cout <<"P(" <<setw(4) <<x <<", " <<setw(4) <<n[i] 
				   <<", " <<setw(4) <<m[i] <<") = " <<y <<endl;
      }
      cout <<endl;
//正弦积分例
      for (i=0; i<=9; i++)
      { 
		  x=0.5+i+i; y=fun.sinn(x);
		  cout <<"x = " <<setw(6) <<x <<"    Si(x) = " <<y <<endl;
      }
      cout <<endl;
//余弦积分例
      for (i=0; i<=9; i++)
      { 
		  x=0.5+i+i; y=fun.coss(x);
		  cout <<"x = " <<setw(6) <<x <<"    Ci(x) = " <<y <<endl;
      }
      cout <<endl;
//指数积分例
      for (i=0; i<=9; i++)
      { 
		  x=0.05+0.2*i; y=fun.expp(x);
		  cout <<"x = " <<setw(6) <<x <<"    Ei(x) = " <<y <<endl;
      }
      cout <<endl;
//第一类椭圆积分例
      for (i=0; i<=10; i++)
      { 
		  x=i*3.1415926/18.0;
          y=fun.elp1(0.5,x);
		  cout <<"F(" <<setw(4) <<0.5 <<", " <<setw(10) <<x
			   <<") = " <<setw(10) <<y <<"        ";
          y=fun.elp1(1.0,x);
		  cout <<"F(" <<setw(4) <<1.0 <<", " <<setw(10) <<x
			   <<") = " <<setw(10) <<y <<endl;
      }
      cout <<endl;
//第二类椭圆积分例
      for (i=0; i<=10; i++)
      { 
		  x=i*3.1415926/18.0;
          y=fun.elp2(0.5,x);
		  cout <<"E(" <<setw(4) <<0.5 <<", " <<setw(10) <<x
			   <<") = " <<setw(10) <<y <<"        ";
          y=fun.elp2(1.0,x);
		  cout <<"E(" <<setw(4) <<1.0 <<", " <<setw(10) <<x
			   <<") = " <<setw(10) <<y <<endl;
      }
      cout <<endl;
  }

