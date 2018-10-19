  //4DHRT.CPP
  //对分法求方程根
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  dhrt
  {
  private: 
           int m, n;
		   double  a, b, h, eps, *x;
  public:
	       dhrt(int mm, double aa, double bb, double hh, double es) //提供信息
		   { 
			   m = mm;  a = aa;  b = bb; 
			   h = hh; eps = es;
			   x = new double[m];
		   }
		   void dh_root ();        //执行对分法
           void output ();         //输出根到文件并显示
		   double func (double);   //计算方程左端函数f(x)值
		   ~dhrt () { delete [] x; }
  };

  #include "math.h"
  void dhrt::dh_root ()      //执行对分法
  { 
	  int js;
      double z,y,z1,y1,z0,y0;
      n=0; z=a; y=func (z);
      while ((z<=b+h/2.0)&&(n!=m))
      { 
		  if (fabs(y)<eps)
          { 
			  n=n+1; x[n-1]=z;
              z=z+h/2.0; y=func (z);
          }
          else
          { 
			  z1=z+h; y1=func (z1);
              if (fabs(y1)<eps)
              { 
				  n=n+1; x[n-1]=z1;
                  z=z1+h/2.0; y=func (z);
              }
              else if (y*y1>0.0)
              { y=y1; z=z1;}
              else
              { 
				  js=0;
                  while (js==0)
                  { 
					  if (fabs(z1-z)<eps)
                      {
						  n=n+1; x[n-1]=(z1+z)/2.0;
                          z=z1+h/2.0; y=func (z);
                          js=1;
                      }
                      else
                      { 
						  z0=(z1+z)/2.0; y0=func (z0);
                          if (fabs(y0)<eps)
                          { 
							  x[n]=z0; n=n+1; js=1;
                              z=z0+h/2.0; y=func (z);
                          }
                          else if ((y*y0)<0.0)
                          { z1=z0; y1=y0;}
                          else { z=z0; y=y0;}
                      }
                  }
              }
          }
      }
  }

  void dhrt::output ()       //输出根到文件并显示
  {
	  int k;
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  fout <<endl;  cout <<endl;
	  for (k=1; k<=n; k++)
	  {
		  fout <<x[k-1] <<endl;
		  cout <<"x( " <<k <<" ) = " <<x[k-1] <<endl;
	  }
	  fout.close ();
  }

  double dhrt::func(double x)    //计算方程左端函数f(x)值
  {
	  double y;
	  y = (((((x - 5.0)*x + 3.0)*x + 1.0)*x - 7.0)*x + 7.0)*x -20.0;
	  return y;
  }

  void main ()      //主函数
  {
	  dhrt  root(6,-2.0,5.0,0.2,0.000001);  //创建对象并提供信息
	  root.dh_root ();    //执行对分法
	  root.output ();     //输出根到文件并显示
  }

//x( 0 ) = -1.40246
//x( 2 ) = 4.33376
