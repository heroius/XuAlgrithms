  //6CBSV.CPP
  //Chebyshev求积法
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  cbsv
  {
  private: 
		   double  a, b, eps, integ;
  public:
	       cbsv (double aa,double bb,double es)  
			                        //顺序提供a,b,eps值的构造函数
		   { a = aa; b = bb; eps = es; }
		   void integration ();    //执行Chebyshev求积法
           void output ();         //输出积分值到文件并显示
		   double func (double);   //计算被积函数值
  };

  void cbsv::integration ()   //执行Chebyshev求积法
  { 
	  int m,i,j;
      double h,d,p,ep,g,aa,bb,s,x;
      static double t[5]={-0.8324975,-0.3745414,0.0,
                                  0.3745414,0.8324975};
      m=1;
      h=b-a; d=fabs(0.001*h);
      p=1.0e+35; ep=1.0+eps;
      while ((ep>=eps)&&(fabs(h)>d))
      { 
		  g=0.0;
          for (i=1;i<=m;i++)
          { 
			  aa=a+(i-1.0)*h; bb=a+i*h;
              s=0.0;
              for (j=0;j<=4;j++)
              { 
				  x=((bb-aa)*t[j]+(bb+aa))/2.0;
                  s=s+func (x);
              }
              g=g+s;
          }
          g=g*h/5.0;
          ep=fabs(g-p)/(1.0+fabs(g));
          p=g; m=m+1; h=(b-a)/m;
      }
      integ = g;
  }

  void cbsv::output ()       //输出积分值到文件并显示
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

  double cbsv::func (double x)
  {  return (x*x+sin(x)); }

  void main ()      //主函数
  {
	  cbsv  solution (2.5, 8.4, 0.000001); 
	                        //创建对象并顺序提供a,b,eps的值
	  solution.integration ();    //执行Chebyshev求积法
	  solution.output ();    //输出积分值到文件并显示
  }

