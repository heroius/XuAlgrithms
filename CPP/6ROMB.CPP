  //6ROMB.CPP
  //Romberg求积法
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  romb 
  {
  private: 
		   double  a, b, eps, integ;
  public:
	       romb (double aa, double bb, double es)//顺序提供a,b,eps值的构造函数
           { a = aa; b = bb; eps = es; }
		   void integration ();    //执行Romberg求积法
           void output ();         //输出积分值到文件并显示
		   double func (double);   //计算被积函数值
  };

  void romb::integration ()     //执行Romberg求积法
  { 
	  int m,n,i,k;
      double y[10],h,ep,p,x,s,q;
      h=b-a;
      y[0]=h*(func (a)+func (b))/2.0;
      m=1; n=1; ep=eps+1.0;
      while ((ep>=eps)&&(m<=9))
      { 
		  p=0.0;
          for (i=0;i<=n-1;i++)
          { 
			  x=a+(i+0.5)*h;
              p=p+func (x);
          }
          p=(y[0]+h*p)/2.0;
          s=1.0;
          for (k=1;k<=m;k++)
          { 
			  s=4.0*s;
              q=(s*p-y[k-1])/(s-1.0);
              y[k-1]=p; p=q;
          }
          ep=fabs(q-y[m-1]);
          m=m+1; y[m-1]=q; n=n+n; h=h/2.0;
      }
      integ = q;
  }

  void romb::output ()       //输出积分值到文件并显示
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

  double romb::func (double x)
  {  return x/(4+x*x); }

  void main ()      //主函数
  {
	  romb  solution(0.0, 1.0, 0.000001); //创建对象并顺序提供a, b, eps值
	  solution.integration ();    //执行Romberg求积法
	  solution.output ();    //输出积分值到文件并显示
  }

