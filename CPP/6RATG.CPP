  //6RATG.CPP
  //连分式求积法
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  ratg 
  {
  private: 
		   double  a, b, eps, integ;
  public:
	       ratg (double aa, double bb, double es)//顺序提供a,b,eps值的构造函数
           { a = aa; b = bb; eps = es; }
		   void integration ();    //执行连分式求积法
           void output ();         //输出积分值到文件并显示
		   double func (double);   //计算被积函数值
  };

  void ratg::integration ()    //执行连分式求积法
  { 
	  int m,n,k,l,j;
      double h[10],bb[10],hh,t1,s1,ep,s,x,t2,g;
      m=1; n=1;
      hh=b-a; h[0]=hh;
      t1=hh*(func (a) + func (b))/2.0;
      s1=t1; bb[0]=s1; ep=1.0+eps;
      while ((ep>=eps)&&(m<=9))
      { 
		  s=0.0;
          for (k=0;k<=n-1;k++)
          { 
			  x=a+(k+0.5)*hh;
              s=s + func (x);
          }
          t2=(t1+hh*s)/2.0;
          m=m+1;
          h[m-1]=h[m-2]/2.0;
          g=t2;
          l=0; j=2;
          while ((l==0)&&(j<=m))
          { 
			  s=g-bb[j-2];
              if (fabs(s)+1.0==1.0) l=1;
              else g=(h[m-1]-h[j-2])/s;
              j=j+1;
          }
          bb[m-1]=g;
          if (l!=0) bb[m-1]=1.0e+35;
          g=bb[m-1];
          for (j=m;j>=2;j--)
             g=bb[j-2]-h[j-2]/g;
          ep=fabs(g-s1);
          s1=g; t1=t2; hh=hh/2.0; n=n+n;
      }
      integ = g;
  }

  void ratg::output ()       //输出积分值到文件并显示
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

  double ratg::func (double x)
  {  return exp(-x*x); }

  void main ()      //主函数
  {
	  ratg  solution(0.0, 4.3, 0.000001); //创建对象并顺序提供a, b, eps值
	  solution.integration ();    //执行连分式求积法
	  solution.output ();    //输出积分值到文件并显示
  }

