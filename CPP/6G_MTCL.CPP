  //6G_MTCL.CPP
  //Monte-Carlo法计算一维积分
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  g_mtcl
  {
  private: 
		   double  a, b, integ, r;
  public:
	       g_mtcl (double aa,double bb)  //顺序提供a,b值的构造函数
		   { a = aa; b = bb; }
		   void integration ();    //执行Monte-Carlo法
           void output ();         //输出积分值到文件并显示
		   double func (double);   //计算被积函数值
           double rnd ();          //产生一个随机数
  };

  void g_mtcl::integration ()  //执行Monte-Carlo法
  { 
	  int m;
      double d,x,s;
      r=1.0; s=0.0; d=10000.0;
      for (m=0; m<=9999; m++)
      { 
		  x=a+(b-a)*rnd();
          s=s+func(x)/d;
      }
      s=s*(b-a);
      integ = s;
  }

  double g_mtcl::rnd ()      //产生一个随机数
  { 
	  int m;
      double s,u,v,p;
      s=65536.0; u=2053.0; v=13849.0;
      m=(int)(r/s); r=r-m*s;
      r=u*r+v; m=(int)(r/s);
      r=r-m*s; p=r/s;
      return(p);
  }

  void g_mtcl::output ()       //输出积分值到文件并显示
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

  double g_mtcl::func (double x)
  {  return (x*x+sin(x)); }

  void main ()      //主函数
  {
	  g_mtcl  solution (2.5, 8.4); 
	                        //创建对象并顺序提供a,b的值
	  solution.integration ();    //执行Monte-Carlo法
	  solution.output ();    //输出积分值到文件并显示
  }

