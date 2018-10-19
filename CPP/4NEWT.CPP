  //4NEWT.CPP
  //Newton法求方程实根
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  newt
  {
  private: 
           int Max;
		   double  x, eps;
		   double  y, dy;
  public:
	       newt(double x0, double es,int m)//提供根的初值、精度要求与最大迭代次数
		   { x = x0; eps = es; Max = m; }
		   void newt_root ();      //执行Newton法
           void output ();         //输出根到文件并显示
		   void func (double);   //计算方程左端函数f(x)的函数值与导数值
  };

  void newt::newt_root ()     //执行Newton法
  { 
	  int k;
      double d,p,x0,x1;
      k=0; x0=x;
      func (x0);
      d=eps+1.0;
      while ((d>=eps)&&(k!=Max))
	  {	
		  if (fabs(dy)+1.0==1.0)
		  { 
			  cout <<"\nf(x)的导数为0！" <<endl;
			  x = x0;
			  return;
		  }
          x1=x0-y/dy;
          func (x1);
          d=fabs(x1-x0); p=fabs(y);
          if (p>d) d=p;
          x0=x1; k = k + 1;
	  }
	  if (k == Max)
	  	  cout <<"\n迭代不收敛！" <<endl;
      x=x1;
  }

  void newt::output ()       //输出根到文件并显示
  {
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  fout <<endl <<x <<endl;
	  cout <<endl <<"x = " <<x <<endl;
	  fout.close ();
  }

  void newt::func(double x)    //计算方程左端函数f(x)的函数值与导数值
  {
	  y = x*x*(x - 1.0) - 1.0;
	  dy = 3.0*x*x - 2.0*x;
  }

  void main ()      //主函数
  {
	  newt  root(1.5,0.000001,60); 
	                //创建对象，并赋根的初值、精度要求与最大迭代次数
	  root.newt_root ();    //执行Newton法
	  root.output ();     //输出根到文件并显示
  }

//x = 1.46557
