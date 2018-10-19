  //4FALS.CPP
  //试位法求方程根
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  fals
  {
  private: 
           
		   double  x_low, x_high, x, eps;
  public:
	       fals(double a, double b, double es)//提供根所在区间与精度要求
		    { x_low = a;  x_high = b; eps = es; }
		   void false_root ();        //执行试位法
           void output ();         //输出根到文件并显示
		   double func (double);   //计算方程左端函数f(x)值
  };

  void fals::false_root ()     //执行试位法
  {
	  int M;
	  double f0, f1, f;
	  M = 0;
	  f0 = func (x_low);  f1 = func (x_high);
	  if (f0*f1 > 0)
	  {
		  cout <<"\n给定区间内可能无根！" <<endl;
		  exit(1);
	  }
	  do
	  {
		  M = M + 1;
		  x = (x_low*f1 - x_high*f0)/(f1 - f0);
          f = func (x);
		  if (f*f0 < 0) { x_high = x; f1 = f; }
		  else  { x_low = x; f0 = f; }
	  } while (fabs(f) >= eps);
      cout <<"\n迭代的次数为：" <<M <<endl;
  }

  void fals::output ()       //输出根到文件并显示
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

  double fals::func(double x)    //计算方程左端函数f(x)值
  {
	  double y;
	  y = x*x*x - 2*x*x + x - 2;
	  return y;
  }

  void main ()      //主函数
  {
	  fals  root(1.0, 3.0, 0.000001); 
	                    //创建对象，并给定根所在区间与精度要求
	  root.false_root ();    //执行试位法
	  root.output ();     //输出根到文件并显示
  }

//x = 2
//迭代24次
