  //6LAG_GAS.CPP
  //Laguerre-Gauss求积法
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  lag_gas
  {
  private: 
		   double  integ;
  public:
		   void integration ();    //执行Laguerre-Gauss求积法
           void output ();         //输出积分值到文件并显示
		   double func (double);   //计算被积函数值
  };

  void lag_gas::integration ()   //执行Laguerre-Gauss求积法
  { 
	  int i;
      static double t[5]={0.26355990,1.41340290,
                3.59642600,7.08580990,12.64080000};
      static double c[5]={0.6790941054,1.638487956,
                 2.769426772,4.315944000,7.104896230};
      integ = 0.0;
      for (i=0; i<=4; i++)
      { 
	      integ = integ + c[i]*func (t[i]); 
	  }
  }

  void lag_gas::output ()       //输出积分值到文件并显示
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

  double lag_gas::func (double x)
  {  return x*exp(-x); }

  void main ()      //主函数
  {
	  lag_gas  solution;     //创建对象
	  solution.integration ();    //执行Laguerre-Gauss求积法
	  solution.output ();    //输出积分值到文件并显示
  }

