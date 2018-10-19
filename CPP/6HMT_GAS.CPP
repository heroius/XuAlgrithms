  //6HMT_GAS.CPP
  //Hermite-Gauss求积法
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  hmt_gas
  {
  private: 
		   double  integ;
  public:
		   void integration ();    //执行Hermite-Gauss求积法
           void output ();         //输出积分值到文件并显示
		   double func (double);   //计算被积函数值
  };

  void hmt_gas::integration ()     //执行Hermite-Gauss求积法
  { 
	  int i;
      static double t[5]={-2.02018200,-0.95857190,
                           0.0,0.95857190,2.02018200};
      static double c[5]={1.181469599,0.9865791417,
               0.9453089237,0.9865791417,1.181469599};
      integ=0.0;
      for (i=0; i<=4; i++)
      { 
		  integ = integ + c[i]*func (t[i]); 
	  }
  }

  void hmt_gas::output ()       //输出积分值到文件并显示
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

  double hmt_gas::func (double x)
  {  return x*x*exp(-x*x); }

  void main ()      //主函数
  {
	  hmt_gas  solution;     //创建对象
	  solution.integration ();    //执行Hermite-Gauss求积法
	  solution.output ();    //输出积分值到文件并显示
  }

