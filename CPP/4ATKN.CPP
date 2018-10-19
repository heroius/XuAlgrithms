  //4ATKN.CPP
  //Aitken迭代法求方程根
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  atkn
  {
  private: 
           int Max;
		   double  x, eps;
  public:
	       atkn(double x0, double es, int m)//提供根的初值、精度要求与最大迭代次数
		    { x = x0;  eps = es; Max = m; }
		   void atkn_root ();        //执行Aitken迭代法
           void output ();         //输出根到文件并显示
		   double func (double);   //计算方程左端函数φ(x)值
  };

  void atkn::atkn_root ()     //执行Aitken迭代法
  { 
	  int flag,k;
      double u,v,x0;
      k=0; x0=x; flag=0;
      while ((flag==0)&&(k!=Max))
      { 
		  k=k+1; 
          u=func (x0); v=func (u);
          if (fabs(u-v)<eps) { x0=v; flag=1; }
          else x0=v-(v-u)*(v-u)/(v-2.0*u+x0);
      }
      x=x0;
	  if (k == Max)
		  cout <<"\n迭代不收敛！" <<endl;
  }

  void atkn::output ()       //输出根到文件并显示
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

  double atkn::func(double x)    //计算方程左端函数φ(x)值
  {
	  double y;
	  y = 6.0 - x*x;
	  return y;
  }

  void main ()      //主函数
  {
	  atkn  root(0.0, 0.000001, 20); 
	                    //创建对象，并赋根的初值、精度要求与最大迭代次数
	  root.atkn_root ();    //执行Aitken迭代法
	  root.output ();     //输出根到文件并显示
  }

//x = 2
