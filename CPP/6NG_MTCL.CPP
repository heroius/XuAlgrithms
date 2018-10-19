  //6NG_MTCL.CPP
  //Monte-Carlo法计算重积分
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  ng_mtcl
  {
  private: 
	       int  n;
		   double  *a, *b, integ, r;
		   double  *x;
  public:
	       ng_mtcl (int nn)
		   {
			   n = nn;
               a = new double[n];  //动态分配内存
	           b = new double[n];
	           x = new double[n];
		   }
	       void input ();   //由文件读入积分各层下限、上限值
		   void integration ();    //执行Monte-Carlo法
           void output ();         //输出积分值到文件并显示
		   double func ();        //计算被积函数值
           double rnd ();          //产生一个随机数
		   ~ng_mtcl ()
		   { delete [] a, b, x;  }
  };

  void ng_mtcl::input ()  //由文件读入积分各层下限、上限值
  {
	  int  k;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
      for (k=0; k<n; k++)  fin >>a[k];  //顺序读入各层下限值
      for (k=0; k<n; k++)  fin >>b[k];  //顺序读入各层上限值
	  fin.close ();
  }

  void ng_mtcl::integration ()  //执行Monte-Carlo法
  { 
	  int m,i;
      double s,d;
      r=1.0; d=10000.0; s=0.0;
      for (m=0; m<=9999; m++)
      { 
		  for (i=0; i<=n-1; i++)
              x[i]=a[i]+(b[i]-a[i])*rnd ();
          s=s+func ()/d;
      }
      for (i=0; i<=n-1; i++) s=s*(b[i]-a[i]);
      integ = s;
  }

  double ng_mtcl::rnd ()      //产生一个随机数
  { 
	  int m;
      double s,u,v,p;
      s=65536.0; u=2053.0; v=13849.0;
      m=(int)(r/s); r=r-m*s;
      r=u*r+v; m=(int)(r/s);
      r=r-m*s; p=r/s;
      return(p);
  }

  void ng_mtcl::output ()       //输出积分值到文件并显示
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

  double ng_mtcl::func ()
  {  return (x[0]*x[0] + x[1]*x[1] + x[2]*x[2]); }

  void main ()      //主函数
  {
	  ng_mtcl  solution(3); 
	  solution.input ();     //由文件读入积分各层下限、上限值
	  solution.integration ();    //执行Monte-Carlo法
	  solution.output ();    //输出积分值到文件并显示
  }

