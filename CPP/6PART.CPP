  //6PART.CPP
  //高振荡函数求积法
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  part 
  {
  private: 
	       int n, m;
		   double  a, b, *fa, *fb, integ_cos, integ_sin;
  public:
	       part (int nn)
		   {
			   n = nn;
               fa = new double[n+1]; //动态分配内存
	           fb = new double[n+1];
		   }
	       void input ();          //由文件读入a, b以及角频率m
		                //f(x)在积分区间左端点x=a处的函数值与各阶导数值、
		                //f(x)在积分区间右端点x=b处的函数值与各阶导数值
		   void integration ();    //执行分部积分法
           void output ();         //输出积分值到文件并显示
		   ~part ()
		   { delete  fa, fb; }
  };

  void part::input ()          //由文件读入a, b以及角频率m
		                //f(x)在积分区间左端点x=a处的函数值与各阶导数值、
		                //f(x)在积分区间右端点x=b处的函数值与各阶导数值
  {
	  int  k;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
	  fin >>a  >>b;                      //读入a,b
      fin >>m;                          //读入m
      for (k=0; k<=n; k++)  fin >>fa[k]; 
	               //读入f(x)在积分区间左端点x=a处的函数值与各阶导数值
      for (k=0; k<=n; k++)  fin >>fb[k]; 
	               //读入f(x)在积分区间右端点x=b处的函数值与各阶导数值
	  fin.close ();
  }

  void part::integration ()         //执行分部积分法
  { 
	  int mm,k,j;
      double sa[4],sb[4],ca[4],cb[4],sma,smb,cma,cmb;
      sma=sin(m*a); smb=sin(m*b);
      cma=cos(m*a); cmb=cos(m*b);
      sa[0]=sma; sa[1]=cma; sa[2]=-sma; sa[3]=-cma;
      sb[0]=smb; sb[1]=cmb; sb[2]=-smb; sb[3]=-cmb;
      ca[0]=cma; ca[1]=-sma; ca[2]=-cma; ca[3]=sma;
      cb[0]=cmb; cb[1]=-smb; cb[2]=-cmb; cb[3]=smb;
      integ_cos = 0.0;  integ_sin = 0.0;
      mm=1;
      for (k=0;k<=n-1;k++)
      { 
		  j=k;
          while (j>=4) j=j-4;
          mm=mm*m;
          integ_cos = integ_cos+(fb[k]*sb[j]-fa[k]*sa[j])/(1.0*mm);
          integ_sin = integ_sin+(fb[k]*cb[j]-fa[k]*ca[j])/(1.0*mm);
      }
      integ_sin = -integ_sin;
  }

  void part::output ()       //输出积分值到文件并显示
  {
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  fout <<integ_cos <<endl;  fout <<integ_sin <<endl;
	  cout <<"integration_cos = " <<integ_cos <<endl;
	  cout <<"integration_sin = " <<integ_sin <<endl;
	  fout.close ();
  }

  void main ()      //主函数
  {
	  part  solution(3); 
	  solution.input ();         //由文件读入a, b以及角频率m
		                //f(x)在积分区间左端点x=a处的函数值与各阶导数值、
		                //f(x)在积分区间右端点x=b处的函数值与各阶导数值
	  solution.integration ();    //执行分部求积法
	  solution.output ();    //输出积分值到文件并显示
  }

