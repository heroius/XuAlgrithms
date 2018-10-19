  //4SNSE.CPP
  //梯度法求解非线性方程组
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  snse
  {
  private: 
           int n, Max;
		   double  *x, eps;
		   double  *f, **df, F, *dF;
  public:
	       snse (int nn)
		   {
			   int i;
			   n = nn;
               x = new double[n];   //动态分配内存
	           f = new double[n];
	           dF = new double[n];
	           df = new double*[n];
	           for (i=0; i<n; i++)  df[i] = new double[n];
		   }
	       void input ();          //由文件读入精度要求eps、最大迭代次数Max
		                           //以及解的一组初值
		   void snse_root ();      //执行梯度法
           void output ();         //输出非线性方程组的一组解到文件并显示
		   void func (double *);   //计算方程组各左端函数值f以及偏导数值df
		   ~snse ()
		   {
			   int k;
			   for (k=0; k<n; k++)  { delete [] df[k]; }
			   delete [] df;
			   delete [] x, f, dF;
		   }
  };

  void snse::input ()      //由文件读入精度要求eps、最大迭代次数Max
	                       //以及解的一组初值
  {
	  int  i;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
	  fin >>eps;  fin >>Max;        //读入eps与Max
	  for (i=0; i<n; i++)  fin>>x[i];     //读入初值
	  fin.close ();
  }

  void snse::snse_root ()         //执行梯度法
  { 
	  int count,j,k;
      double d,s;
      count=Max;
	  func (x);
	  for (k=0; k<n; k++)
	  {
		  dF[k] = 0.0;
		  for (j=0; j<n; j++)
			  dF[k] = dF[k] + f[j]*df[j][k];
		  dF[k] = 2*dF[k];
	  }
      F = 0.0;
	  for (k=0; k<n; k++)  F = F+f[k]*f[k];
      while ((F>=eps)&&(count!=0))
      { 
		  count=count-1;
          d=0.0;
          for (k=0; k<=n-1; k++) d=d+dF[k]*dF[k];
          if (d+1.0==1.0) 
		  { 
			  cout <<"\n可能遇到局部极值点！改变初值再试." <<endl;
			  return;
		  }
          s=F/d;
          for (j=0; j<=n-1; j++) x[j]=x[j]-s*dF[j];
          func (x);
	      for (k=0; k<n; k++)
		  {
		      dF[k] = 0.0;
		      for (j=0; j<n; j++)
			      dF[k] = dF[k] + f[j]*df[j][k];
		      dF[k] = 2*dF[k];
		  }
          F = 0.0;
	      for (k=0; k<n; k++)  F = F+f[k]*f[k];
      }
	  cout <<"\n实际迭代次数为：" <<Max-count <<endl;
      if (count == 0)
		  cout <<"\n没有满足精度要求或可能不收敛！" <<endl;
  }

  void snse::output ()       //输出非线性方程组的一组解到文件并显示
  {
	  int k;
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  fout <<endl;  cout <<endl;
	  for (k=0; k<n; k++)
	  {
		  fout <<x[k] <<"   ";
		  cout <<x[k] <<"   ";
	  }
	  fout <<endl;  cout <<endl;
	  fout.close ();
  }

  void snse::func (double *x)//计算方程组各左端函数值f以及偏导数值df
  {
	  f[0] = x[0] - 5*x[1]*x[1] + 7*x[2]*x[2] + 12;
	  f[1] = 3*x[0]*x[1] +x[0]*x[2] - 11*x[0];
	  f[2] = 2*x[1]*x[2] + 40*x[0];
	  df[0][0] = 1;
	  df[0][1] = 10*x[1];
	  df[0][2] = 14*x[2];
	  df[1][0] = 3*x[1] + x[2] -11;
	  df[1][1] = 3*x[0];
	  df[1][2] = x[0];
	  df[2][0] = 40;
	  df[2][1] = 2*x[2];
	  df[2][2] = 2*x[1];
  }

  void main ()      //主函数
  {
	  snse  root(3); 
	  root.input ();      //由文件读入精度要求eps、最大迭代次数Max
	                       //以及解的一组初值
	  root.snse_root ();    //执行梯度法
	  root.output ();      //输出非线性方程组的一组解到文件并显示
  }

//1.00019   5.00043   -4.00039
//实际迭代281次
