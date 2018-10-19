  //3GSDL.CPP
  //Gauss-Seidel迭代法
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  gsdl
  {
  private: 
           int n;
		   double  **a, *b, *x, eps;
  public:
	       gsdl (int nn)
		   {
			   int i;
			   n = nn;
	           a = new double*[n];   //动态分配内存空间
	           for (i=0; i<n; i++) a[i] = new double[n];
	           b = new double[n];
	           x = new double[n];
		   }
	       void input ();    //从文件读入系数矩阵A、常数向量B以及eps
		   void a_gsdl ();    //执行Gauss-Seidel迭代法
           void output ();   //输出结果到文件并显示
		   ~gsdl ()
		   {
			   int i;
			   for (i=0; i<n; i++) { delete [] a[i]; }
			   delete [] a;
			   delete [] b, x;
		   }
  };

  void gsdl::input ()      //从文件读入系数矩阵A、常数向量B以及eps
  {
	  int  i, j;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
	  for (i=0; i<n; i++)                       //读入矩阵A
		  for (j=0; j<n; j++)  fin >>a[i][j];
	  for (i=0; i<n; i++)  fin >>b[i];          //读入常数向量B
	  fin >>eps;                              //读入eps
	  fin.close ();
  }

  void gsdl::a_gsdl ()       //执行Gauss-Seidel迭代法
  { 
	  int i,j;
      double p,t,s,q;
      for (i=0; i<=n-1; i++)
      { 
		  p=0.0; x[i]=0.0;
          for (j=0; j<=n-1; j++)
            if (i!=j)   p=p+fabs(a[i][j]);
          if (p>=fabs(a[i][i]))
          {
			  cout <<"\n程序工作失败！" <<endl;
			  return;
		  }
      }
      p=eps+1.0;
      while (p>=eps)
      { 
		  p=0.0;
          for (i=0; i<=n-1; i++)
          { 
			  t=x[i]; s=0.0;
              for (j=0; j<=n-1; j++)
              if (j!=i) s=s+a[i][j]*x[j];
              x[i]=(b[i]-s)/a[i][i];
              q=fabs(x[i]-t)/(1.0+fabs(x[i]));
              if (q>p) p=q;
          }
      }
  }

  void gsdl::output ()       //输出结果到文件并显示
  {
	  int  i;
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  fout <<endl;  cout <<endl;
	  for (i=0; i<n; i++)
	  {
		  fout <<x[i] <<"   ";
		  cout <<x[i] <<"   ";
	  }
	  fout <<endl;  cout <<endl;
	  fout.close ();
  }

  void main ()      //主函数
  {
	  gsdl  c(4);
	  c.input ();         //从文件读入系数矩阵A、常数向量B以及eps
	  c.a_gsdl ();        //执行Gauss-Seidel迭代法
	  c.output ();        //输出结果到文件并显示
  }

