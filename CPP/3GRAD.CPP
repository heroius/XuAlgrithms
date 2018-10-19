  //3GRAD.CPP
  //共轭梯度法求解对称正定方程组
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  grad
  {
  private: 
           int n;
		   double  **a, *b, *x, eps;
  public:
	       grad (int nn)
		   {
			   int i;
			   n = nn;
	           a = new double*[n];   //动态分配内存空间
	           for (i=0; i<n; i++) a[i] = new double[n];
	           b = new double[n];
	           x = new double[n];
		   }
	       void input ();    //从文件读入系数矩阵A、常数向量B以及eps
		   void a_grad ();   //执行共轭梯度法
           void output ();   //输出结果到文件并显示
		   ~grad ()
		   {
			   int i;
			   for (i=0; i<n; i++) { delete [] a[i]; }
			   delete [] a;
			   delete [] b, x;
		   }
  };

  void grad::input ()      //从文件读入系数矩阵A、常数向量B以及eps
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
	  for (i=0; i<n; i++)  fin >>b[i];           //读入常数向量B
	  fin >>eps;
	  fin.close ();
  }

  void grad::a_grad ()     //执行共轭梯度法
  { 
	  int i, k, j;
      double *p,*r,*s,*q,alpha,beta,d,e;
      p = new double[n];
      r = new double[n];
      s = new double[n];
      q = new double[n];
      for (i=0; i<=n-1; i++)
      { x[i]=0.0; p[i]=b[i]; r[i]=b[i]; }
      i=0;
      while (i<=n-1)
      { 
	      for (k=0; k<n; k++)
		  {
		      s[k] = 0.0;
		      for (j=0; j<n; j++)
			       s[k] = s[k]+a[k][j]*p[j];
		  }
          d=0.0; e=0.0;
          for (k=0; k<=n-1; k++)
          { 
			  d=d+p[k]*b[k]; 
			  e=e+p[k]*s[k]; 
		  }
          alpha=d/e;
          for (k=0; k<=n-1; k++)
                 x[k]=x[k]+alpha*p[k];
	      for (k=0; k<n; k++)
		  {
		      q[k] = 0.0;
		      for (j=0; j<n; j++)
			       q[k] = q[k]+a[k][j]*x[j];
		  }
          d=0.0;
          for (k=0; k<=n-1; k++)
          { 
			  r[k]=b[k]-q[k]; 
			  d=d+r[k]*s[k]; 
		  }
          beta=d/e; d=0.0;
          for (k=0; k<=n-1; k++) d=d+r[k]*r[k];
          d=sqrt(d);
          if (d<eps) 
          { 
			  delete [] p, r, s, q;
			  return;
		  }
          for (k=0; k<=n-1; k++)
              p[k]=r[k]-beta*p[k];
          i=i+1;
      }
      delete [] p, r, s, q;
  }

  void grad::output ()       //输出结果到文件并显示
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
	  grad  g(4); 
	  g.input ();         //从文件读入系数矩阵A、常数向量B以及eps
	  g.a_grad ();        //执行共轭梯度法
	  g.output ();        //输出结果到文件并显示
  }

