  //1TRMUL.CPP
  //实矩阵相乘
  #include  <iostream>
  #include  <fstream>
  using namespace std;
  class  trmul
  {
  private: 
           int m, n, k;
		   double  **a, **b, **c;
  public:
	       trmul (int mm, int nn, int kk)
		   {
			   int i;
			   m = mm; n = nn; k = kk;
	           a = new double*[m];      //动态分配内存空间
	           for (i=0; i<m; i++) a[i] = new double[n];
	           b = new double*[n];
	           for (i=0; i<n; i++) b[i] = new double[k];
	           c = new double*[m];
	           for (i=0; i<m; i++) c[i] = new double[k];
		   }
	       void input ();    //从文件读入矩阵A与B
		   void mul ();      //执行C=AB
           void output ();   //矩阵C写到文件
		   ~trmul ()
		   {
			   int i;
			   for (i=0; i<m; i++) { delete [] a[i]; }
			   delete [] a;
			   for (i=0; i<n; i++) { delete [] b[i]; }
			   delete [] b;
			   for (i=0; i<m; i++) { delete [] c[i]; }
			   delete [] c;
		   }
  };

  void trmul::input ()      //从文件读入矩阵A与B
  {
	  int  i, j;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
	  for (i=0; i<m; i++)                       //读入矩阵A
		  for (j=0; j<n; j++)  fin >>a[i][j];
	  for (i=0; i<n; i++)                       //读入矩阵B
		  for (j=0; j<k; j++)  fin >>b[i][j];
	  fin.close ();
  }

  void trmul::mul ()         //执行C=AB
  { 
	  int i, j, t;
      for (i=0; i<m; i++)
      for (j=0; j<k; j++)
      { 
		  c[i][j] = 0.0;
		  for (t=0; t<n; t++)
			  c[i][j] = c[i][j] + a[i][t]*b[t][j];
      }
  }

  void trmul::output ()       //矩阵C写到文件
  {
	  int  i, j;
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  for (i=0; i<m; i++)
	  {
		  for (j=0; j<k; j++)
		  {
			  fout <<"    " <<c[i][j];
			  cout <<"    " <<c[i][j];
		  }
		  fout <<endl; cout <<endl;
	  }
	  fout.close ();
  }

  void main ()      //主函数
  {
	  trmul  c(4, 5, 3);
	  c.input ();     //从文件读入矩阵A和B
	  c.mul ();       //执行C=AB
	  c.output ();    //输出矩阵C到文件
  }

