  //1CHOL.CPP
  //Cholesky分解与求行列式值
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  chol
  {
  private: 
           int n;
		   double  **a, det;
  public:
	       chol (int nn)
		   {
			   int i;
			   n = nn;
	           a = new double*[n];  //动态分配内存空间
	           for (i=0; i<n; i++) a[i] = new double[n];
		   }
	       void input ();       //由文件读入矩阵A的元素
		   void a_chol ();      //执行Cholesky分解与求行列式值
           void output ();      //输出L矩阵与行列式值到文件并显示
		   ~chol ()
		   {
			   int i;
			   for (i=0; i<n; i++) { delete [] a[i]; }
			   delete [] a;
		   }
  };

  void chol::input ()      //由文件读入矩阵A的元素
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
	  fin.close ();
  }

  void chol::a_chol ()           //执行Cholesky分解与求行列式值
  { 
	  int i,j,k;
      double d;
      if ((a[0][0]+1.0==1.0)||(a[0][0]<0.0))
      { 
		  cout <<"\n分解失败！" <<endl;
		  exit(1);
	  }
      a[0][0]=sqrt(a[0][0]);
      d=a[0][0];
      for (i=1; i<=n-1; i++)  a[i][0]=a[i][0]/a[0][0];
      for (j=1; j<=n-1; j++)
      { 
          for (k=0; k<=j-1; k++)
              a[j][j]=a[j][j]-a[j][k]*a[j][k];
          if ((a[j][j]+1.0==1.0)||(a[j][j]<0.0))
		  { 
		      cout <<"\n分解失败！" <<endl;
		      exit(1);
		  }
          a[j][j]=sqrt(a[j][j]);
          d=d*a[j][j];
          for (i=j+1; i<=n-1; i++)
          { 
              for (k=0; k<=j-1; k++)
                 a[i][j]=a[i][j]-a[i][k]*a[j][k];
              a[i][j]=a[i][j]/a[j][j];
          }
      }
      det=d*d;
      for (i=0; i<=n-2; i++)
      for (j=i+1; j<=n-1; j++)  a[i][j]=0.0;
  }

  void chol::output ()       //输出L矩阵与行列式值到文件并显示
  {
	  int  i, j;
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  for (i=0; i<n; i++)
	  {
		  for (j=0; j<n; j++)
		  {
			  fout <<"    " <<a[i][j];
			  cout <<"    " <<a[i][j];
		  }
		  fout <<endl;  cout <<endl;
	  }
	  fout <<det <<endl; cout <<"det = " <<det <<endl;
	  fout.close ();
  }

  void main ()      //主函数
  {
	  chol  c(4); 
	  c.input ();            //由文件读入矩阵A的元素
	  c.a_chol ();           //Cholesky分解与求行列式值
	  c.output ();        //输出L矩阵与行列式值到文件并显示
  }

