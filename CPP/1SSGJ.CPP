  //1SSGJ.CPP
  //对称正定矩阵求逆
  #include  <iostream>
  #include  <cmath>
  #include  <fstream>
  using namespace std;
  class  ssgj
  {
  private: 
           int n;
		   double  **a;
  public:
	       ssgj (int nn)
		   {
			   int i;
			   n = nn;
	           a = new double*[n];   //动态分配内存空间
	           for (i=0; i<n; i++) a[i] = new double[n];
		   }
	       void input ();    //从文件读入矩阵A的元素
		   void inv ();      //执行求逆运算
           void output ();   //输出逆矩阵到文件
		   ~ssgj ()          //释放内存空间
		   {
			   int i;
			   for (i=0; i<n; i++) { delete [] a[i]; }
			   delete [] a;
		   }
  };

  void ssgj::input ()      //从文件读入矩阵A的元素
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

  void ssgj::inv ()         //执行求逆运算
  { 
	  int i,j,k,m;
      double w,g,*b;
      b = new double[n];
      for (k=0; k<=n-1; k++)
	  { 
		  w=a[0][0];
          if (fabs(w)+1.0==1.0)
          { 
			  delete  [] b;
			  cout <<"\n矩阵奇异！"; 
			  exit(1);
		  }
          m=n-k-1;
          for (i=1; i<=n-1; i++)
          { 
			  g=a[i][0]; b[i]=g/w;
              if (i<=m) b[i]=-b[i];
              for (j=1; j<=i; j++)
                  a[i-1][j-1]=a[i][j]+g*b[j];
          }
          a[n-1][n-1]=1.0/w;
          for (i=1; i<=n-1; i++)
              a[n-1][i-1]=b[i];
      }
      for (i=0; i<=n-2; i++)
      for (j=i+1; j<=n-1; j++)
          a[i][j]=a[j][i];
      delete [] b;
  }

  void ssgj::output ()       //输出逆矩阵到文件
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
	  fout.close ();
  }

  void main ()      //主函数
  {
	  ssgj  c(4);
	  c.input ();         //从文件读入矩阵A的元素
	  c.inv ();           //执行求逆运算
	  c.output ();        //输出逆矩阵到文件
  }
