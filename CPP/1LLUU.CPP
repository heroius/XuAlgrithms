  //1LLUU.CPP
  //矩阵三角分解
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  lluu
  {
  private: 
           int n;
		   double  **a, **l, **u;
  public:
	       lluu (int nn)
		   {
			   int i;
			   n = nn;
	           a = new double*[n];   //动态分配内存空间
	           for (i=0; i<n; i++) a[i] = new double[n];
	           l = new double*[n];
	           for (i=0; i<n; i++) l[i] = new double[n];
	           u = new double*[n];
	           for (i=0; i<n; i++) u[i] = new double[n];
		   }
	       void input ();    //由文件读入矩阵A的元素
		   void lu ();       //LU分解
           void output ();   //输出L矩阵与U矩阵到文件并显示
		   ~lluu ()
		   {
			   int i;
			   for (i=0; i<n; i++) { delete [] a[i]; }
			   delete [] a;
			   for (i=0; i<n; i++) { delete [] l[i]; }
			   delete [] l;
			   for (i=0; i<n; i++) { delete [] u[i]; }
			   delete [] u;
		   }
  };

  void lluu::input ()      //由文件读入矩阵A的元素
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

  void lluu::lu ()              //LU分解
  { 
	  int i,j,k;
      for (k=0; k<=n-2; k++)
      { 
	      if (fabs(a[k][k])+1.0==1.0)
          { 
			  cout <<"\n分解失败！" <<endl;
			  exit(1);
		  }
          for (i=k+1; i<=n-1; i++)  a[i][k]=a[i][k]/a[k][k];
          for (i=k+1; i<=n-1; i++)
          {
              for (j=k+1; j<=n-1; j++)
              {
                a[i][j]=a[i][j]-a[i][k]*a[k][j];
              }
          }
      }
      for (i=0; i<=n-1; i++)
      {
		  for (j=0; j<i; j++)
          { l[i][j]=a[i][j]; u[i][j]=0.0; }
          l[i][i]=1.0; u[i][i]=a[i][i];
          for (j=i+1; j<=n-1; j++)
          { l[i][j]=0.0; u[i][j]=a[i][j]; }
      }
  }

  void lluu::output ()       //输出L矩阵与U矩阵到文件并显示
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
			  fout <<"    " <<l[i][j];
			  cout <<"    " <<l[i][j];
		  }
		  fout <<endl;  cout <<endl;
	  }
	  fout <<endl;  cout <<endl;
	  for (i=0; i<n; i++)
	  {
		  for (j=0; j<n; j++)
		  {
			  fout <<"    " <<u[i][j];
			  cout <<"    " <<u[i][j];
		  }
		  fout <<endl;  cout <<endl;
	  }
	  fout.close ();
  }

  void main ()      //主函数
  {
	  lluu  c(4); 
	  c.input ();         //由文件读入矩阵A的元素
	  c.lu ();            //LU分解
	  c.output ();        //输出L矩阵与U矩阵到文件并显示
  }

