  //1SDET.CPP
  //求行列式值
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  sdet
  {
  private: 
           int n;
		   double  **a, det;
  public:
	       sdet (int  nn)
		   {
			   int i;
			   n = nn;
	           a = new double*[n];  //动态分配内存空间
	           for (i=0; i<n; i++) a[i] = new double[n];
		   }
	       void input ();      //从文件读入矩阵元素
		   void a_det ();      //执行求行列式值
           void output ();     //输出行列式值到文件并显示
		   ~sdet ()
		   {
			   int i;
			   for (i=0; i<n; i++) { delete [] a[i]; }
			   delete [] a;
		   }
  };

  void sdet::input ()      //从文件读入矩阵元素
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

  void sdet::a_det ()                //执行求行列式值
  { 
	  int i,j,k,is,js;
      double f,q,d;
      f=1.0; det=1.0;
      for (k=0; k<=n-2; k++)
      { 
		  q=0.0;
          for (i=k; i<=n-1; i++)
          for (j=k; j<=n-1; j++)
          { 
			  d=fabs(a[i][j]);
	          if (d>q)  { q=d; is=i; js=j;}
          }
          if (q+1.0==1.0)  det=0.0;
          else
		  {
			  if (is!=k)
			  {
				  f=-f;
                  for (j=k; j<=n-1; j++)
				  { 
                     d=a[k][j]; a[k][j]=a[is][j]; a[is][j]=d;
				  }
			  }
              if (js!=k)
			  {
				 f=-f;
                 for (i=k; i<=n-1; i++)
				 { 
                     d=a[i][js]; a[i][js]=a[i][k]; a[i][k]=d;
				 }
			  }
              det=det*a[k][k];
              for (i=k+1; i<=n-1; i++)
			  {
				 d=a[i][k]/a[k][k];
                 for (j=k+1; j<=n-1; j++)
				 {
                     a[i][j]=a[i][j]-d*a[k][j];
				 }
			  }
		  }
	  }
	  det=f*det*a[n-1][n-1];
  }

  void sdet::output ()       //输出行列式值到文件并显示
  {
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  fout <<endl <<det <<endl;
	  cout <<endl <<det <<endl;
	  fout.close ();
  }

  void main ()      //主函数
  {
	  sdet  d(4);
	  d.input ();         //从文件读入矩阵元素
	  d.a_det ();         //求行列式值
	  d.output ();        //输出行列式值到文件并显示
  }

