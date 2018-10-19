  //1RINV.CPP
  //实矩阵求逆
  #include  <iostream>
  #include  <cmath>
  #include  <fstream>
  using namespace std;
  class  rinv
  {
  private: 
           int n;
		   double  **a;
  public:
	       rinv (int nn)
		   {
			   int i;
			   n = nn;
	           a = new double*[n];   //动态分配内存空间
	           for (i=0; i<n; i++) a[i] = new double[n];
		   }
	       void input ();    //从文件读入矩阵A的元素
		   void inv ();      //执行求逆运算
           void output ();   //输出逆矩阵到文件
		   ~rinv ()          //释放内存空间
		   {
			   int i;
			   for (i=0; i<n; i++) { delete [] a[i]; }
			   delete [] a;
		   }
  };

  void rinv::input ()      //从文件读入矩阵A的元素
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

  void rinv::inv ()          //执行求逆运算
  { 
	  int *is,*js,i,j,k;
      double d,p;
      is = new int[n];
      js = new int[n];
      for (k=0; k<=n-1; k++)
      { 
		  d=0.0;
          for (i=k; i<=n-1; i++)
          for (j=k; j<=n-1; j++)
          { 
			  p=fabs(a[i][j]);
              if (p>d) { d=p; is[k]=i; js[k]=j;}
          }
          if (d+1.0==1.0)
          { 
			  delete [] is,js;
              cout <<"\nA为奇异矩阵！没有逆矩阵. " <<endl;
			  exit(1);
          }
          if (is[k]!=k)
            for (j=0; j<=n-1; j++)
            { 
                p=a[k][j]; a[k][j]=a[is[k]][j]; a[is[k]][j]=p;
            }
          if (js[k]!=k)
            for (i=0; i<=n-1; i++)
            { 
                p=a[i][k]; a[i][k]=a[i][js[k]]; a[i][js[k]]=p;
            }
          a[k][k]=1.0/a[k][k];
          for (j=0; j<=n-1; j++)
            if (j!=k)  a[k][j]=a[k][j]*a[k][k];
          for (i=0; i<=n-1; i++)
            if (i!=k)
              for (j=0; j<=n-1; j++)
                if (j!=k) a[i][j]=a[i][j]-a[i][k]*a[k][j];
          for (i=0; i<=n-1; i++)
            if (i!=k)  a[i][k]=-a[i][k]*a[k][k];
      }
      for (k=n-1; k>=0; k--)
      { 
		  if (js[k]!=k)
            for (j=0; j<=n-1; j++)
            { 
                p=a[k][j]; a[k][j]=a[js[k]][j]; a[js[k]][j]=p;
            }
          if (is[k]!=k)
            for (i=0; i<=n-1; i++)
            { 
                p=a[i][k]; a[i][k]=a[i][is[k]]; a[i][is[k]]=p;
            }
      }
      delete [] is, js;
  }

  void rinv::output ()       //输出逆矩阵到文件
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
	  rinv  c(4);
	  c.input ();         //从文件读入矩阵A的元素
	  c.inv ();           //执行求逆运算
	  c.output ();        //输出逆矩阵到文件
  }
