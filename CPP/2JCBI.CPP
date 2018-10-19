  //2JCBI.CPP
  //用Jacobi方法求实对称矩阵A的特征值与特征向量
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  jcbi
  {
  private: 
           int n, Max;
		   double  **a, **v, eps;
  public:
	       jcbi (int nn)
		   {
			   int i;
			   n = nn;
	           a = new double*[n];  //动态分配内存空间
	           for (i=0; i<n; i++) a[i] = new double[n];
	           v = new double*[n];
	           for (i=0; i<n; i++) v[i] = new double[n];
		   }
	       void input ();      //从文件读入实对称矩阵A以及eps,Max
		   void jacobi ();     //执行Jacobi算法
           void output ();     //输出特征值与特征向量到文件并显示
		   ~jcbi ()
		   {
			   int i;
			   for (i=0; i<n; i++) { delete [] a[i]; }
			   delete [] a;
			   for (i=0; i<n; i++) { delete [] v[i]; }
			   delete [] v;
		   }
  };

  void jcbi::input ()   //从文件读入实对称矩阵A以及eps,Max
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
	  fin >>eps;  fin >>Max;
	  fin.close ();
  }

  void jcbi::jacobi ()      //执行Jacobi算法
  { 
	  int i,j,p,q,l;
      double fm,cn,sn,omega,x,y,d;
      l=1;
      for (i=0; i<=n-1; i++)
      { 
		  v[i][i]=1.0;
          for (j=0; j<=n-1; j++)
              if (i!=j) v[i][j]=0.0;
      }
      while (1==1)
      { 
		  fm=0.0;
          for (i=1; i<=n-1; i++)
          for (j=0; j<=i-1; j++)
          { 
			  d=fabs(a[i][j]);
              if ((i!=j)&&(d>fm))
              { 
				  fm=d; p=i; q=j;
			  }
          }
          if (fm<eps)  return;
          if (l>Max)
		  {
			  cout <<"\n程序工作失败！" <<endl;
			  return;
		  }
          l=l+1;
          x=-a[p][q]; y=(a[q][q]-a[p][p])/2.0;
          omega=x/sqrt(x*x+y*y);
          if (y<0.0) omega=-omega;
          sn=1.0+sqrt(1.0-omega*omega);
          sn=omega/sqrt(2.0*sn);
          cn=sqrt(1.0-sn*sn);
          fm=a[p][p];
          a[p][p]=fm*cn*cn+a[q][q]*sn*sn+a[p][q]*omega;
          a[q][q]=fm*sn*sn+a[q][q]*cn*cn-a[p][q]*omega;
          a[p][q]=0.0; a[q][p]=0.0;
          for (j=0; j<=n-1; j++)
          if ((j!=p)&&(j!=q))
          { 
              fm=a[p][j];
              a[p][j]=fm*cn+a[q][j]*sn;
              a[q][j]=-fm*sn+a[q][j]*cn;
          }
          for (i=0; i<=n-1; i++)
            if ((i!=p)&&(i!=q))
            { 
                fm=a[i][p];
                a[i][p]=fm*cn+a[i][q]*sn;
                a[i][q]=-fm*sn+a[i][q]*cn;
            }
          for (i=0; i<=n-1; i++)
          { 
              fm=v[i][p];
              v[i][p]=fm*cn+v[i][q]*sn;
              v[i][q]=-fm*sn+v[i][q]*cn;
          }
      }
  }

  void jcbi::output ()       //输出特征值与特征向量到文件并显示
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
		  fout <<a[i][i] <<endl;
		  cout <<a[i][i] <<endl;
	  }
	  fout <<endl;  cout <<endl;
	  for (i=0; i<n; i++)
	  {
		  for (j=0; j<n; j++)
		  {
			  fout <<"    " <<v[i][j];
			  cout <<"    " <<v[i][j];
		  }
		  fout <<endl;  cout <<endl;
	  }
	  fout <<endl;
	  fout.close ();
  }

  void main ()      //主函数
  {
	  jcbi  c(3); 
	  c.input ();         //从文件读入实对称矩阵A以及eps,Max
	  c.jacobi ();        //执行Jacobi算法
	  c.output ();        //输出特征值与特征向量到文件并显示
  }

