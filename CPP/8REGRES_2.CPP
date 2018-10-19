  //8REGRES_2.CPP
  //多元线性回归分析程序
  #include  <iostream>
  #include  <iomanip>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  regres_2
  {
  private: 
	       int m, n, nn;
		   double  **x, *y, *a;
		   double  q, s, r, *v, u;
		   double  **b;
  public:
	       regres_2 (int mm, int nnn)
		   {
			   int i;
			   m = mm; n = nnn;
	           nn = m + 1;
               x = new double*[m];    //动态分配内存
	           for (i=0; i<m; i++)  x[i] = new double[n];
	           y = new double[n];
	           a = new double[m+1];
	           v = new double[m];
	           b = new double*[m+1];
	           for (i=0; i<=m; i++) b[i] = new double[m+1];
		   }
	       void  input (); //由文件读入n组m个自变量值X,n个随机变量观测值y
		   void analysis ();    //执行多元线性回归分析
           void chol ();       //执行Cholesky分解法
           void output ();  //输出回归系数a以及各统计量到文件并显示
		   ~regres_2 ()
		   {
			   int i;
			   for (i=0; i<m; i++) { delete [] x[i]; }
			   delete [] x;
			   for (i=0; i<=m; i++) { delete [] b[i]; }
			   delete [] b;
			   delete [] y, a, v; 
		   }
  };

  void regres_2::input ()//由文件读入n组m个自变量值X,n个随机变量观测值y
  {
	  int  i, j;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
      for (i=0; i<m; i++)  //读入n组m个自变量值X
      for (j=0; j<n; j++)  fin >>x[i][j];
      for (i=0; i<n; i++)  fin >>y[i]; //读入n个随机变量观测值
	  fin.close ();
  }

  void regres_2::analysis() //执行多元线性回归分析
  { 
	  int i,j,k;
      double e,p,yy,pp;
      b[nn-1][nn-1] = n;
      for (j=0; j<=m-1; j++)
      { 
		  p=0.0;
          for (i=0; i<=n-1; i++) p=p+x[j][i];
          b[m][j]=p;
          b[j][m]=p;
      }
      for (i=0; i<=m-1; i++)
      for (j=i; j<=m-1; j++)
      { 
		  p=0.0;
          for (k=0; k<=n-1; k++)  p=p+x[i][k]*x[j][k];
          b[j][i]=p;
          b[i][j]=p;
      }
      a[m]=0.0;
      for (i=0; i<=n-1; i++) a[m]=a[m]+y[i];
      for (i=0; i<=m-1; i++)
      { 
		  a[i]=0.0;
          for (j=0; j<=n-1; j++)  a[i]=a[i]+x[i][j]*y[j];
      }
      chol();
      yy=0.0;
      for (i=0; i<=n-1; i++)  yy=yy+y[i]/n;
      q=0.0; e=0.0; u=0.0;
      for (i=0; i<=n-1; i++)
      { 
		  p=a[m];
          for (j=0; j<=m-1; j++) p=p+a[j]*x[j][i];
          q=q+(y[i]-p)*(y[i]-p);
          e=e+(y[i]-yy)*(y[i]-yy);
          u=u+(yy-p)*(yy-p);
      }
      s=sqrt(q/n);
      r=sqrt(1.0-q/e);
      for (j=0; j<=m-1; j++)
      { 
		  p=0.0;
          for (i=0; i<=n-1; i++)
          { 
			  pp=a[m];
              for (k=0; k<=m-1; k++)
                  if (k!=j) pp=pp+a[k]*x[k][i];
              p=p+(y[i]-pp)*(y[i]-pp);
          }
          v[j]=sqrt(1.0-q/p);
      }
  }

  void regres_2::chol ()       //执行Cholesky分解法
  { 
	  int i,j,k;
      b[0][0]=sqrt(b[0][0]);
      for (j=1; j<=nn-1; j++)  b[0][j]=b[0][j]/b[0][0];
      for (i=1; i<=nn-1; i++)
      { 
          for (j=1; j<=i; j++)
          { 
              b[i][i]=b[i][i]-b[j-1][i]*b[j-1][i];
          }
          b[i][i]=sqrt(b[i][i]);
          if (i!=(nn-1))
          { 
			  for (j=i+1; j<=nn-1; j++)
              { 
                  for (k=1; k<=i; k++)
                      b[i][j]=b[i][j]-b[k-1][i]*b[k-1][j];
                  b[i][j]=b[i][j]/b[i][i];
              }
          }
      } 
      a[0]=a[0]/b[0][0];
      for (i=1; i<=nn-1; i++)
      { 
          for (k=1; k<=i; k++)
              a[i]=a[i]-b[k-1][i]*a[k-1];
          a[i]=a[i]/b[i][i];
      } 
      a[nn-1]=a[nn-1]/b[nn-1][nn-1];
      for (k=nn-1; k>=1; k--)
      { 
          for (i=k; i<=nn-1; i++)
              a[k-1]=a[k-1]-b[k-1][i]*a[i];
          a[k-1]=a[k-1]/b[k-1][k-1];
      }
  }

  void regres_2::output () //输出回归系数a以及各统计量到文件并显示
  {
	  int i;
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  for (i=0; i<=m; i++)
	  {
		  fout <<a[i] <<"    ";
		  cout <<"a(" <<i <<") = " <<a[i] <<"   ";
	  }
	  fout <<endl;  cout <<endl;
	  cout <<"偏差平方和 = " <<q <<endl;    fout <<q <<endl;
	  cout <<"平均标准偏差 = " <<s <<endl;  fout <<s <<endl;
	  cout <<"复相关系数 = " <<r <<endl;    fout <<r <<endl;
	  cout <<"回归平方和 = " <<u <<endl;    fout <<u <<endl;
	  cout <<"m个自变量的偏相关系数:" <<endl;   fout <<endl;
	  for (i=0; i<m; i++)
	  {
		  fout <<v[i] <<"    ";
		  cout <<"v(" <<i <<") = " <<v[i] <<"   ";
	  }
	  fout <<endl;  cout <<endl;
	  fout.close ();
  }

  void main ()      //主函数
  {
	  regres_2  s(3, 5);           //创建对象
	  s.input ();  //由文件读入n组m个自变量值X,n个随机变量观测值y
	  s.analysis ();     //执行多元线性回归分析
	  s.output (); //输出回归系数a以及各统计量到文件并显示
  }

