  //9LPLQ.CPP
  //求解不等式线性规划问题
  #include  <iostream>
  #include  <iomanip>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  lplq
  {
  private: 
	       int m, n, flag;
		   double  **a, *b, *c, *x;
		   double  **p, **d;
		   int *js, *jjs, *is;
  public:
	       lplq (int mm, int nn)
		   {
			   int i;
			   m = mm; n= nn;
               b = new double[m];      //动态分配内存
	           c = new double[m+n];
	           x = new double[m+n];
	           js = new int[m];
	           jjs = new int[m];
	           is = new int[m];
	           a = new double*[m];
	           for (i=0; i<m; i++) { a[i] = new double[m+n]; }
	           p = new double*[m];
	           for (i=0; i<m; i++) { p[i] = new double[m]; }
	           d = new double*[m];
	           for (i=0; i<m; i++) { d[i] = new double[m+n]; }
		   }
	       void  input ();      //由文件读入矩阵a与向量b,c
		   void solution ();    //求解不等式线性规划问题
		   int inv ();      //矩阵求逆
           void output ();  //输出标志以及极小值点与极小值到文件并显示
		   ~lplq ()
		   {
			   int i;
			   for (i=0; i<m; i++) { delete [] a[i]; }
			   delete [] a;
			   for (i=0; i<m; i++) { delete [] p[i]; }
			   delete [] p;
			   for (i=0; i<m; i++) { delete [] d[i]; }
			   delete [] d;
			   delete [] b, c, x, js, jjs, is; 
		   }
  };

  void lplq::input ()     //由文件读入矩阵a与向量b,c
  {
	  int  i, j;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
	  for (i=0; i<m; i++)                  //读入矩阵a
	  for (j=0; j<m+n; j++)  fin >>a[i][j];
      for (i=0; i<m; i++)  fin >>b[i];     //读入向量b
      for (i=0; i<m+n; i++)  fin >>c[i];   //读入向量c
	  fin.close ();
  }

  void lplq::solution ()  //求解不等式线性规划问题
  { 
	  int i,mn,k,j;
      double s,z,dd,y;
      for (i=0; i<=m-1; i++) jjs[i]=n+i;
      mn=m+n; s=0.0;
      while (1==1)
      { 
		  for (i=0; i<=m-1; i++)
          for (j=0; j<=m-1; j++)  p[i][j]=a[i][jjs[j]];
          if (inv () == 0)
          { 
			  x[n]=s;   flag = 0;
              return;
          }
		  for (i=0; i<m; i++)
		  {
			  for (j=0; j<m+n; j++)
			  {
				  d[i][j] = 0.0;
				  for (k=0; k<m; k++)
					  d[i][j] = d[i][j] + p[i][k]*a[k][j];
			  }
		  }
          for (i=0; i<=mn-1; i++) x[i]=0.0;
          for (i=0; i<=m-1; i++)
          { 
			  s=0.0;
              for (j=0; j<=m-1; j++) s=s+p[i][j]*b[j];
              x[jjs[i]]=s;
          }
          k=-1; dd=1.0e-35;
          for (j=0; j<=mn-1; j++)
          { 
			  z=0.0;
              for (i=0; i<=m-1; i++) z=z+c[jjs[i]]*d[i][j];
              z=z-c[j];
              if (z>dd) { dd=z; k=j;}
          }
          if (k==-1)
          { 
			  s=0.0;
              for (j=0; j<=n-1; j++) s=s+c[j]*x[j];
	          x[n]=s; 
			  flag = 1;
              return;
          }
          j=-1;
          dd=1.0e+20;
          for (i=0; i<=m-1; i++)
              if (d[i][k]>=1.0e-20)
			  { 
				  y=x[jjs[i]]/d[i][k];
                  if (y<dd) { dd=y; j=i;}
			  }
          if (j==-1) 
          { 
			  x[n]=s;   flag = 0; 
			  return;
          }
          jjs[j]=k;
      }
  }

  int lplq::inv ()          //执行求逆运算
  { 
	  int i,j,k;
      double d,t;
      for (k=0; k<=m-1; k++)
      { 
		  d=0.0;
          for (i=k; i<=m-1; i++)
          for (j=k; j<=m-1; j++)
          { 
			  t=fabs(p[i][j]);
              if (t>d) { d=t; is[k]=i; js[k]=j;}
          }
          if (d+1.0==1.0)
          { 
              cout <<"\nA为奇异矩阵！没有逆矩阵. " <<endl;
			  return 0;
          }
          if (is[k]!=k)
            for (j=0; j<=m-1; j++)
            { 
                t=p[k][j]; p[k][j]=p[is[k]][j]; p[is[k]][j]=t;
            }
          if (js[k]!=k)
            for (i=0; i<=m-1; i++)
            { 
                t=p[i][k]; p[i][k]=p[i][js[k]]; p[i][js[k]]=t;
            }
          p[k][k]=1.0/p[k][k];
          for (j=0; j<=m-1; j++)
            if (j!=k)  p[k][j]=p[k][j]*p[k][k];
          for (i=0; i<=m-1; i++)
            if (i!=k)
              for (j=0; j<=m-1; j++)
                if (j!=k) p[i][j]=p[i][j]-p[i][k]*p[k][j];
          for (i=0; i<=m-1; i++)
            if (i!=k)  p[i][k]=-p[i][k]*p[k][k];
      }
      for (k=m-1; k>=0; k--)
      { 
		  if (js[k]!=k)
            for (j=0; j<=m-1; j++)
            { 
                t=p[k][j]; p[k][j]=p[js[k]][j]; p[js[k]][j]=t;
            }
          if (is[k]!=k)
            for (i=0; i<=m-1; i++)
            { 
                t=p[i][k]; p[i][k]=p[i][is[k]]; p[i][is[k]]=t;
            }
      }
      return 1;
  }



  void lplq::output () //输出标志以及极小值点与极小值到文件并显示
  {
	  int i;
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
      cout <<"flag = " <<flag <<endl;
	  fout <<flag <<endl;
	  for (i=0; i<n; i++)
	  {
		  cout <<"x(" <<i <<")=" <<x[i] <<endl;
		  fout <<x[i] <<endl;
	  }
	  cout <<"f = " <<x[n] <<endl;
	  fout <<x[n] <<endl;
	  fout.close ();
  }

  void main ()      //主函数
  {
	  lplq  s(3, 3);           //创建对象
	  s.input ();        //由文件读入矩阵a与向量b,c
	  s.solution ();     //求解不等式线性规划问题
	  s.output (); //输出标志以及极小值点与极小值到文件并显示
  }

