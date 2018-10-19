  //1MAQR.CPP
  //QR分解
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  maqr
  {
  private: 
           int m, n;
		   double  **a, **q;
  public:
	       maqr (int mm, int nn)
		   {
			   int i;
			   m = mm; n = nn;
	           a = new double*[m];   //动态分配内存空间
	           for (i=0; i<m; i++) a[i] = new double[n];
	           q = new double*[m];
	           for (i=0; i<m; i++) q[i] = new double[m];
		   }
	       void input ();    //由文件读入矩阵A的元素
		   void qr ();       //QR分解
           void output ();   //输出Q和R矩阵到文件并显示
		   ~maqr ()
		   {
			   int i;
			   for (i=0; i<m; i++) { delete [] a[i]; }
			   delete [] a;
			   for (i=0; i<m; i++) { delete [] q[i]; }
			   delete [] q;
		   }
  };

  void maqr::input ()      //由文件读入矩阵A的元素
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
	  fin.close ();
  }

  void maqr::qr ()      //QR分解
  { 
	  int i,j,k,nn,jj;
      double u,alpha,w,t;
      if (m<n)
      { 
		  cout <<"\nQR分解失败！" <<endl;
		  exit(1);
	  }
      for (i=0; i<=m-1; i++)
      for (j=0; j<=m-1; j++)
      { 
		  q[i][j]=0.0;
          if (i==j) q[i][j]=1.0;
      }
      nn=n;
      if (m==n) nn=m-1;
      for (k=0; k<=nn-1; k++)
	  {
		  
		  u=0.0;
          for (i=k; i<=m-1; i++)
          { 
			  w=fabs(a[i][k]);
              if (w>u) u=w;
          }
          alpha=0.0;
          for (i=k; i<=m-1; i++)
          { t=a[i][k]/u; alpha=alpha+t*t; }
          if (a[k][k]>0.0) u=-u;
          alpha=u*sqrt(alpha);
          if (fabs(alpha)+1.0==1.0)
          { 
			  cout <<"\nQR分解失败！" <<endl;
			  exit(1);
		  }
          u=sqrt(2.0*alpha*(alpha-a[k][k]));
          if ((u+1.0)!=1.0)
          { 
			  a[k][k]=(a[k][k]-alpha)/u;
              for (i=k+1; i<=m-1; i++) a[i][k]=a[i][k]/u;
              for (j=0; j<=m-1; j++)
              { 
				  t=0.0;
                  for (jj=k; jj<=m-1; jj++)
                      t=t+a[jj][k]*q[jj][j];
                  for (i=k; i<=m-1; i++)
                      q[i][j]=q[i][j]-2.0*t*a[i][k];
              }
              for (j=k+1; j<=n-1; j++)
              { 
				  t=0.0;
                  for (jj=k; jj<=m-1; jj++)
                      t=t+a[jj][k]*a[jj][j];
                  for (i=k; i<=m-1; i++)
                      a[i][j]=a[i][j]-2.0*t*a[i][k];
              }
              a[k][k]=alpha;
              for (i=k+1; i<=m-1; i++)  a[i][k]=0.0;
          }
      }
      for (i=0; i<=m-2; i++)
      for (j=i+1; j<=m-1;j++)
      { 
          t=q[i][j]; q[i][j]=q[j][i]; q[j][i]=t;
      }
  }

  void maqr::output ()       //输出Q和R矩阵到文件并显示
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
		  for (j=0; j<m; j++)
		  {
			  fout <<"    " <<q[i][j];
			  cout <<"    " <<q[i][j];
		  }
		  fout <<endl;  cout <<endl;
	  }
	  fout <<endl;  cout <<endl;
	  for (i=0; i<m; i++)
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
	  maqr  q(4, 3); 
	  q.input ();         //由文件读入矩阵A的元素
	  q.qr ();            //QR分解
	  q.output ();        //输出Q和R矩阵到文件并显示
  }

