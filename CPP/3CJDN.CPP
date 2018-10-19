  //3CJDN.CPP
  //全选主元Gauss-Jordan消去法求解复系数方程组
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  cjdn
  {
  private: 
           int n, m;
		   double  **ar, **ai, **br, **bi;
  public:
	       cjdn (int nn, int mm)
		   {
			   int i;
			   n = nn; m= mm;
	           ar = new double*[n];   //动态分配内存空间
	           for (i=0; i<n; i++) ar[i] = new double[n];
	           ai = new double*[n];
	           for (i=0; i<n; i++) ai[i] = new double[n];
	           br = new double*[n]; 
	           for (i=0; i<n; i++) br[i] = new double[m];
	           bi = new double*[n];
	           for (i=0; i<n; i++) bi[i] = new double[m];
		   }
	       void input ();        //从文件读入系数矩阵A以及m组常数向量B
		   void c_gauss_jordan ();      //执行Gauss-Jordan消去法
           void output ();          //输出结果到文件并显示
		   ~cjdn ()
		   {
			   int i;
			   for (i=0; i<n; i++) { delete [] ar[i]; }
			   delete [] ar;
			   for (i=0; i<n; i++) { delete [] ai[i]; }
			   delete [] ai;
			   for (i=0; i<n; i++) { delete [] br[i]; }
			   delete [] br;
			   for (i=0; i<n; i++) { delete [] bi[i]; }
			   delete [] bi;
		   }
  };

  void cjdn::input ()  //从文件读入系数矩阵A以及m组常数向量B
  {
	  int  i, j;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
	  for (i=0; i<n; i++)                       //读入矩阵A
		  for (j=0; j<n; j++)  fin >>ar[i][j];
	  for (i=0; i<n; i++)
		  for (j=0; j<n; j++)  fin >>ai[i][j];
	  for (i=0; i<m; i++)                       //读入m组常数向量B
		  for (j=0; j<n; j++)  fin >>br[j][i];
	  for (i=0; i<m; i++)
		  for (j=0; j<n; j++)  fin >>bi[j][i];
	  fin.close ();
  }

  void cjdn::c_gauss_jordan ()        //执行Gauss-Jordan消去法
  { 
	  int *js,k,i,j,is;
      double p,q,s,d;
      js = new int[n];
      for (k=0;k<=n-1;k++)
      { 
		  d=0.0;
          for (i=k;i<=n-1;i++)
          for (j=k;j<=n-1;j++)
          { 
              p=ar[i][j]*ar[i][j]+ai[i][j]*ai[i][j];
              if (p>d) {d=p;js[k]=j;is=i;}
          }
          if (d+1.0==1.0)
          { 
			  delete [] js;
			  cout <<"\n系数矩阵奇异！无解. " <<endl;
              return;
          }
          if (is!=k)
          { 
			  for (j=k;j<=n-1;j++)
              { 
                  p=ar[k][j]; ar[k][j]=ar[is][j]; ar[is][j]=p;
                  p=ai[k][j]; ai[k][j]=ai[is][j]; ai[is][j]=p;
              }
              for (j=0;j<=m-1;j++)
              { 
                  p=br[k][j]; br[k][j]=br[is][j]; br[is][j]=p;
                  p=bi[k][j]; bi[k][j]=bi[is][j]; bi[is][j]=p;
              }
          }
          if (js[k]!=k)
            for (i=0;i<=n-1;i++)
            { 
                p=ar[i][k]; ar[i][k]=ar[i][js[k]]; ar[i][js[k]]=p;
                p=ai[i][k]; ai[i][k]=ai[i][js[k]]; ai[i][js[k]]=p;
            }
          for (j=k+1;j<=n-1;j++)
          { 
              p=ar[k][j]*ar[k][k]; 
			  q=-ai[k][j]*ai[k][k];
              s=(ar[k][k]-ai[k][k])*(ar[k][j]+ai[k][j]);
              ar[k][j]=(p-q)/d; 
			  ai[k][j]=(s-p-q)/d;
          }
          for (j=0;j<=m-1;j++)
          { 
              p=br[k][j]*ar[k][k]; 
			  q=-bi[k][j]*ai[k][k];
              s=(ar[k][k]-ai[k][k])*(br[k][j]+bi[k][j]);
              br[k][j]=(p-q)/d; 
			  bi[k][j]=(s-p-q)/d;
          }
          for (i=0;i<=n-1;i++)
            if (i!=k)
            { 
                for (j=k+1;j<=n-1;j++)
                { 
                    p=ar[i][k]*ar[k][j]; 
					q=ai[i][k]*ai[k][j];
                    s=(ar[i][k]+ai[i][k])*(ar[k][j]+ai[k][j]);
                    ar[i][j]=ar[i][j]-p+q;
                    ai[i][j]=ai[i][j]-s+p+q;
                }
                for (j=0;j<=m-1;j++)
                { 
                    p=ar[i][k]*br[k][j]; 
					q=ai[i][k]*bi[k][j];
                    s=(ar[i][k]+ai[i][k])*(br[k][j]+bi[k][j]);
                    br[i][j]=br[i][j]-p+q; 
					bi[i][j]=bi[i][j]-s+p+q;
                }
            }
      }
      for (k=n-1;k>=0;k--)
        if (js[k]!=k)
          for (j=0;j<=m-1;j++)
          { 
              p=br[k][j]; br[k][j]=br[js[k]][j]; br[js[k]][j]=p;
              p=bi[k][j]; bi[k][j]=bi[js[k]][j]; bi[js[k]][j]=p;
          }
      delete [] js;
      return;
  }

  void cjdn::output ()       //输出结果到文件并显示
  {
	  int  i, j;
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  fout <<endl;  cout <<endl;
	  for (i=0; i<m; i++)
	  {
		  for (j=0; j<n; j++)
		  {
			  fout <<br[j][i] <<"   ";  cout <<br[j][i] <<"   ";
		  }
		  fout <<endl;  cout <<endl;
		  for (j=0; j<n; j++)
		  {
			  fout <<bi[j][i] <<"   ";  cout <<bi[j][i] <<"   ";
		  }
		  fout <<endl;  cout <<endl;
		  fout <<endl;  cout <<endl;
	  }
	  fout.close ();
  }

  void main ()      //主函数
  {
	  cjdn  c(4, 2);
	  c.input ();      //从文件读入系数矩阵A以及m组常数向量B
	  c.c_gauss_jordan ();           //执行Gauss-Jordan消去法
	  c.output ();            //输出结果到文件并显示
  }

