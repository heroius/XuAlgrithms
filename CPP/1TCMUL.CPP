  //1TCMUL.CPP
  //复矩阵相乘
  #include  <iostream>
  #include  <fstream>
  using namespace std;
  class  tcmul
  {
  private: 
           int m, n, k;
		   double  **ar, **ai, **br, **bi, **cr, **ci;
  public:
	       tcmul (int mm, int nn, int kk)
		   {
			   int i;
			   m = mm;  n = nn;  k = kk;
	           ar = new double*[m];   //动态分配内存空间
	           for (i=0; i<m; i++) ar[i] = new double[n];
	           br = new double*[n];
	           for (i=0; i<n; i++) br[i] = new double[k];
	           cr = new double*[m];
	           for (i=0; i<m; i++) cr[i] = new double[k];
	           ai = new double*[m];
	           for (i=0; i<m; i++) ai[i] = new double[n];
	           bi = new double*[n];
	           for (i=0; i<n; i++) bi[i] = new double[k];
	           ci = new double*[m];
	           for (i=0; i<m; i++) ci[i] = new double[k];
		   }
	       void input ();    //由文件读入矩阵A和B的元素
		   void mul ();      //执行C=AB
           void output ();   //矩阵C输出到文件
		   ~tcmul ()         //释放内存空间
		   {
			   int i;
			   for (i=0; i<m; i++) { delete [] ar[i]; }
			   delete [] ar;
			   for (i=0; i<n; i++) { delete [] br[i]; }
			   delete [] br;
			   for (i=0; i<m; i++) { delete [] cr[i]; }
			   delete [] cr;
			   for (i=0; i<m; i++) { delete [] ai[i]; }
			   delete [] ai;
			   for (i=0; i<n; i++) { delete [] bi[i]; }
			   delete [] bi;
			   for (i=0; i<m; i++) { delete [] ci[i]; }
			   delete [] ci;
		   }
  };

  void tcmul::input ()      //由文件读入矩阵A和B的元素
  {
	  int  i, j;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
	  for (i=0; i<m; i++)                       //读入矩阵A
		  for (j=0; j<n; j++)  fin >>ar[i][j];
	  for (i=0; i<m; i++)
		  for (j=0; j<n; j++)  fin >>ai[i][j];
	  for (i=0; i<n; i++)                       //读入矩阵B
		  for (j=0; j<k; j++)  fin >>br[i][j];
	  for (i=0; i<n; i++)
		  for (j=0; j<k; j++)  fin >>bi[i][j];
	  fin.close ();
  }

  void tcmul::mul ()            //执行C=AB
  { int i, j, l;
    double p, q, s;
    for (i=0; i<m; i++)
    for (j=0; j<k; j++)
      {
		cr[i][j] = 0.0; ci[i][j] = 0.0;
        for (l=0; l<n; l++)
          { 
            p = ar[i][l]*br[l][j];
            q = ai[i][l]*bi[l][j];
            s = (ar[i][l] + ai[i][l])*(br[l][j] + bi[l][j]);
            cr[i][j] = cr[i][j] + p - q;
            ci[i][j] = ci[i][j] + s - p - q;
          }
      }
  }

  void tcmul::output ()       //矩阵C写到文件
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
		  for (j=0; j<k; j++)
		  {
			  fout <<"    " <<cr[i][j];
			  cout <<"    " <<cr[i][j];
		  }
		  fout <<endl; cout <<endl;
	  }
	  fout <<endl; cout <<endl;
	  for (i=0; i<m; i++)
	  {
		  for (j=0; j<k; j++)
		  {
			  fout <<"    " <<ci[i][j];
			  cout <<"    " <<ci[i][j];
		  }
		  fout <<endl; cout <<endl;
	  }
	  fout.close ();
  }

  void main ()      //主函数
  {
	  tcmul  c(3, 4, 4);
	  c.input ();         //由文件读入矩阵A和B的元素
	  c.mul ();           //执行C=AB
	  c.output ();        //矩阵C写到文件
  }
