  //3TRDE.CPP
  //求解三对角方程组
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  trde
  {
  private: 
           int n;
		   double *b, *d;
  public:
	       trde (int nn)
		   {
			   n = nn;
	           b = new double[3*n-2];    //动态分配内存空间
	           d = new double[n];
		   }
	       void input ();  //从文件读入存放三对角矩阵的向量B以及常数向量D
		   void a_trde ();      //执行追赶法
           void output ();     //输出结果到文件并显示
		   ~trde ()
		   {
			   delete [] b;
			   delete [] d;
		   }
  };

  void trde::input ()      //从文件读入存放三对角矩阵的向量B以及常数向量D
  {
	  int  i;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
	  for (i=0; i<3*n-2; i++)  fin >>b[i];  //读入三对角矩阵A中的非0元素
	  for (i=0; i<n; i++)  fin >>d[i];       //读入常数向量d
	  fin.close ();
  }

  void trde::a_trde ()        //执行追赶法
  { 
	  int k,j,m;
      double s;
      m= 3*n-2;
      for (k=0;k<=n-2;k++)
      { 
		  j=3*k; s=b[j];
          if (fabs(s)+1.0==1.0)
          { 
			  cout <<"\n程序工作失败！无解. " <<endl;
			  return;
		  }
          b[j+1]=b[j+1]/s;
          d[k]=d[k]/s;
          b[j+3]=b[j+3]-b[j+2]*b[j+1];
          d[k+1]=d[k+1]-b[j+2]*d[k];
      }
      s=b[3*n-3];
      if (fabs(s)+1.0==1.0)
      { 
		  cout <<"\n程序工作失败！无解. " <<endl;
		  return;
	  }
      d[n-1]=d[n-1]/s;
      for (k=n-2;k>=0;k--)
          d[k]=d[k]-b[3*k+1]*d[k+1];
  }

  void trde::output ()       //输出结果到文件并显示
  {
	  int  i;
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  fout <<endl;  cout <<endl;
	  for (i=0; i<n; i++)
	  {
		  fout <<d[i] <<"   ";
		  cout <<d[i] <<"   ";
	  }
	  fout <<endl;  cout <<endl;
	  fout.close ();
  }

  void main ()      //主函数
  {
	  trde  c(5); 
	  c.input ();  //从文件读入存放三对角矩阵的向量B以及常数向量D
	  c.a_trde ();        //执行追赶法
	  c.output ();        //输出结果到文件并显示
  }

