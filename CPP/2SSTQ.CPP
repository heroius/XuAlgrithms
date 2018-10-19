  //2SSTQ.CPP
  //计算实对称三对角矩阵全部特征值和特征向量
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  sstq
  {
  private: 
           int n, Max;
		   double  eps, **q, *b, *c;
  public:
	       sstq (int nn)
		   {
			   int i;
			   n = nn;
	           q = new double*[n];  //动态分配内存空间
	           for (i=0; i<n; i++) q[i] = new double[n];
	           b = new double[n];
	           c = new double[n];
		   }
	       void input ();   //从文件读入对称三对角矩阵元素以及eps,Max
		   void a_sstq ();      //计算全部特征值和特征向量
		   void output ();     //输出全部特征值和特征向量到文件并显示
		   ~sstq ()             //释放内存空间
		   {
			   int i;
			   for (i=0; i<n; i++) { delete [] q[i]; }
			   delete [] q;
			   delete [] b, c;
		   }
  };

  void sstq::input ()  //从文件读入对称三对角矩阵元素以及eps,Max
  {
	  int  i, j;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
      for (i=0; i<n; i++)  fin >>b[i];      //读入主对角元素
	  for (i=0; i<n-1; i++)  fin >>c[i];    //读入次对角元素
	  fin >>eps;  fin >>Max;
	  for (i=0; i<n; i++)                       //置Q为单位矩阵
	  {
		  for (j=0; j<n; j++)  q[i][j] = 0.0;
		  q[i][i] = 1.0;
	  }
	  fin.close ();
  }

  void sstq::a_sstq ()       //计算全部特征值和特征向量
  { 
	  int i,j,k,m,it;
      double d,f,h,g,p,r,e,s;
      c[n-1]=0.0; d=0.0; f=0.0;
      for (j=0; j<=n-1; j++)
      { 
		  it=0;
          h=eps*(fabs(b[j])+fabs(c[j]));
          if (h>d) d=h;
          m=j;
          while ((m<=n-1)&&(fabs(c[m])>d)) m=m+1;
          if (m!=j)
          { 
			  do
              { 
				  if (it==Max)
                  { 
					  cout <<"\n程序工作失败！" <<endl;
                      return;
                  }
                  it=it+1;
                  g=b[j];
                  p=(b[j+1]-g)/(2.0*c[j]);
                  r=sqrt(p*p+1.0);
                  if (p>=0.0) b[j]=c[j]/(p+r);
                  else b[j]=c[j]/(p-r);
                  h=g-b[j];
                  for (i=j+1; i<=n-1; i++)  b[i]=b[i]-h;
                  f=f+h; p=b[m]; e=1.0; s=0.0;
                  for (i=m-1; i>=j; i--)
                  { 
					  g=e*c[i]; h=e*p;
                      if (fabs(p)>=fabs(c[i]))
                      { 
						  e=c[i]/p; r=sqrt(e*e+1.0);
                          c[i+1]=s*p*r; s=e/r; e=1.0/r;
                      }
                      else
					  { 
						  e=p/c[i]; r=sqrt(e*e+1.0);
                          c[i+1]=s*c[i]*r;
                          s=1.0/r; e=e/r;
                      }
                      p=e*b[i]-s*g;
                      b[i+1]=h+s*(e*g+s*b[i]);
                      for (k=0; k<=n-1; k++)
                      { 
                          h=q[k][i+1]; q[k][i+1]=s*q[k][i]+e*h;
                          q[k][i]=e*q[k][i]-s*h;
                      }
                  }
                  c[j]=s*p; b[j]=e*p;
              }
              while (fabs(c[j])>d);
          }
          b[j]=b[j]+f;
      }
      for (i=0; i<=n-1; i++)
      { 
		  k=i; p=b[i];
          if (i+1<=n-1)
          { 
			  j=i+1;
              while ((j<=n-1)&&(b[j]<=p))
              {  k=j; p=b[j]; j=j+1; }
          }
          if (k!=i)
          { 
			  b[k]=b[i]; b[i]=p;
              for (j=0; j<=n-1; j++)
              { 
                  p=q[j][i]; q[j][i]=q[j][k]; q[j][k]=p;
              }
          }
      }
  }

  void sstq::output ()  //输出全部特征值和特征向量到文件并显示
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
		  fout <<b[i] <<"   ";
		  cout <<b[i] <<"   ";
	  }
	  fout <<endl <<endl;  cout <<endl <<endl;
	  for (i=0; i<n; i++)
	  {
		  for (j=0; j<n; j++)
		  {
			  fout <<q[i][j] <<"   ";
			  cout <<q[i][j] <<"   ";
		  }
		  fout <<endl;  cout <<endl;
	  }
	  fout.close ();
  }

  void main ()      //主函数
  {
	  sstq  c(5);
	  c.input ();   //从文件读入对称三对角矩阵元素以及eps,Max
	  c.a_sstq();             //计算全部特征值和特征向量
	  c.output ();      //输出全部特征值和特征向量到文件并显示
  }

