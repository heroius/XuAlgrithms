  //2STRQ.CPP
  //执行Householder变换,计算全部特征值和特征向量
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  strq
  {
  private: 
           int n, Max;
		   double  **a, eps, **q, *b, *c;
  public:
	       strq (int nn)
		   {
			   int i;
			   n = nn;
	           a = new double*[n];   //动态分配内存空间
	           for (i=0; i<n; i++) a[i] = new double[n];
	           q = new double*[n];
	           for (i=0; i<n; i++) q[i] = new double[n];
	           b = new double[n];
	           c = new double[n];
		   }
	       void input ();       //从文件读入实对称矩阵A以及eps,Max
		   void a_strq ();      //执行Householder变换
           void output1 ();     //输出对称三对角矩阵主对角
		                       //与次对角元素以及Q矩阵到文件并显示
		   void sstq ();        //计算全部特征值和特征向量
		   void output2 ();     //输出全部特征值和特征向量到文件并显示
		   ~strq ()             //释放内存空间
		   {
			   int i;
			   for (i=0; i<n; i++) { delete [] a[i]; }
			   delete [] a;
			   for (i=0; i<n; i++) { delete [] q[i]; }
			   delete [] q;
			   delete [] b, c;
		   }
  };

  void strq::input ()      //从文件读入实对称矩阵A以及eps,Max
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

  void strq::a_strq ()          //执行Householder变换
  { 
	  int i,j,k;
      double h,f,g,h2;
      for (i=0; i<=n-1; i++)
      for (j=0; j<=n-1; j++)  q[i][j]=a[i][j];
      for (i=n-1; i>=1; i--)
	  {
		  h=0.0;
          if (i>1)
          for (k=0; k<=i-1; k++)  h=h+q[i][k]*q[i][k];
          if (h+1.0==1.0)
          { 
			  c[i]=0.0;
              if (i==1)  c[i]=q[i][i-1];
              b[i]=0.0;
          }
          else
          { 
			  c[i]=sqrt(h);
              if (q[i][i-1]>0.0) c[i]=-c[i];
              h=h-q[i][i-1]*c[i];
              q[i][i-1]=q[i][i-1]-c[i];
              f=0.0;
              for (j=0; j<=i-1; j++)
              { 
				  q[j][i]=q[i][j]/h;
                  g=0.0;
                  for (k=0; k<=j; k++)  g=g+q[j][k]*q[i][k];
                  if (j+1<=i-1)
                     for (k=j+1; k<=i-1; k++)
                           g=g+q[k][j]*q[i][k];
                  c[j]=g/h;
                  f=f+g*q[j][i];
              }
              h2=f/(h+h);
              for (j=0; j<=i-1; j++)
              { 
				  f=q[i][j];
                  g=c[j]-h2*f;
                  c[j]=g;
                  for (k=0; k<=j; k++)
                  { 
                      q[j][k]=q[j][k]-f*c[k]-g*q[i][k];
                  }
              }
              b[i]=h;
          }
      }
      for (i=0; i<=n-2; i++) c[i]=c[i+1];
      c[n-1]=0.0;
      b[0]=0.0;
      for (i=0; i<=n-1; i++)
      { 
		  if ((b[i]!=0.0)&&(i-1>=0))
            for (j=0; j<=i-1; j++)
            { 
				g=0.0;
                for (k=0; k<=i-1; k++)  g=g+q[i][k]*q[k][j];
                for (k=0; k<=i-1; k++)
                { 
                    q[k][j]=q[k][j]-g*q[k][i];
                }
            }
          b[i]=q[i][i]; q[i][i]=1.0;
          if (i-1>=0)
            for (j=0; j<=i-1; j++)
            { 
				q[i][j]=0.0; q[j][i]=0.0;
			}
      }
  }

  void strq::output1 ()  //输出对称三对角矩阵主对角
	                    //与次对角元素以及Q矩阵到文件并显示
  {
	  int  i, j;
	  char str2[20];
	  cout <<"\n三对角矩阵输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  for (i=0; i<n; i++)
	  {
		  fout <<b[i] <<"   ";
		  cout <<b[i] <<"   ";
	  }
	  fout <<endl;  cout <<endl;
	  for (i=0; i<n-1; i++)
	  {
		  fout <<c[i] <<"   ";
		  cout <<c[i] <<"   ";
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

  void strq::sstq ()       //计算全部特征值和特征向量
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

  void strq::output2 ()  //输出全部特征值和特征向量到文件并显示
  {
	  int  i, j;
	  char str2[20];
	  cout <<"\n特征值和特征向量输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  for (i=0; i<n; i++)
	  {
		  fout <<"   " <<b[i];
		  cout <<"   " <<b[i];
	  }
	  fout <<endl <<endl;  cout <<endl <<endl;
	  for (i=0; i<n; i++)
	  {
		  for (j=0; j<n; j++)
		  {
			  fout <<"    " <<q[i][j];
			  cout <<"    " <<q[i][j];
		  }
		  fout <<endl;  cout <<endl;
	  }
	  fout.close ();
  }

  void main ()      //主函数
  {
	  strq  c(5); 
	  c.input ();             //从文件读入实对称矩阵A以及eps,Max
	  c.a_strq ();            //执行Householder变换
	  c.output1 ();           //输出Q矩阵以及对称三对角矩阵主对角
	                          //与次对角元素到文件并显示
	  c.sstq();             //计算全部特征值和特征向量
	  c.output2 ();           //输出全部特征值和特征向量到文件并显示
  }

