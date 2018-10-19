  //2HHQR.CPP
  //用QR方法计算上H矩阵全部特征值
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  hhqr
  {
  private: 
           int n, Max;
		   double  **a, eps, *u, *v;
  public:
	       hhqr (int nn)
		   {
			   int i;
			   n = nn;
	           a = new double*[n];   //动态分配内存空间
	           for (i=0; i<n; i++) a[i] = new double[n];
	           u = new double[n];
	           v = new double[n];
		   }
	       void input ();       //从文件读入上H矩阵A以及eps,Max
		   void a_hhqr ();        //用QR方法计算全部特征值
		   void output ();     //输出全部特征值到文件并显示
		   ~hhqr ()             //释放内存空间
		   {
			   int i;
			   for (i=0; i<n; i++) { delete [] a[i]; }
			   delete [] a;
			   delete [] u, v;
		   }
  };

  void hhqr::input ()      //从文件读入上H矩阵A以及eps,Max
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

  void  hhqr::a_hhqr ()         //用QR方法计算全部特征值
  { 
	  int m,it,i,j,k,l;
      double b,c,w,g,xy,p,q,r,x,s,e,f,z,y;
      it=0; m=n;
      while (m!=0)
      { 
		  l=m-1;
          while ((l>0)&&(fabs(a[l][l-1])>eps*
	            (fabs(a[l-1][l-1])+fabs(a[l][l])))) l=l-1;
          if (l==m-1)
          { 
			  u[m-1]=a[m-1][m-1]; v[m-1]=0.0;
              m=m-1; it=0;
          }
          else if (l==m-2)
          { 
			  b=-(a[m-1][m-1]+a[m-2][m-2]);
              c=a[m-1][m-1]*a[m-2][m-2]-a[m-1][m-2]*a[m-2][m-1];
              w=b*b-4.0*c;
              y=sqrt(fabs(w));
              if (w>0.0)
              { 
				  xy=1.0;
                  if (b<0.0) xy=-1.0;
                  u[m-1]=(-b-xy*y)/2.0;
                  u[m-2]=c/u[m-1];
                  v[m-1]=0.0; v[m-2]=0.0;
              }
              else
              { 
				  u[m-1]=-b/2.0; u[m-2]=u[m-1];
                  v[m-1]=y/2.0; v[m-2]=-v[m-1];
              }
              m=m-2; it=0;
          }
          else
          { 
			  if (it>=Max)
              { 
				  cout <<"\n程序工作失败！" <<endl;
                  return;
              }
              it=it+1;
              for (j=l+2; j<=m-1; j++)  a[j][j-2]=0.0;
              for (j=l+3; j<=m-1; j++)  a[j][j-3]=0.0;
              for (k=l; k<=m-2; k++)
              { 
				  if (k!=l)
                  { 
					  p=a[k][k-1]; q=a[k+1][k-1];
                      r=0.0;
                      if (k!=m-2)  r=a[k+2][k-1];
                  }
                  else
                  { 
					  x=a[m-1][m-1]+a[m-2][m-2];
                      y=a[m-2][m-2]*a[m-1][m-1]-a[m-2][m-1]*a[m-1][m-2];
                      p=a[l][l]*(a[l][l]-x)+a[l][l+1]*a[l+1][l]+y;
                      q=a[l+1][l]*(a[l][l]+a[l+1][l+1]-x);
                      r=a[l+1][l]*a[l+2][l+1];
                  }
                  if ((fabs(p)+fabs(q)+fabs(r))!=0.0)
                  { 
					  xy=1.0;
                      if (p<0.0) xy=-1.0;
                      s=xy*sqrt(p*p+q*q+r*r);
                      if (k!=l) a[k][k-1]=-s;
                      e=-q/s; f=-r/s; x=-p/s;
                      y=-x-f*r/(p+s);
                      g=e*r/(p+s);
                      z=-x-e*q/(p+s);
                      for (j=k; j<=m-1; j++)
                      { 
                          p=x*a[k][j]+e*a[k+1][j];
                          q=e*a[k][j]+y*a[k+1][j];
                          r=f*a[k][j]+g*a[k+1][j];
                          if (k!=m-2)
                          {
                              p=p+f*a[k+2][j];
                              q=q+g*a[k+2][j];
                              r=r+z*a[k+2][j]; a[k+2][j]=r;
                          }
                          a[k+1][j]=q; a[k][j]=p;
                      }
                      j=k+3;
                      if (j>=m-1) j=m-1;
                      for (i=l; i<=j; i++)
                      { 
                          p=x*a[i][k]+e*a[i][k+1];
                          q=e*a[i][k]+y*a[i][k+1];
                          r=f*a[i][k]+g*a[i][k+1];
                          if (k!=m-2)
                          { 
                              p=p+f*a[i][k+2];
                              q=q+g*a[i][k+2];
                              r=r+z*a[i][k+2]; a[i][k+2]=r;
                          }
                          a[i][k+1]=q; a[i][k]=p;
                      }
                  }
              }
          }
      }
  }

  void hhqr::output ()  //输出全部特征值到文件并显示
  {
	  int  i;
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  for (i=0; i<n; i++)
	  {
		  fout <<u[i] <<"    " <<v[i] <<endl;
		  cout <<u[i] <<"    " <<v[i] <<endl;
	  }
	  fout <<endl;  cout <<endl;
	  fout.close ();
  }

  void main ()      //主函数
  {
	  hhqr  c(5);
	  c.input ();             //从文件读入上H矩阵A以及eps,Max
	  c.a_hhqr();               //用QR方法计算全部特征值
	  c.output ();           //输出全部特征值到文件并显示
  }

