  //4PQRT.CPP
  //连分式法求方程根
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  pqrt
  {
  private: 
           
		   double  x, eps;
  public:
	       pqrt(double x0, double es)//提供根的初值、精度要求
		    { x = x0;  eps = es;  }
		   void pq_root ();        //执行连分式法
           void output ();         //输出根到文件并显示
		   double func (double);   //计算方程左端函数f(x)值
  };

  void pqrt::pq_root ()     //执行连分式法
  { 
	  int i,j,m,k;
      double a[10],y[10],z,h,x0,q;
      k=10; q=1.0e+35; x0=x; h=0.0;
      while (k!=0)
      { 
		  k=k-1; j=0; 
          while (j<=7)
          { 
			  if (j<=2) z=x0+0.1*j;
              else z=h;
              y[j]=func (z);
              h=z;
              if (j==0) a[0]=z;
              else
              { 
				  m=0; i=0;
                  while ((m==0)&&(i<=j-1))
                  { 
					  if (fabs(h-a[i])+1.0==1.0) m=1;
                      else h=(y[j]-y[i])/(h-a[i]);
                      i=i+1;
                  }
                  a[j]=h;
                  if (m!=0) a[j]=q;
                  h=0.0;
                  for (i=j-1; i>=0; i--)
                  { 
					  if (fabs(a[i+1]+h)+1.0==1.0) h=q;
                      else h=-y[i]/(a[i+1]+h);
                  }
                  h=h+a[0];
              }
              if (fabs(y[j])>=eps) j=j+1;
              else { j=10; k=0;}
           }
           x0=h;
       }
       x=h;
   }

  void pqrt::output ()       //输出根到文件并显示
  {
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  fout <<endl <<x <<endl;
	  cout <<endl <<"x = " <<x <<endl;
	  fout.close ();
  }

  double pqrt::func(double x)    //计算方程左端函数f(x)值
  {
	  double y;
	  y = x*x*(x - 1.0) - 1.0;
	  return y;
  }

  void main ()      //主函数
  {
	  pqrt  root(1.0, 0.000001);    //创建对象，并赋根的初值、精度要求
	  root.pq_root ();    //执行连分式法
	  root.output ();     //输出根到文件并显示
  }

//x = 1.46557
