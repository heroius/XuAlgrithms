  //9RAT_MAX1.CPP
  //一维极值连分式法
  #include  <iostream>
  #include  <iomanip>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  rat_max1
  {
  private: 
	       int flag;
		   double  x, y, dy, eps;
  public:
	       rat_max1 (double x0, double es) 
			                  //顺序提供x0,eps值的构造函数
		   { x = x0;  eps = es;  }
		   void solution ();    //执行一维极值连分式法
		   void func (double);      //计算目标函数值与导数值
           void output ();  //输出标志与极值点以及函数值到文件并显示
  };

  void rat_max1::solution ()   //执行一维极值连分式法
  { 
	  int i,j,m,jt,k;
      double xx,h1,h2,dx,yy[10],b[10];
      k=0; jt=1; h2=0.0;
      while (jt==1)
      { 
		  j=0;
          while (j<=7)
          { 
			  if (j<=2) xx=x+0.01*j;
              else xx=h2;
              func (xx);
              if (fabs(dy)<eps) { jt=0; j=10; }
              else
              { 
				  h1=dy; h2=xx;
                  if (j==0) { yy[0]=h1; b[0]=h2;}
                  else
                  { 
					  yy[j]=h1; m=0; i=0;
                      while ((m==0)&&(i<=j-1))
					  { 
						  if (fabs(h2-b[i])+1.0==1.0) m=1;
			              else h2=(h1-yy[i])/(h2-b[i]);
                          i=i+1;
					  }
                      b[j]=h2;
                      if (m!=0) b[j]=1.0e+35;
                      h2=0.0;
                      for (i=j-1; i>=0; i--) h2=-yy[i]/(b[i+1]+h2);
                      h2=h2+b[0];
                  }
                  j=j+1;
              }
          }
          x=h2;
          k = k + 1;
          if (k == 10)  jt = 0;
      }
      xx=x;  yy[0] = y;
      func (xx); 
      if (fabs(x)<=1.0) dx=1.0e-05;
      else dx=fabs(x*1.0e-05);
      xx=x-dx;
      func (xx); h1=y;
      xx=x+dx;
      func (xx); h2=y;
      flag = -1;
      if ((h1+h2-2.0*yy[0])<=0.0) flag = 1;
	  y = yy[0];
  }

  void rat_max1::output () //输出标志与极值点以及函数值到文件并显示
  {
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
      cout <<"flag = " <<flag <<endl;
	  fout <<flag <<endl;
	  if (flag == 1)  cout <<"MAX:  ";
	  else  cout <<"MIN:  ";
	  cout <<"x = " <<x <<"     ";
	  fout <<x <<"    ";
	  cout <<"f(x) = "  <<y <<endl;
	  fout <<y <<endl;
	  fout.close ();
  }

  void rat_max1::func (double x)  //计算目标函数值与导数值
  {
	  y = (x - 1.0)*(10.0 - x);
	  dy = -2.0*x + 11.0;
  }

  void main ()      //主函数
  {
	  rat_max1  s(1.0, 1.0e-10);  //创建对象并顺序提供x0,eps值
	  s.solution ();     //执行一维极值连分式法
	  s.output ();   //输出标志与极值点以及函数值到文件并显示
  }

