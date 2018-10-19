  //9RAT_MAXN.CPP
  //n维极值连分式法
  #include  <iostream>
  #include  <iomanip>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  rat_maxn
  {
  private: 
	       int n, flag;
		   double  *x, z, eps;
  public:
	       rat_maxn (int nn)
		   {
			   n = nn;
               x = new double[n];    //动态分配内存
		   }
	       void input ();       //由文件读入eps以及初值
		   void solution ();    //执行n维极值连分式法
		   double func (int);   //计算目标函数值或偏导数值
           void output ();  //输出标志与极值点以及函数值到文件并显示
		   ~rat_maxn ()
		   {   delete [] x; }
  };

  void rat_maxn::input ()     //由文件读入eps以及初值
  {
	  int  i;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
      fin >>eps;                //读入eps
      for (i=0; i<n; i++)  fin >>x[i];    //读入n个初值
	  fin.close ();
  }

  void rat_maxn::solution ()  //执行n维极值连分式法
  { 
	  int i,j,m,l,jt,il,k;
      double y[10],b[10],p,zz,t,h1,h2,f,dx;
      k=0; jt=1; h2=0.0;
      while (jt==1)
      { 
		  t=0.0; k =k + 1;
          for (i=1; i<=n; i++)
          { 
			  f=func (i);
              t=t+fabs(f);
          }
          if (t<eps) jt=0;
          else
          { 
			  for (i=0; i<=n-1; i++)
              { 
				  il=5;
                  while (il!=0)
                  { 
					  j=0; t=x[i]; il=il-1;
                      while (j<=7)
                      { 
						  if (j<=2) zz=t+j*0.01;
                          else zz=h2;
                          x[i]=zz;
                          f=func (i+1);
                          if (fabs(f)+1.0==1.0) { j=10; il=0;}
                          else
                          { 
							  h1=f; h2=zz;
                              if (j==0) { y[0]=h1; b[0]=h2;}
                              else
                              { 
								  y[j]=h1; m=0; l=0;
                                  while ((m==0)&&(l<=j-1))
                                  { 
									  p=h2-b[l];
                                      if (fabs(p)+1.0==1.0) m=1;
                                      else h2=(h1-y[l])/p;
                                      l=l+1;
                                  }
                                  b[j]=h2;
                                  if (m!=0) b[j]=1.0e+35;
                                  h2=0.0;
                                  for (l=j-1; l>=0; l--)
                                      h2=-y[l]/(b[l+1]+h2);
                                  h2=h2+b[0];
                              }
                              j=j+1;
                          }
                      }
                      x[i]=h2; 
                  }
                  x[i]=zz;
              }
              if (k == 10) jt=0;
          }
      }
    flag = 1;
    z = func (0); 
    dx=0.00001; t=x[0];
    x[0]=t+dx; h1=func (0);
    x[0]=t-dx; h2=func (0);
    x[0]=t;
    t=h1+h2-2.0*f;
    if (t>0.0) flag = -1;
    j=1; jt=1;
    while (jt==1)
      { j=j+1; dx=0.00001; jt=0;
        t=x[j-1];
        x[j-1]=t+dx; h2=func (0);
        x[j-1]=t-dx; h1=func (0);
        x[j-1]=t; t=h1+h2-2.0*f;
        if ((t*flag<0.0)&&(j<n)) jt=1;
      }
    if (t*flag>0.0) flag = 0;
  }

  void rat_maxn::output () //输出标志与极值点以及函数值到文件并显示
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
	  if (flag == 1)  cout <<"MAX:  " <<endl;
	  else  cout <<"MIN:  " <<endl;
      for (i=0; i<n; i++)
	  {
		  cout <<"x(" <<i <<")=" <<x[i] <<endl;
		  fout <<x[i] <<endl;
	  }
	  cout <<"z = " <<z <<endl;
	  fout <<z <<endl;
	  fout.close ();
  }

  double rat_maxn::func (int k)  //计算目标函数值或各偏导数值
  {
	  double y;
	  switch(k)
	  {
		  case 0: y = (x[0]-1)*(x[0]-1)+(x[1]+2)*(x[1]+2)+2; break;
		  case 1: y = 2.0*(x[0]-1); break;
		  case 2: y = 2.0*(x[1]+2); break;
          default:;
	  }
	  return y;
  }

  void main ()      //主函数
  {
	  rat_maxn  s(2);      //创建对象
	  s.input ();       //由文件读入eps以及初值
	  s.solution ();     //执行n维极值连分式法
	  s.output ();   //输出标志与极值点以及函数值到文件并显示
  }
