  //7WITTY.CPP
  //Witty方法
  #include  <iostream>
  #include  <iomanip>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  witty
  {
  private: 
	       int n, k;
		   double  t, h, *y, **z;
		   double  *d, *q, *a;
  public:
	       witty (int nn, int kk)
		   {
			   int i;
			   n = nn; k = kk;
               y = new double[n];       //动态分配内存
	           d = new double[n];
	           q = new double[n];
	           a = new double[n];
	           z = new double*[n];
	           for (i=0; i<n; i++) { z[i] = new double[k]; }
		   }
	       void  input ();         //由文件读入数据成员t,h
		                      //以及n个未知函数在起始点t处的函数值
		   void solution ();    //执行Witty方法
           void output ();        //输出k个积分点上
		                          //的未知函数值到文件并显示
		   void func (double,double*,double*);   
		                         //计算微分方程组中各方程右端函数值
		   ~witty ()
		   { 
			   int i;
			   for (i=0; i<n; i++) { delete [] z[i]; }
			   delete [] z;
			   delete [] y, d, q, a;
		   }
  };

  void witty::input ()     //由文件读入数据成员t,h
		                      //以及n个未知函数在起始点t处的函数值
  {
	  int  i;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
      fin >>t >>h;        //读入t,h
      for (i=0; i<n; i++)  fin >>y[i];  
	                            //读入n个未知函数在起始点t处的函数值	  
	  fin.close ();
  }

  void witty::solution ()  //执行Witty方法
  { 
	  int i,j;
      double x;
//      for (i=0; i<=n-1; i++) z[i*k]=y[i];
      func (t,y,d);
      for (j=0; j<k; j++)
      { 
		  for (i=0; i<=n-1; i++)
              q[i]=y[i]+h*d[i]/2.0;
          x=t+(j-0.5)*h;
          func (x,q,a);
          for (i=0; i<=n-1; i++)
          { 
			  d[i]=2.0*a[i]-d[i];
              z[i][j]=y[i]+h*a[i];
			  y[i] = z[i][j];
          }
      }
  }

  void witty::output ()   //输出k个积分点上
	                      //的未知函数值到文件并显示
  {
	  int i, j;
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
      for (i=0; i<k; i++)
	  {
		  cout <<"t = " <<t+(i+1)*h <<endl;
		  for (j=0; j<n; j++)
		  {
			  fout <<z[j][i] <<"   ";
			  cout <<"y(" <<j <<")=" <<setw(12) <<z[j][i] <<"  ";
		  }
		  fout <<endl;  cout <<endl;
	  }
	  fout.close ();
  }

  void witty::func (double t, double y[], double d[])
  {
	  d[0] = y[1];
	  d[1] = -y[0];
	  d[2] = -y[2];
  }

  void main ()      //主函数
  {
	  witty  s(3, 10);           //创建对象
	  s.input ();        //由文件读入数据成员t,h
		                 //以及n个未知函数在起始点t处的函数值
	  s.solution ();    //执行Witty方法
	  s.output (); //输出k个积分点上的未知函数值到文件并显示
  }

