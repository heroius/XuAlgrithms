  //8RHIS.CPP
  //随机样本分析程序
  #include  <iostream>
  #include  <iomanip>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  rhis
  {
  private: 
	       int n, m, k, *g, *q;
		   double  *x, x0, h, dt[3];
  public:
	       rhis (int nn, int mm)
		   {
			   n = nn; m = mm;
               g = new int[m];   //动态分配内存
	           q = new int[m];
	           x = new double[n];
		   }
	       void  input ();         //由文件读入x0,h
		                      //以及n个随机变量的样本点值
		   void analysis ();    //执行分析计算
           void output ();        //输出直方图以及各参数到文件并显示
		   ~rhis ()
		   {  delete [] g, q, x;  }
  };

  void rhis::input () //由文件读入x0,h以及n个随机变量的样本点值
  {
	  int  i;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
      fin >>x0 >>h;        //读入x0,h
      for (i=0; i<n; i++)  fin >>x[i];  //读入n个随机变量的样本点值	  
	  fin.close ();
  }

  void rhis::analysis()  //执行分析计算
  { 
	  int i,j,max;
      double s;
      dt[0]=0.0;
      for (i=0; i<=n-1; i++) dt[0]=dt[0]+x[i]/n;
      dt[1]=0.0;
      for (i=0; i<=n-1; i++)
          dt[1]=dt[1]+(x[i]-dt[0])*(x[i]-dt[0]);
      dt[1]=dt[1]/n;
      dt[2]=sqrt(dt[1]);
      for (i=0; i<=m-1; i++)
      { 
		  q[i]=0;
          s=x0+(i+0.5)*h-dt[0];
          s=exp(-s*s/(2.0*dt[1]));
          g[i]=n*s*h/(dt[2]*2.5066);
      }
      s=x0+m*h;
      for (i=0; i<=n-1; i++)
          if ((x[i]-x0)>=0.0)
              if ((s-x[i])>=0.0)
			  { j=(x[i]-x0)/h; q[j]=q[j]+1; }
      k=1; max=0;
      for (i=0; i<=m-1; i++)
          if (q[i]>max) max=q[i];
      while (max>50) { max=max/2; k=2*k;}
  }

  void rhis::output ()   //输出直方图以及各参数到文件并显示
  {
	  int i, j, t;
      double s;
      char a[50];
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
      fout <<endl;
	  cout <<endl;
	  fout <<n <<endl <<endl;
	  cout <<"n = " <<n <<endl <<endl;
      fout <<x0 <<"   " <<h <<"   " <<m <<endl <<endl;
	  cout <<"x0 = " <<x0 <<"   " <<"h = " <<h <<"   " 
		   <<"m = " <<m <<endl <<endl;
      fout <<dt[0] <<endl <<dt[1] <<endl <<dt[2] <<endl <<endl;
	  cout <<"算术平均值 = " <<dt[0] <<endl
		   <<"方差 = "<<dt[1] <<endl 
		   <<"标准差 = " <<dt[2] <<endl <<endl;
      for (i=0; i<=m-1; i++)
      { 
		  s=x0+(i+0.5)*h;
          for (j=0; j<=49; j++) a[j]=' ';
          j=q[i]/k;
          for (t=0; t<=j-1; t++) a[t]='X';
          j=g[i]/k;
          if ((j>0)&&(j<=50)) a[j]='*';
		  fout <<setw(12) <<s <<setw(5) <<q[i] <<"   ";
		  cout <<setw(12) <<s <<setw(5) <<q[i] <<"   ";
          for (j=0; j<=49; j++)
		  { fout <<a[j];  cout <<a[j]; }
          fout <<endl;   cout <<endl;
      }
      fout <<endl;   cout <<endl;
	  fout <<"1 : " <<k <<endl;
      cout <<"1 : " <<k <<endl;
      fout <<endl;   cout <<endl;
	  fout.close ();
  }

  void main ()      //主函数
  {
	  rhis  s(100, 10);           //创建对象
	  s.input ();   //由文件读入x0,h以及n个随机变量的样本点值
		                 //以及n个未知函数在起始点t处的函数值
	  s.analysis ();    //执行分析计算
	  s.output (); //
  }

