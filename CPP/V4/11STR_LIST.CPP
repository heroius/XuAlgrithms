  //11STR_LIST.CPP
  //结构表的查找与排序
  #include  <iostream>
  #include  <iomanip>
  #include  <fstream>
  #include  <cmath>
  #include  <string>
//  #define  SORT_KEY  name       //定义排序用关键字成员
//  #define  SERCH_KEY  name       //定义查找用关键字成员
  #define  SORT_KEY  age       //定义排序用关键字成员
  #define  SERCH_KEY  age       //定义查找用关键字成员
  using namespace std;
  template <typename T, typename KT>
  class  str_list
  {
  private: 
	       int n;    //表长度
		   T  *L;    //结构表
  public:
	       str_list(int nn)
		   {  n = nn;   L = new T[n];  }
	       void  input ();      //由文件读入结构表元素
		   void  heap ();       //结构表堆排序
		   void  sift(int,int);  //堆排序中的调整
		   void ibsh (KT);        //结构表对分查找
		   void serch (KT);       //结构表顺序查找
           void  output ();     //输出结构表中的元素到文件并显示
		   void  PRINT(int,ofstream&); //输出结构表中的一个元素到文件并显示
		   ~str_list ()
		   {  delete [] L;  }
  };

  template <typename T, typename KT>
  void str_list<T,KT>::input ()     //由文件读入结构表元素
  {
	  int  i;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
      for (i=0; i<n; i++)   //读入结构表元素的各成员
	  {
		  fin >>L[i].num;
		  fin >>L[i].name;
		  fin >>L[i].sex;
		  fin >>L[i].age;
		  fin >>L[i].score;
	  }
	  fin.close ();
  }

  template <typename T, typename KT>  
  void  str_list<T,KT>::heap ()       //结构表堆排序
  {
	  int i,mm;
      T t;
      mm=n/2;
      for (i=mm-1; i>=0; i--)  sift (i,n-1);
      for (i=n-1; i>=1; i--)
      { 
		  t=L[0]; L[0]=L[i]; L[i]=t;
          sift (0,i-1);
      }
      return;
  }

  template <typename T, typename KT>
  void str_list<T,KT>::sift (int i, int n)  //堆排序中的调整
  { 
	  int j;
      T t;
      t=L[i]; j=2*(i+1)-1;
      while (j<=n)
      { if ((j<n)&&(L[j].SORT_KEY<L[j+1].SORT_KEY)) j=j+1;
        if (t.SORT_KEY<L[j].SORT_KEY)  
		{ L[i]=L[j]; i=j; j=2*(i+1)-1; }
        else j=n+1;
      }
      L[i]=t;
      return;
  }

  template <typename T, typename KT>
  void str_list<T,KT>::output () //输出结构表中的元素到文件并显示
  {
	  int k;
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  for (k=0; k<n; k++)  PRINT (k, fout); //输出结构表中一个元素
	  fout.close ();
  }

  template <typename T, typename KT>
  void str_list<T,KT>::ibsh (KT x)       //结构表对分查找
  {
	  char str2[20];
	  int i,j,k;
      i=1; j=n;
      while (i<=j)
      { 
		  k=(i+j)/2;
          if (L[k-1].SERCH_KEY == x)
		  {
			  cout <<"\n输出文件名:  ";
	          cin >>str2;
	          ofstream fout (str2);
	          if (!fout)
			  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
			  j = k - 1;
			  while ((j<n)&&(L[j].SERCH_KEY == x))
			  {
                  PRINT (j, fout);      //输出结构表中一个元素
				  j = j + 1;
			  }
	          fout.close ();
			  return;
		  }
          if (L[k-1].SERCH_KEY > x) j=k-1;
          else i=k+1;
      }
      cout <<"该结构中没有这个记录！" <<endl;
  }

  template <typename T, typename KT>
  void str_list<T,KT>::serch (KT x)       //结构表顺序查找
  {
	  int k;
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  for (k=0; k<n; k++)
	  {
		  if (L[k].SERCH_KEY == x) PRINT (k, fout);
	  }
	  fout.close ();	  
  }

  template <typename T, typename KT>
  void  str_list<T,KT>::PRINT (int k, ofstream& fout) 
                         //输出结构表中的一个元素到文件并显示
  {
	  fout <<setw(5) <<L[k].num;     //输出到已经打开的文件
	  fout <<setw(10) <<L[k].name;
	  fout <<setw(4) <<L[k].sex;
	  fout <<setw(5) <<L[k].age;
	  fout <<setw(8) <<L[k].score <<endl;
	  cout <<setw(5) <<L[k].num;     //显示
	  cout <<setw(10) <<L[k].name;
	  cout <<setw(4) <<L[k].sex;
	  cout <<setw(5) <<L[k].age;
	  cout <<setw(8) <<L[k].score <<endl;
	  return;
  }

  struct Student
  {
	  int   num;
	  string  name;
	  char  sex;
	  int   age;
	  double  score;
  };
  void main ()      //主函数
  {
//	  str_list<struct Student, string>  p(10);       //创建对象
//	  p.input ();
//	  p.serch ("Mao");
//	  p.heap ();
//	  p.output ();
//	  p.ibsh ("Zhao");
	  str_list<struct Student, int>  p(10);       //创建对象
	  p.input ();
	  p.serch (18);
	  p.heap ();
	  p.output ();
	  p.ibsh (19);
  }


