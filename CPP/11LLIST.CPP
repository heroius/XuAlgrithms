  //11LLIST.CPP
  //顺序表的查找与排序
  #include  <iostream>
  #include  <iomanip>
  #include  <fstream>
  #include  <cmath>
  #include  <string>
  using namespace std;
  template <typename T>
  class  llist
  {
  private: 
	       int n;    //表长度
		   T  *L;    //顺序表
  public:
	       llist(int nn)
		   {  n = nn;   L = new T[n];  }
	       void  input ();      //由文件读入元素
		   void  bubble ();     //冒泡排序
		   void  bub (T *, int);  //执行冒泡排序
		   void  quick ();      //快速排序
           void  qck (T *, int);  //快速排序递归
           void  split (T *, int, int*);  //快速排序分割
		   void  shell ();      //Shell排序
		   void  heap ();       //堆排序
		   void  sift (int,int); //堆排序中的调整
		   int ibsh (T);        //对分查找
           void  output ();     //输出顺序表中的元素到文件并显示
		   ~llist ()
		   {  delete [] L;  }
  };

  template <typename T>
  void llist<T>::input ()     //由文件读入元素
  {
	  int  i;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
      for (i=0; i<n; i++)  fin >>L[i];   //读入元素
	  fin.close ();
  }

  template <typename T>
  void  llist<T>::bubble ()    //冒泡排序
  { bub (L, n);  return;  }

  template <typename T>
  void  llist<T>::bub (T *p, int n)    //执行冒泡排序
  {
	  int m,k,j,i;
      T d;
      k=0; m=n-1;
      while (k<m)
      { 
		  j=m-1; m=0;
          for (i=k; i<=j; i++)
              if (p[i]>p[i+1])
			  { d=p[i]; p[i]=p[i+1]; p[i+1]=d; m=i; }
          j=k+1; k=0;
          for (i=m; i>=j; i--)
              if (p[i-1]>p[i])
			  { d=p[i]; p[i]=p[i-1]; p[i-1]=d; k=i; }
      }
      return;
  }

  template <typename T>
  void  llist<T>::quick ()      //快速排序
  { qck(L, n);  return; }

  template <typename T>  
  void  llist<T>::qck (T *p, int n) //快速排序递归
  {
	  int m,i0,*i;
      T *s;
      i=&i0;
      if (n>10)
      { 
		  split(p,n,i);
          m=i0;
          qck(p,m);
	      s=p+(i0+1);
          m=n-(i0+1);
          qck(s,m);
      }
      else bub(p,n);
      return;
  }

  template <typename T>
  void llist<T>::split (T *p, int n, int *m)  //快速排序分割
  { 
	  int i,j,k,l;
      T t;
      i=0; j=n-1;
      k=(i+j)/2;
      if ((p[i]>=p[j])&&(p[j]>=p[k])) l=j;
      else if ((p[i]>=p[k])&&(p[k]>=p[j])) l=k;
      else l=i;
      t=p[l]; p[l]=p[i];
      while (i!=j)
      { 
		  while ((i<j)&&(p[j]>=t)) j=j-1;
          if (i<j)
          { 
			  p[i]=p[j]; i=i+1;
              while ((i<j)&&(p[i]<=t)) i=i+1;
              if (i<j)
              { p[j]=p[i]; j=j-1; }
          }
      }
      p[i]=t; *m=i;
      return;
  }

  template <typename T>  
  void  llist<T>::shell ()     //Shell排序
  {
	  int k,j,i;
      T t;
      k=n/2;
      while (k>0)
      { 
		  for (j=k; j<=n-1; j++)
          { 
			  t=L[j]; i=j-k;
              while ((i>=0)&&(L[i]>t))
              { L[i+k]=L[i]; i=i-k; }
              L[i+k]=t;
          }
          k=k/2;
      }
      return;
  }

  template <typename T>  
  void  llist<T>::heap ()       //堆排序
  {
	  int i,mm;
      T t;
      mm=n/2;
      for (i=mm-1; i>=0; i--) sift(i,n-1);
      for (i=n-1; i>=1; i--)
      { 
		  t=L[0]; L[0]=L[i]; L[i]=t;
          sift(0,i-1);
      }
      return;
  }

  template <typename T>
  void llist<T>::sift(int i, int n)  //堆排序中的调整
  { 
	  int j;
      T t;
      t=L[i]; j=2*(i+1)-1;
      while (j<=n)
      { if ((j<n)&&(L[j]<L[j+1])) j=j+1;
        if (t<L[j])  { L[i]=L[j]; i=j; j=2*(i+1)-1;}
        else j=n+1;
      }
      L[i]=t;
      return;
  }

  template <typename T>
  void llist<T>::output () //输出顺序表中的元素到文件并显示
  {
	  int i, k;
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  k = 0;
	  while (k<n)
	  {
		  i=0;
		  while ((i < 8) && (k < n))
		  {
			  cout <<L[k] <<"    ";   fout <<L[k] <<"    ";
			  k = k + 1;   i = i + 1;
		  }
		  cout <<endl;   fout <<endl;
	  }
	  fout.close ();
  }

  template <typename T>
  int llist<T>::ibsh (T x)       //对分查找
  {
	  int i,j,k;
      i=1; j=n;
      while (i<=j)
      { 
		  k=(i+j)/2;
          if (L[k-1] == x)  return k-1;
          if (L[k-1] > x) j=k-1;
          else i=k+1;
      }
      return (-1);
  }

  void main ()      //主函数
  {
	  llist<int>  p(12);       //创建对象
      llist<int>  p1(100);
	  llist<double>  p2(100);
	  llist<char>    p3(22);
	  llist<string>  p4(10);
	  p.input ();              //读入有序表p123.txt
      p.output ();
      cout <<p.ibsh (19) <<endl;
	  cout <<p.ibsh (34) <<endl;
	  p1.input ();             //读入100个随机整数llist_i.txt
	  p2.input ();             //读入100个随机实数list_r.txt
	  p3.input ();             //读入22个字符llist_c.txt
	  p4.input ();             //读入10个字符串llist_h.txt
//冒泡排序
	  p1.bubble ();   //整数排序
	  p1.output ();
	  p2.bubble ();   //实数排序
	  p2.output ();
	  p3.bubble ();   //字符排序
	  p3.output ();
	  p4.bubble ();   //字符串排序
	  p4.output ();
//快速排序
	  p1.quick ();     //整数排序
	  p1.output ();
	  p2.quick ();     //实数排序
	  p2.output ();
	  p3.quick ();     //字符排序
	  p3.output ();
	  p4.quick ();     //字符串排序
	  p4.output ();
//Shell排序
	  p1.shell ();     //整数排序
	  p1.output ();
	  p2.shell ();     //实数排序
	  p2.output ();
	  p3.shell ();     //字符排序
	  p3.output ();
	  p4.shell ();     //字符串排序
	  p4.output ();
//堆排序
	  p1.heap ();       //整数排序
	  p1.output ();
	  p2.heap ();       //实数排序
	  p2.output ();
	  p3.heap ();       //字符排序
	  p3.output ();
	  p4.heap ();       //字符串排序
	  p4.output ();
  }


