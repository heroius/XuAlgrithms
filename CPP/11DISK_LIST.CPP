  //11DISK_LIST.CPP
  //磁盘文本文件的查找与排序
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
  template <typename T>
  class  disk_list
  {
  private: 
	       long int n;    //记录个数
  public:
	       disk_list(long int nn)  {  n = nn; }
		   void  quick ();       //快速排序
           void  disk(fstream&, long int, long int);//递归快速排序		   
		   void  ibsh (char *);        //对分查找
		   void  serch (char*);       //顺序查找
           void  output ();       //显示文件内容
		   void  PRINT(T);        //显示一个记录
  };

  template <typename T>
  void  disk_list<T>::quick ()       //快速排序
  {
	  char str2[20];
	  cout <<"\n输入排序文件名:  ";
	  cin >>str2;
	  fstream f (str2);
	  if (!f)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
      disk(f, 0, n-1);
  }

  template <typename T>
  void disk_list<T>::disk(fstream& f, long int k, long int m)
  { 
	  long int i,j,kk;
	  T w, u, t;
      i=k; j=m; kk=(long)(i+j)/2;
      f.seekg(kk*sizeof(T));
      f.read((char *)&w, sizeof(T));
      t = w;
      do
	  { 
		  f.seekg(i*sizeof(T));
		  f.read((char *)&w, sizeof(T));		
          while ((strcmp(w.SORT_KEY,t.SORT_KEY)<0)&&(i<=j))
		  { 
			  i=i+1; 
			  f.seekg(i*sizeof(T));
		      f.read((char *)&w, sizeof(T));
		  }
		  f.seekg(j*sizeof(T));
		  f.read((char *)&w, sizeof(T));
          while ((strcmp(w.SORT_KEY,t.SORT_KEY)>0)&&(i<=j))
		  { 
			  j=j-1; 
		      f.seekg(j*sizeof(T));
		      f.read((char *)&w, sizeof(T));
		  }
          if (i<=j)
		  {
		      f.seekg(i*sizeof(T));
		      f.read((char *)&w, sizeof(T));
		      f.seekg(j*sizeof(T));
		      f.read((char *)&u, sizeof(T));

		      f.seekg(i*sizeof(T));
		      f.write((char *)&u, sizeof(T));
		      f.seekg(j*sizeof(T));
		      f.write((char *)&w, sizeof(T));
			 
			  i=i+1; j=j-1;
		  }
	  }while (i<=j);
      if (k<j) disk(f,k,j);
      if (i<m) disk(f,i,m);
      return;
  }

  template <typename T>
  void disk_list<T>::output () //显示文件内容
  {
	  int k;
	  T  x;
	  char str2[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str2;
	  ifstream fin (str2);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  for (k=0; k<n; k++)
	  {
		  fin.read((char *)&x, sizeof(x));
		  PRINT (x);
	  }
	  fin.close ();
  }

  template <typename T>
  void disk_list<T>::ibsh (char *x)       //对分查找
  {
	  T w;
	  char str2[20];
	  int i,j,k;
	  cout <<"\n输入对分查找文件名:  ";
	  cin >>str2;
	  ifstream fin (str2);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }	  
	  i=1; j=n;
      while (i<=j)
      { 
		  k=(i+j)/2;
		  fin.seekg((k-1)*sizeof(T));
		  fin.read((char *)&w, sizeof(T));

          if (strcmp(w.SERCH_KEY,x) == 0)
		  {
			  j = k - 1;
			  while ((j<n)&&(strcmp(w.SERCH_KEY,x) == 0))
			  {
                  PRINT (w);      //输出一个记录
				  j = j + 1;
		          fin.read((char *)&w, sizeof(T));
			  }
	          fin.close ();
			  return;
		  }
          if (strcmp(w.SERCH_KEY,x) > 0) j=k-1;
          else i=k+1;
      }
      cout <<"该磁盘文件中没有这个记录！" <<endl;
  }

  template <typename T>
  void disk_list<T>::serch (char *key)   //顺序查找
  {
	  int k;
	  T  x;
	  char str2[20];
	  cout <<"\n输入顺序查找文件名:  ";
	  cin >>str2;
	  ifstream fin (str2);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  for (k=0; k<n; k++)
	  {
		  fin.read((char *)&x, sizeof(x));
		  if (strcmp(x.SERCH_KEY,key) == 0)  PRINT (x);
	  }
	  fin.close ();
  }

  template <typename T>
  void  disk_list<T>::PRINT (T x)   //显示一个记录
  {
	  cout <<x.num;
	  cout <<x.name;
	  cout <<x.sex;
	  cout <<x.age;
	  cout <<x.score <<endl;
  }

  struct Student
  {
	  char  num[6];
	  char  name[9];
	  char  sex[4];
	  char  age[5];
	  char  score[8];
  };
  void main ()      //主函数
  {
/*
	  int i;
	  struct Student x[10] ={ 
	  {"101  ", "Zhang   ", "M  ", "19  ", "095.6  "},
	  {"102  ", "Wang    ", "F  ", "18  ", "092.4  "},
	  {"103  ", "Zhao    ", "M  ", "19  ", "085.7  "},
	  {"104  ", "Li      ", "M  ", "20  ", "096.3  "},
	  {"105  ", "Gou     ", "M  ", "19  ", "090.2  "},
	  {"106  ", "Lin     ", "M  ", "18  ", "091.5  "},
	  {"107  ", "Ma      ", "F  ", "17  ", "098.7  "},
	  {"108  ", "Zhen    ", "M  ", "21  ", "090.1  "},
	  {"109  ", "Xu      ", "M  ", "19  ", "100.0  "},
	  {"110  ", "Mao     ", "F  ", "18  ", "094.9  "}};
	  
	          ofstream fout ("disk_list.dat");
	          if (!fout)
			  { cout <<"\n不能打开这个文件 " <<endl; exit(1); }
      for (i=0; i<10; i++)
		  fout.write((char *)&x[i], sizeof(struct Student));
*/

	  disk_list<struct Student>  p(10);       //创建对象
	  p.output ();              //显示磁盘文件内容
	  p.serch ("18  ");         //顺序查找成员age为"18  "的记录
	  p.quick ();               //对磁盘文件按成员age进行快速排序
	  p.output ();              //显示排序后的磁盘文件内容
	  p.ibsh ("19  ");          //对分查找成员age为"19  "的记录
  }


