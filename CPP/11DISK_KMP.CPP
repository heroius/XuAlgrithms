  //11DISK_KMP.CPP
  //字符串匹配
  #include  <iostream>
  #include  <iomanip>
  #include  <fstream>
  #include  <cmath>
  #include  <string>
  using namespace std;
  class  disk_kmp
  {
  private: 
	       
  public:
		   int kmp (char *, char *);       //执行KMP算法
  };

  int disk_kmp::kmp(char *p, char *str)  //执行KMP算法
  { 
	  int i,j,k,jt,*flink;
      char s;
      k=strlen(p);
      flink = new int[k];

      flink[0]=-1; i=1;
      while (i<=k-1)   //构造失败链接数组
      { 
		  j=flink[i-1];
          while ((j!=-1)&&(p[j]!=p[i-1])) j=flink[j];
          flink[i]=j+1;
          i=i+1;
      }      

	  ifstream fin (str);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str <<endl; exit(1); }	  
	  i=1; j=0;
      jt=1;
      if ((s=fin.get())==EOF) { i = 0; jt = 0; }
      while (jt==1)
      { 
		  while ((j!=-1)&&(p[j]!=s)) j=flink[j];
          if (j!=k-1) 
          { 
			  j=j+1; i=i+1;
              if ((s=fin.get())==EOF) { jt=0; i=0;}
          }
          else  { i = i - k + 1; jt = 0; }
      }
	  fin.close ();
      delete [] flink; 
	  return i;
  }

  void main ()      //主函数
  {
	  disk_kmp p;               //创建对象
	  cout <<p.kmp("gpq", "disk_kmp.txt") <<endl;

  }
