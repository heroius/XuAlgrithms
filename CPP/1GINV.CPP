  //1GINV.CPP
  //用奇异值分解法求矩阵广义逆
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  ginv
  {
  private: 
           int m, n, ka;
		   double  **a, **u, **v, **aa, eps;
		   double  *s, *e, *w, fg[2], cs[2];
  public:
	       ginv (int mm, int nn)
		   {
			   int i;
			   m = mm; n= nn;
	           a = new double*[m];   //动态分配内存空间
	           for (i=0; i<m; i++) a[i] = new double[n];
	           u = new double*[m];
	           for (i=0; i<m; i++) u[i] = new double[m];
	           v = new double*[n];
	           for (i=0; i<n; i++) v[i] = new double[n];
	           aa = new double*[n];
	           for (i=0; i<n; i++) aa[i] = new double[m];
               ka = m + 1;
	           if (m < n)  ka = n + 1;
               s = new double[ka];
               e = new double[ka];
               w = new double[ka];
		   }
	       void input ();        //由文件读入矩阵A的元素以及eps
		   void uav ();          //奇异值分解
           void ppp();            //奇异值分解的子程序
           void sss();            //奇异值分解的子程序
		   void a_ginv ();         //求A的广义逆
		   void output ();        //输出U矩阵和V矩阵以及奇异值
		                          //与A的广义逆到文件并显示
		   ~ginv ()
		   {
			   int i;
			   for (i=0; i<m; i++) { delete [] a[i]; }
			   delete [] a;
			   for (i=0; i<m; i++) { delete [] u[i]; }
			   delete [] u;
			   for (i=0; i<n; i++) { delete [] v[i]; }
			   delete [] v;
			   for (i=0; i<n; i++) { delete [] aa[i]; }
			   delete [] aa;
			   delete [] s, e, w;
		   }
  };

  void ginv::input ()      //由文件读入矩阵A的元素以及eps
  {
	  int  i, j;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
	  for (i=0; i<m; i++)                       //读入矩阵A
		  for (j=0; j<n; j++)  fin >>a[i][j];
	  fin >>eps;                                //读入eps
	  for (i=0; i<m; i++)
		  for (j=0; j<m; j++)  u[i][j] = 0;
	  for (i=0; i<n; i++)
		  for (j=0; j<n; j++)  v[i][j] = 0;
	  fin.close ();
  }

  void ginv::uav()   //奇异值分解
  { 
	  int i,j,k,l,it,ll,kk,mm,nn,m1,ks;
      double d,dd,t,sm,sm1,em1,sk,ek,b,c,shh;
      it=60; k=n;
      if (m-1<n) k=m-1;
      l=m;
      if (n-2<m) l=n-2;
      if (l<0) l=0;
      ll=k;
      if (l>k) ll=l;
      if (ll>=1)
      {
		  for (kk=1; kk<=ll; kk++)
          { 
			  if (kk<=k)
              { 
				  d=0.0;
                  for (i=kk; i<=m; i++)
                      d=d+a[i-1][kk-1]*a[i-1][kk-1];
                  s[kk-1]=sqrt(d);
                  if (s[kk-1]!=0.0)
                  { 
                      if (a[kk-1][kk-1]!=0.0)
                      {
						  s[kk-1]=fabs(s[kk-1]);
                          if (a[kk-1][kk-1]<0.0)  s[kk-1]=-s[kk-1];
                      }
                      for (i=kk; i<=m; i++)
                      { 
                          a[i-1][kk-1]=a[i-1][kk-1]/s[kk-1];
                      }
                      a[kk-1][kk-1]=1.0+a[kk-1][kk-1];
                  }
                  s[kk-1]=-s[kk-1];
              }
              if (n>=kk+1)
              { 
				  for (j=kk+1; j<=n; j++)
                  { 
					  if ((kk<=k)&&(s[kk-1]!=0.0))
                      { 
						  d=0.0;
                          for (i=kk; i<=m; i++)
                          { 
                              d=d+a[i-1][kk-1]*a[i-1][j-1];
                          }
                          d=-d/a[kk-1][kk-1];
                          for (i=kk; i<=m; i++)
                          { 
                              a[i-1][j-1]=a[i-1][j-1]+d*a[i-1][kk-1];
                          }
                      }
                    e[j-1]=a[kk-1][j-1];
                  }
              }
              if (kk<=k)
              { 
				  for (i=kk; i<=m; i++)
                  { 
                      u[i-1][kk-1]=a[i-1][kk-1];
                  }
              }
              if (kk<=l)
              { 
				  d=0.0;
                  for (i=kk+1; i<=n; i++)
                      d=d+e[i-1]*e[i-1];
                  e[kk-1]=sqrt(d);
                  if (e[kk-1]!=0.0)
                  { 
					  if (e[kk]!=0.0)
                      { 
						  e[kk-1]=fabs(e[kk-1]);
                          if (e[kk]<0.0) e[kk-1]=-e[kk-1];
                      }
                      for (i=kk+1; i<=n; i++)
                          e[i-1]=e[i-1]/e[kk-1];
                      e[kk]=1.0+e[kk];
                  }
                  e[kk-1]=-e[kk-1];
                  if ((kk+1<=m)&&(e[kk-1]!=0.0))
                  { 
					  for (i=kk+1; i<=m; i++) w[i-1]=0.0;
                      for (j=kk+1; j<=n; j++)
                      for (i=kk+1; i<=m; i++)
                           w[i-1]=w[i-1]+e[j-1]*a[i-1][j-1];
                      for (j=kk+1; j<=n; j++)
                      for (i=kk+1; i<=m; i++)
                      { 
                          a[i-1][j-1]=a[i-1][j-1]-w[i-1]*e[j-1]/e[kk];
                      }
                  }
                  for (i=kk+1; i<=n; i++)
                       v[i-1][kk-1]=e[i-1];
              }
          }
      }
      mm=n;
      if (m+1<n) mm=m+1;
      if (k<n) s[k]=a[k][k];
      if (m<mm) s[mm-1]=0.0;
      if (l+1<mm) e[l]=a[l][mm-1];
      e[mm-1]=0.0;
      nn=m;
      if (m>n) nn=n;
      if (nn>=k+1)
      { 
		  for (j=k+1; j<=nn; j++)
          { 
			  for (i=1; i<=m; i++)   u[i-1][j-1]=0.0;
              u[j-1][j-1]=1.0;
          }
      }
      if (k>=1)
      { 
		  for (ll=1; ll<=k; ll++)
          { 
			  kk=k-ll+1; 
              if (s[kk-1]!=0.0)
              { 
				  if (nn>=kk+1)
                    for (j=kk+1; j<=nn; j++)
                    { 
						d=0.0;
                        for (i=kk; i<=m; i++)
                        { 
                            d=d+u[i-1][kk-1]*u[i-1][j-1]/u[kk-1][kk-1];
                        }
                        d=-d;
                        for (i=kk; i<=m; i++)
                        {
                            u[i-1][j-1]=u[i-1][j-1]+d*u[i-1][kk-1];
                        }
                    }
                    for (i=kk; i<=m; i++)
                        u[i-1][kk-1]=-u[i-1][kk-1];
                    u[kk-1][kk-1]=1.0+u[kk-1][kk-1];
                    if (kk-1>=1)
                      for (i=1; i<=kk-1; i++)  u[i-1][kk-1]=0.0;
              }
              else
              { 
				  for (i=1; i<=m; i++)  u[i-1][kk-1]=0.0;
                  u[kk-1][kk-1]=1.0;
              }
          }
      }
      for (ll=1; ll<=n; ll++)
      { 
		  kk=n-ll+1;
          if ((kk<=l)&&(e[kk-1]!=0.0))
          { 
			  for (j=kk+1; j<=n; j++)
              { 
				  d=0.0;
                  for (i=kk+1; i<=n; i++)
                  { 
                      d=d+v[i-1][kk-1]*v[i-1][j-1]/v[kk][kk-1];
                  }
                  d=-d;
                  for (i=kk+1; i<=n; i++)
                  { 
                      v[i-1][j-1]=v[i-1][j-1]+d*v[i-1][kk-1];
                  }
              }
          }
          for (i=1; i<=n; i++)  v[i-1][kk-1]=0.0;
          v[kk-1][kk-1]=1.0;
      }
      for (i=1; i<=m; i++)
      for (j=1; j<=n; j++)  a[i-1][j-1]=0.0;
      m1=mm; it=60;
      while (1==1)
      { 
		  if (mm==0)
          { 
			  ppp();
			  return;
          }
          if (it==0)
          { 
			  ppp();
              cout <<"\n程序工作失败！" <<endl;  return;
          }
          kk=mm-1;
	      while ((kk!=0)&&(fabs(e[kk-1])!=0.0))
          { 
			  d=fabs(s[kk-1])+fabs(s[kk]);
              dd=fabs(e[kk-1]);
              if (dd>eps*d) kk=kk-1;
              else e[kk-1]=0.0;
          }
          if (kk==mm-1)
          { 
			  kk=kk+1;
              if (s[kk-1]<0.0)
              { 
				  s[kk-1]=-s[kk-1];
                  for (i=1; i<=n; i++)
                     v[i-1][kk-1]=-v[i-1][kk-1];
              }
              while ((kk!=m1)&&(s[kk-1]<s[kk]))
              { 
				  d=s[kk-1]; s[kk-1]=s[kk]; s[kk]=d;
                  if (kk<n)
                    for (i=1; i<=n; i++)
                    { 
                        d=v[i-1][kk-1]; 
						v[i-1][kk-1]=v[i-1][kk]; 
						v[i-1][kk]=d;
                    }
                  if (kk<m)
                    for (i=1; i<=m; i++)
                    { 
                        d=u[i-1][kk-1]; 
						u[i-1][kk-1]=u[i-1][kk]; 
						u[i-1][kk]=d;
                    }
                  kk=kk+1;
              }
              it=60;
              mm=mm-1;
          }
          else
          { 
			  ks=mm;
              while ((ks>kk)&&(fabs(s[ks-1])!=0.0))
              { 
				  d=0.0;
                  if (ks!=mm) d=d+fabs(e[ks-1]);
                  if (ks!=kk+1) d=d+fabs(e[ks-2]);
                  dd=fabs(s[ks-1]);
                  if (dd>eps*d) ks=ks-1;
                  else s[ks-1]=0.0;
              }
              if (ks==kk)
              { 
				  kk=kk+1;
                  d=fabs(s[mm-1]);
                  t=fabs(s[mm-2]);
                  if (t>d) d=t;
                  t=fabs(e[mm-2]);
                  if (t>d) d=t;
                  t=fabs(s[kk-1]);
                  if (t>d) d=t;
                  t=fabs(e[kk-1]);
                  if (t>d) d=t;
                  sm=s[mm-1]/d; sm1=s[mm-2]/d;
                  em1=e[mm-2]/d;
                  sk=s[kk-1]/d; ek=e[kk-1]/d;
                  b=((sm1+sm)*(sm1-sm)+em1*em1)/2.0;
                  c=sm*em1; c=c*c; shh=0.0;
                  if ((b!=0.0)||(c!=0.0))
                  { 
					  shh=sqrt(b*b+c);
                      if (b<0.0) shh=-shh;
                      shh=c/(b+shh);
                  }
                  fg[0]=(sk+sm)*(sk-sm)-shh;
                  fg[1]=sk*ek;
                  for (i=kk; i<=mm-1; i++)
                  { 
					  sss();
                      if (i!=kk) e[i-2]=fg[0];
                      fg[0]=cs[0]*s[i-1]+cs[1]*e[i-1];
                      e[i-1]=cs[0]*e[i-1]-cs[1]*s[i-1];
                      fg[1]=cs[1]*s[i];
                      s[i]=cs[0]*s[i];
                      if ((cs[0]!=1.0)||(cs[1]!=0.0))
                        for (j=1; j<=n; j++)
                        { 
                            d=cs[0]*v[j-1][i-1]+cs[1]*v[j-1][i];
                            v[j-1][i]=-cs[1]*v[j-1][i-1]+cs[0]*v[j-1][i];
                            v[j-1][i-1]=d;
                        }
                      sss ();
                      s[i-1]=fg[0];
                      fg[0]=cs[0]*e[i-1]+cs[1]*s[i];
                      s[i]=-cs[1]*e[i-1]+cs[0]*s[i];
                      fg[1]=cs[1]*e[i];
                      e[i]=cs[0]*e[i];
                      if (i<m)
                        if ((cs[0]!=1.0)||(cs[1]!=0.0))
                          for (j=1; j<=m; j++)
                          { 
                              d=cs[0]*u[j-1][i-1]+cs[1]*u[j-1][i];
                              u[j-1][i]=-cs[1]*u[j-1][i-1]+cs[0]*u[j-1][i];
                              u[j-1][i-1]=d;
                          }
                  }
                  e[mm-2]=fg[0];
                  it=it-1;
              }
              else
              { 
				  if (ks==mm)
                  { 
					  kk=kk+1;
                      fg[1]=e[mm-2]; e[mm-2]=0.0;
                      for (ll=kk; ll<=mm-1; ll++)
                      { 
						  i=mm+kk-ll-1;
                          fg[0]=s[i-1];
                          sss ();
                          s[i-1]=fg[0];
                          if (i!=kk)
                          { 
							  fg[1]=-cs[1]*e[i-2];
                              e[i-2]=cs[0]*e[i-2];
                          }
                          if ((cs[0]!=1.0)||(cs[1]!=0.0))
                            for (j=1; j<=n; j++)
                            { 
                                d=cs[0]*v[j-1][i-1]+cs[1]*v[j-1][mm-1];
                                v[j-1][mm-1]=-cs[1]*v[j-1][i-1]+cs[0]*v[j-1][mm-1];
                                v[j-1][i-1]=d;
                            }
                      }
                  }
                  else
                  { 
					  kk=ks+1;
                      fg[1]=e[kk-2];
                      e[kk-2]=0.0;
                      for (i=kk; i<=mm; i++)
                      { 
						  fg[0]=s[i-1];
                          sss ();
                          s[i-1]=fg[0];
                          fg[1]=-cs[1]*e[i-1];
                          e[i-1]=cs[0]*e[i-1];
                          if ((cs[0]!=1.0)||(cs[1]!=0.0))
                            for (j=1; j<=m; j++)
                            { 
                                d=cs[0]*u[j-1][i-1]+cs[1]*u[j-1][kk-2];
                                u[j-1][kk-2]=-cs[1]*u[j-1][i-1]+cs[0]*u[j-1][kk-2];
                                u[j-1][i-1]=d;
                            }
                      }
                  }
              }
          }
      }
  }

  void ginv::ppp ()
  { 
	  int i,j;
      double d;
      if (m>=n) i=n;
      else i=m;
      for (j=1; j<=i-1; j++)
      { 
		  a[j-1][j-1]=s[j-1];
          a[j-1][j]=e[j-1];
      }
      a[i-1][i-1]=s[i-1];
      if (m<n) a[i-1][i]=e[i-1];
      for (i=1; i<=n-1; i++)
      for (j=i+1; j<=n; j++)
      { 
          d=v[i-1][j-1]; 
		  v[i-1][j-1]=v[j-1][i-1]; 
		  v[j-1][i-1]=d;
      }
  }

  void ginv::sss ()
  { 
	  double r,d;
      if ((fabs(fg[0])+fabs(fg[1]))==0.0)
	  { 
		  cs[0]=1.0; cs[1]=0.0; d=0.0;
	  }
      else 
      { 
		  d=sqrt(fg[0]*fg[0]+fg[1]*fg[1]);
          if (fabs(fg[0])>fabs(fg[1]))
          { 
			  d=fabs(d);
              if (fg[0]<0.0) d=-d;
          }
          if (fabs(fg[1])>=fabs(fg[0]))
          { 
			  d=fabs(d);
              if (fg[1]<0.0) d=-d;
          }
          cs[0]=fg[0]/d; cs[1]=fg[1]/d;
      }
      r=1.0;
      if (fabs(fg[0])>fabs(fg[1])) r=cs[1];
      else
          if (cs[0]!=0.0) r=1.0/cs[0];
      fg[0]=d; fg[1]=r;
  }

  void ginv::a_ginv ()         //求A的广义逆
  { 
	  int i,j,k,l;
	  uav ();           //奇异值分解
      j=n;
      if (m<n) j=m;
      j=j-1;
      k=0;
      while ((k<=j)&&(a[k][k]!=0.0)) k=k+1;
      k=k-1;
      for (i=0; i<=n-1; i++)
      for (j=0; j<=m-1; j++)
      { 
		  aa[i][j]=0.0;
          for (l=0; l<=k; l++)
          { 
              aa[i][j]=aa[i][j]+v[l][i]*u[j][l]/a[l][l];
          }
      }
  }

  void ginv::output () //输出U矩阵和V矩阵以及奇异值
	                   //与广义逆到文件并显示
  {
	  int  i, j;
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  for (i=0; i<m; i++)
	  {
		  for (j=0; j<m; j++)
		  {
			  fout <<"    " <<u[i][j];
			  cout <<"    " <<u[i][j];
		  }
		  fout <<endl;  cout <<endl;
	  }
	  fout <<endl;  cout <<endl;
	  for (i=0; i<n; i++)
	  {
		  for (j=0; j<n; j++)
		  {
			  fout <<"    " <<v[i][j];
			  cout <<"    " <<v[i][j];
		  }
		  fout <<endl;  cout <<endl;
	  }
	  fout <<endl;  
	  cout <<endl;
	  for (i=0; i<m; i++)
	  {
		  for (j=0; j<n; j++)
		  {
			  fout <<"    " <<a[i][j];
			  cout <<"    " <<a[i][j];
		  }
		  fout <<endl;  cout <<endl;
	  }	  
	  fout <<endl;
	  cout <<endl;
	  for (i=0; i<n; i++)
	  {
		  for (j=0; j<m; j++)
		  {
			  fout <<"    " <<aa[i][j];
			  cout <<"    " <<aa[i][j];
		  }
		  fout <<endl;  cout <<endl;
	  }
	  fout.close ();
  }

  void main ()      //主函数
  {
	  ginv  uv(5, 4);
	  uv.input ();         //由文件读入矩阵A的元素以及eps
	  uv.a_ginv ();           //计算广义逆
	  uv.output ();        //输出U矩阵和V矩阵以及奇异值
	                       //与广义逆到文件并显示
  }
