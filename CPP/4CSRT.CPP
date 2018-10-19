  //4CSRT.CPP
  //牛顿下山法求复系数代数方程全部实根
  #include  <iostream>
  #include  <fstream>
  #include  <cmath>
  using namespace std;
  class  csrt
  {
  private: 
           int n;
		   double  *ar, *ai, *xr, *xi;
  public:
	       csrt (int nn)
		   {
			   n = nn;
	           ar = new double[n+1];   //动态分配内存
	           ai = new double[n+1];
	           xr = new double[n];
	           xi = new double[n];
		   }
	       void input ();          //由文件读入代数方程系数
		   void csrt_root ();      //执行牛顿下山法
           void output ();     //输出根到文件并显示
		   ~csrt ()
		   {  delete [] ar, ai, xr, xi;   }
  };

  void csrt::input ()      //由文件读入代数方程系数
  {
	  int  i;
	  char str1[20];
	  cout <<"\n输入文件名:  ";
	  cin >>str1;
	  ifstream  fin (str1);
	  if (!fin)
	  { cout <<"\n不能打开这个文件 " <<str1 <<endl; exit(1); }
	  for (i=n; i>=0; i--)  fin >>ar[i];    //读入代数方程系数
	  for (i=n; i>=0; i--)  fin >>ai[i];
	  fin.close ();
  }

  void csrt::csrt_root ()        //执行牛顿下山法
  { 
	  int m,i,jt,k,is,it;
      double t,x,y,x1,y1,dx,dy,p,q,w,dd,dc,c;
      double g,u,v,pq,g1,u1,v1;
      m=n;
      p=sqrt(ar[m]*ar[m]+ai[m]*ai[m]);
      while ((m>0)&&(p+1.0==1.0))
      {  
		  m=m-1;
          p=sqrt(ar[m]*ar[m]+ai[m]*ai[m]);
      }
      if (m<=0)
      { 
		  cout <<"\n程序工作失败！" <<endl;
		  return;
	  }
      for (i=0; i<=m; i++)
      { 
		  ar[i]=ar[i]/p; 
		  ai[i]=ai[i]/p;
	  }
      for (i=0; i<=m/2; i++)
      { 
		  w=ar[i]; ar[i]=ar[m-i]; ar[m-i]=w;
          w=ai[i]; ai[i]=ai[m-i]; ai[m-i]=w;
      }
      k=m; is=0; w=1.0;
      jt=1;
      while (jt==1)
      { 
		  pq=sqrt(ar[k]*ar[k]+ai[k]*ai[k]);
	      while (pq<1.0e-12)
          { 
			  xr[k-1]=0.0; xi[k-1]=0.0; k=k-1;
              if (k==1)
              { 
				  p=ar[0]*ar[0]+ai[0]*ai[0];
                  xr[0]=-w*(ar[0]*ar[1]+ai[0]*ai[1])/p;
                  xi[0]=w*(ar[1]*ai[0]-ar[0]*ai[1])/p;
                  return;
              }
              pq=sqrt(ar[k]*ar[k]+ai[k]*ai[k]);
          }
	      q=log(pq); q=q/(1.0*k); q=exp(q);
          p=q; w=w*p;
          for (i=1; i<=k; i++)
          { 
			  ar[i]=ar[i]/q; 
			  ai[i]=ai[i]/q; 
			  q=q*p;
		  }
          x=0.0001; x1=x; y=0.2; y1=y; dx=1.0;
          g=1.0e+37; 
  l40:
          u=ar[0]; v=ai[0];
          for (i=1; i<=k; i++)
          { 
			  p=u*x1; q=v*y1;
              pq=(u+v)*(x1+y1);
              u=p-q+ar[i]; v=pq-p-q+ai[i];
          }
          g1=u*u+v*v;
          if (g1>=g)
          { 
			  if (is!=0)
              { 
				  it=1;
	              if (it==0)
				  { 
		              is=1;
                      dd=sqrt(dx*dx+dy*dy);
                      if (dd>1.0) dd=1.0;
                      dc=6.28/(4.5*k); c=0.0;
				  }
                  while(1==1)
				  { 
		              c=c+dc;
                      dx=dd*cos(c); dy=dd*sin(c);
                      x1=x+dx; y1=y+dy;
                      if (c<=6.29) { it=0; break; }
                      dd=dd/1.67;
                      if (dd<=1.0e-07) { it=1; break; }
					  c=0.0;
				  }
                  if (it==0) goto l40;
              }
              else
              {
	              it=1;
                  while (it==1)
				  { 
		              t=t/1.67; it=0;
                      x1=x-t*dx;  y1=y-t*dy;
                      if (k>=30)
					  { 
			               p=sqrt(x1*x1+y1*y1);
                           q=exp(75.0/k);
                           if (p>=q) it=1;
					  }
				  }
                  if (t>=1.0e-03) goto l40;
                  if (g>1.0e-18)
                  { 
					  it=0;
	                  if (it==0)
					  { 
		                  is=1;
                          dd=sqrt(dx*dx+dy*dy);
                          if (dd>1.0) dd=1.0;
                          dc=6.28/(4.5*k); c=0.0;
					  }
                      while(1==1)
					  { 
		                  c=c+dc;
                          dx=dd*cos(c); dy=dd*sin(c);
                          x1=x+dx; y1=y+dy;
                          if (c<=6.29) { it=0; break; }
                          dd=dd/1.67;
                          if (dd<=1.0e-07) { it=1; break; }
						  c=0.0;
					  }
                      if (it==0) goto l40;
                  }
              }
              for (i=1; i<=k; i++)
			  { 
		           ar[i]=ar[i]+ar[i-1]*x-ai[i-1]*y;
                   ai[i]=ai[i]+ar[i-1]*y+ai[i-1]*x;
			  }
              xr[k-1]=x*w; xi[k-1]=y*w;
              k=k-1;
              if (k==1)
			  { 
		           p=ar[0]*ar[0]+ai[0]*ai[0];
                   xr[0]=-w*(ar[0]*ar[1]+ai[0]*ai[1])/p;
                   xi[0]=w*(ar[1]*ai[0]-ar[0]*ai[1])/p;
			  }
          }
          else
          { 
			  g=g1; x=x1; y=y1; is=0;
              if (g<=1.0e-22)
			  {
				  for (i=1; i<=k; i++)
				  { 
		               ar[i]=ar[i]+ar[i-1]*x-ai[i-1]*y;
                       ai[i]=ai[i]+ar[i-1]*y+ai[i-1]*x;
				  }
                  xr[k-1]=x*w; xi[k-1]=y*w;
                  k=k-1;
                  if (k==1)
				  { 
		               p=ar[0]*ar[0]+ai[0]*ai[0];
                       xr[0]=-w*(ar[0]*ar[1]+ai[0]*ai[1])/p;
                       xi[0]=w*(ar[1]*ai[0]-ar[0]*ai[1])/p;
				  }
			  }
              else
              { 
				  u1=k*ar[0]; v1=ai[0];
                  for (i=2; i<=k; i++)
                  { 
					  p=u1*x; q=v1*y; pq=(u1+v1)*(x+y);
                      u1=p-q+(k-i+1)*ar[i-1];
                      v1=pq-p-q+(k-i+1)*ai[i-1];
                  }
                  p=u1*u1+v1*v1;
                  if (p<=1.0e-20)
                  { 
					  it=0;
	                  if (it==0)
					  { 
		                  is=1;
                          dd=sqrt(dx*dx+dy*dy);
                          if (dd>1.0) dd=1.0;
                          dc=6.28/(4.5*k); c=0.0;
					  }
                      while(1==1)
					  { 
		                  c=c+dc;
                          dx=dd*cos(c); dy=dd*sin(c);
                          x1=x+dx; y1=y+dy;
                          if (c<=6.29) { it=0; break; }
                          dd=dd/1.67;
                          if (dd<=1.0e-07) { it=1; break; }
						  c=0.0;
					  }
                      if (it==0) goto l40;
                      for (i=1; i<=k; i++)
					  { 
		                  ar[i]=ar[i]+ar[i-1]*x-ai[i-1]*y;
                          ai[i]=ai[i]+ar[i-1]*y+ai[i-1]*x;
					  }
                      xr[k-1]=x*w; xi[k-1]=y*w;
                      k=k-1;
                      if (k==1)
					  { 
		                  p=ar[0]*ar[0]+ai[0]*ai[0];
                          xr[0]=-w*(ar[0]*ar[1]+ai[0]*ai[1])/p;
                          xi[0]=w*(ar[1]*ai[0]-ar[0]*ai[1])/p;
					  }
                  }
                  else
                  { 
					  dx=(u*u1+v*v1)/p;
                      dy=(u1*v-v1*u)/p;
                      t=1.0+4.0/k;
	                  it=1;
                      while (it==1)
					  { 
		                  t=t/1.67; it=0;
                          x1=x-t*dx;  y1=y-t*dy;
                          if (k>=30)
						  { 
			                  p=sqrt(x1*x1+y1*y1);
                              q=exp(75.0/k);
                              if (p>=q) it=1;
						  }
					  }
                      if (t>=1.0e-03) goto l40;
                      if (g>1.0e-18)
                      { 
						  it=0;
	                      if (it==0)
						  { 
		                      is=1;
                              dd=sqrt(dx*dx+dy*dy);
                              if (dd>1.0) dd=1.0;
                              dc=6.28/(4.5*k); c=0.0;
						  }
                          while(1==1)
						  { 
		                      c=c+dc;
                              dx=dd*cos(c); dy=dd*sin(c);
                              x1=x+dx; y1=y+dy;
                              if (c<=6.29) { it=0; break; }
                              dd=dd/1.67;
                              if (dd<=1.0e-07) { it=1; break; }
							  c=0.0;
						  }
                          if (it==0) goto l40;
                      }
                      for (i=1; i<=k; i++)
					  { 
		                  ar[i]=ar[i]+ar[i-1]*x-ai[i-1]*y;
                          ai[i]=ai[i]+ar[i-1]*y+ai[i-1]*x;
					  }
                      xr[k-1]=x*w; xi[k-1]=y*w;
                      k=k-1;
                      if (k==1)
					  { 
		                  p=ar[0]*ar[0]+ai[0]*ai[0];
                          xr[0]=-w*(ar[0]*ar[1]+ai[0]*ai[1])/p;
                          xi[0]=w*(ar[1]*ai[0]-ar[0]*ai[1])/p;
					  }
                  }
              }
          }
          if (k==1) jt=0;
          else jt=1;
      }
  }

  void csrt::output ()       //输出根到文件并显示
  {
	  int  k;
	  char str2[20];
	  cout <<"\n输出文件名:  ";
	  cin >>str2;
	  ofstream fout (str2);
	  if (!fout)
	  { cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	  for (k=0; k<n; k++)
	  {
		  fout <<xr[k] <<"    " <<xi[k] <<endl;
		  cout <<xr[k] <<"  +j  " <<xi[k] <<endl;
	  }
	  fout.close ();
  }

  void main ()      //主函数
  {
	  csrt  root(5); 
	  root.input ();      //由文件读入代数方程系数
	  root.csrt_root ();    //执行牛顿下山法
	  root.output ();     //输出根到文件并显示
  }


