using System;
using System.Collections.Generic;
using System.Text;

namespace Heroius.XuAlgrithms
{
    /// <summary>
    /// 数据处理
    /// </summary>
    public static class DataProcessing
    {
        /// <summary>
        /// 根据给定的一维随机样本， 要求：
        /// <para>1 计算算术平均值、方差与标准差。</para>
        /// <para>2 按高斯分布计算出在各给定区间上近似的理论样本点数.</para>
        /// <para>3 输出经验直方图。</para>
        /// </summary>
        /// <param name="x">x[n]存放随机变量的n个样本值</param>
        /// <param name="n">随机样本点数</param>
        /// <param name="x0">直方图中随机变扯的起始值</param>
        /// <param name="h">直方图中随机变址等区间长度值</param>
        /// <param name="m">直方图中区间总数</param>
        /// <param name="dt">dt[3]: 
        /// <para>dt(0)返回随机样本的算术平均值</para>
        /// <para>dt(1)返回随机样本的方差</para>
        /// <para>dt(2)返回随机样本的标准差</para>
        ///</param>
        /// <param name="g">返回m个区间的按高斯分布所应有的近似理论样本点数</param>
        /// <param name="q">返回落在m个区间中每一个区间上的随机样本实际点数</param>
        public static void RHIS(double[] x, int n, double x0, double h, int m, out double[] dt, out int[] g, out int[] q)
        {
            int j;
            double s;
            string a;
            dt = new double[3];
            g = new int[m];
            q = new int[m];
            dt[0] = 0;
            for (int i = 0; i < n; i++)
            {
                dt[0] = dt[0] + x[i] / n;
            }
            dt[1] = 0;
            for (int i = 0; i < n; i++)
            {
                dt[1] = dt[1] + (x[i] - dt[0]) * (x[i] - dt[0]);
            }
            dt[1] = dt[1] / n;
            dt[2] = Math.Sqrt(dt[1]);
            for (int i = 0; i < m; i++)
            {
                q[i] = 0;
                s = x0 + (i + 0.5) * h - dt[0];
                s = Math.Exp(-s * s / (2 * dt[1]));
                g[i] = Convert.ToInt32(Math.Floor(n * s * h / (dt[2] * 2.5066)));
            }
            s = x0 + m * h;
            for (int i = 0; i < n; i++)
            {
                if (x[i] >= x0)
                {
                    if (s >= x[i])
                    {
                        j = Convert.ToInt32(Math.Floor((x[i] - x0) / h));
                        q[j] = q[j] + 1;
                    }
                }
            }
        }

        /// <summary>
        /// 给定n个数据点(Xk,Yk)(k=0,1,...,n-1)，用直线y=ax+b进行回归分析
        /// </summary>
        /// <param name="x">存放自变量的n个取值</param>
        /// <param name="y">存放与自变量的n个取值相对应的随机变量y的观测值</param>
        /// <param name="n">观测点数</param>
        /// <param name="a">a[2]: 
        /// <para>a(0) 返回回归系数b</para><para>a(l) 返回回归系数a</para></param>
        /// <param name="dt">dt[6]:
        /// <para>dt(0) 返回偏差平方和q</para>
        /// <para>dt(1) 返回平均标准偏差s</para>
        /// <para>dt(2) 返回回归平方和p</para>
        /// <para>dt(3) 返回最大偏差Umax</para>
        /// <para>dt(4) 返回最小偏差Umin</para>
        /// <para>dt(5) 返回偏差平均值u</para>
        /// </param>
        public static void SQT1(double[] x, double[] y, int n, out double[] a, out double[] dt)
        {
            a = new double[2];
            double xx = 0, yy = 0, e = 0, f = 0, q = 0, u = 0, p = 0, umax = 0, umin = 1e30, s;
            for (int i = 0; i < n; i++)
            {
                xx += x[i] / n;
                yy += y[i] / n;
            }
            for (int i = 0; i < n; i++)
            {
                q = x[i] - xx;
                e += q * q;
                f += q * (y[i] - yy);
            }
            a[1] = f / e;
            a[0] = yy - a[1] * xx;
            for (int i = 0; i < n; i++)
            {
                s = a[1] * x[i] + a[0];
                q += (y[i] - s) * (y[i] - s);
                p += (s - yy) * (s - yy);
                e = Math.Abs(y[i] - s);
                if (e > umax) umax = e;
                if (e < umin) umin = e;
                u += e / n;
            }
            dt = new double[6] {
                q,
                Math.Sqrt(q/n),
                p,
                umax,
                umin,
                u
            };
        }

        /// <summary>
        /// 根据随机变量及自变量X0,X1,...,Xm-的n组观测值(X0k,X1k,...,Xm-1k,Yk)(k=0,1,...,n-1) 作线性回归分析。
        /// </summary>
        /// <param name="x">x[m,n]: 每一列存放m 个自变量的观测值</param>
        /// <param name="y">y[n]: 存放随机变量y的n个观测值</param>
        /// <param name="m">自变量个数</param>
        /// <param name="n">观测数据的组数</param>
        /// <param name="a">a[m+1]: 返回回归系数a0,a1,...,am-1,am</param>
        /// <param name="dt">dt[4]: 
        /// <para>dt(0)返回 偏差平方和q</para>
        /// <para>dt(1)返回 平均标准偏差s</para>
        /// <para>dt(2)返回 复相关系数r</para>
        /// <para>dt(3)返回 回归平方和u</para>
        /// </param>
        /// <param name="v">v[m]：返回m个自变量的偏相关系数</param>
        public static void SQT2(double[,] x, double[] y, int m, int n, out double[] a, out double[] dt, out double[] v)
        {
            double[] xx = Utility.C.Convert(x);
            double q, e, u, p, yy, s, r, pp;
            int mm = m + 1;
            double[] b = new double[mm * mm];
            b[mm * mm - 1] = n;
            a = new double[mm];
            v = new double[m];

            for (int j = 0; j < m; j++)
            {
                p = 0;
                for (int i = 0; i < n; i++)
                {
                    p += xx[j * n + i];
                }
                b[m * mm + j] = p;
                b[j * mm + m] = p;
            }
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    p = 0;
                    for (int k = 0; k < n; k++)
                    {
                        p += xx[i * n + k] * xx[j * n + k];
                    }
                    b[j * mm + i] = p;
                    b[i * mm + j] = p;
                }
            }
            a[m] = 0;
            for (int i = 0; i < n; i++)
            {
                a[m] += y[i];
            }
            for (int i = 0; i < m; i++)
            {
                a[i] = 0;
                for (int j = 0; j < n; j++)
                {
                    a[i] += xx[i * n + j] * y[j];
                }
            }

            var temp_b = Utility.C.Convert(b, mm, mm);
            var temp_a = Utility.C.Convert(a, mm, 1);
            LinearEquations.CHLK(ref temp_b, mm, 1, ref temp_a);
            b = Utility.C.Convert(temp_b);
            a = Utility.C.Convert(temp_a);

            yy = 0;
            for (int i = 0; i < n; i++)
            {
                yy += y[i] / n;
            }
            q = 0; e = 0; u = 0;
            for (int i = 0; i < n; i++)
            {
                p = a[m];
                for (int j = 0; j < m; j++)
                {
                    p += a[j] * xx[j * n + i];
                }
                q += (y[i] - p) * (y[i] - p);
                e += (y[i] - yy) * (y[i] - yy);
                u += (yy - p) * (yy - p);
            }
            s = Math.Sqrt(q / n);
            r = Math.Sqrt(1.0 - q / e);
            for (int j = 0; j < m; j++)
            {
                p = 0;
                for (int i = 0; i < n; i++)
                {
                    pp = a[m];
                    for (int k = 0; k < m; k++)
                    {
                        if (k != j) pp += a[k] * xx[k * n + i];
                    }
                    p += (y[i] - pp) * (y[i] - pp);
                }
                v[j] = Math.Sqrt(1.0 - q / p);
            }
            dt = new double[4] {
                q,
                s,
                r,
                u
            };
        }

        /// <summary>
        /// 对多元线性回归进行因子筛选，
        /// 最后给出一定显著性水平下各因子均为显著的回归方程中的诸回归系数、
        /// 偏回归平方和、估计的标准偏差、复相关系数以及F-检验值、
        /// 各回归系数的标准偏差、因变量条件期望值的估计值与残差。
        /// </summary>
        /// <param name="n">自变量x的个数</param>
        /// <param name="k">观测数据的点数</param>
        /// <param name="x">x[k,n+1]: 其中前n列存放自变量因子Xi(i=0,1,...,n-1)的k次观测值；最后一列存放因变量y的k次观测值</param>
        /// <param name="f1">欲选入因子时显著性检验的F-分布值</param>
        /// <param name="f2">欲剔除因子时显著性检验的F-分布值</param>
        /// <param name="eps">防止系数相关矩阵退化的判据</param>
        /// <param name="xx">xx[n+1]: 前n个分量返回n个自变量因子的算术平均值Xi(i=0,1,...,n-1)；最后一个分量返回因变量y的算术平均值y</param>
        /// <param name="b">b[n+1]: 返回回归方程中各因子的回归系数及常数项B0,B1,...,Bn</param>
        /// <param name="v">v[n+1]: 前n个分量返回各因子的偏回归平方和Vi(i=0,1,...,n-1)；最后一个分量返回残差平方和q</param>
        /// <param name="s">s[n+1]: 前n个分量返回各因子回归系数的标准偏差Si(i=0,1,...,n-1); 最后一个分量返回估计的标准偏差s</param>
        /// <param name="dt">dt[2]: dt(0)返回复相关系数；dt(1)返回F-检验值</param>
        /// <param name="ye">ye[k]: 返回对应于k个观测值的因变量条件期望值的k个估计值Ei(i=0,1,...,k-1)</param>
        /// <param name="yr">yr[k]: 返回因变量的k个观测值的残差δi(i=0,1,...,k-1)</param>
        /// <param name="r">r[n+1,n+1]: 返回最终的规格化的系数相关矩阵R</param>
        public static void SQT3(int n, int k, double[,] x, double f1, double f2, double eps, out double[] xx, out double[] b, out double[] v, out double[] s, out double[] dt, out double[] ye, out double[] yr, out double[,] r)
        {
            double[] temp_x = Utility.C.Convert(x);
            double[] temp_r = new double[(n + 1) * (n + 1)];
            xx = new double[n + 1];
            b = new double[n + 1];
            v = new double[n + 1];
            s = new double[n + 1];
            dt = new double[2];
            ye = new double[k];
            yr = new double[k];
            int m, imi, imx, l, it;
            double z, phi, sd, vmi, vmx, q, fmi, fmx;
            m = n + 1;
            q = 0;
            for (int j = 0; j <= n; j++)
            {
                z = 0;
                for (int i = 0; i < k; i++)
                {
                    z += temp_x[i * m + j] / k;
                }
                xx[j] = z;
            }
            for (int i = 0; i <= n; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    z = 0;
                    for (int ii = 0; ii < k; ii++)
                    {
                        z = z + (temp_x[ii * m + i] - xx[i]) * (temp_x[ii * m + j] - xx[j]);
                    }
                    temp_r[i * m + j] = z;
                }
            }
            for (int i = 0; i <= n; i++)
            {
                ye[i] = Math.Sqrt(temp_r[i * m + i]);
            }
            for (int i = 0; i <= n; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    temp_r[i * m + j] /= ye[i] * ye[j];
                    temp_r[j * m + i] = temp_r[i * m + j];
                }
            }
            phi = k - 1;
            sd = ye[n] / Math.Sqrt(k - 1);
            it = 1;
            while (it==1)
            {
                it = 0;
                vmi = 1e35;
                vmx = 0;
                imi = -1;
                imx = -1;
                for (int i = 0; i <= n; i++)
                {
                    v[i] = temp_r[i * m + n] * temp_r[n * m + i] / temp_r[i * m + i];
                    if (v[i]>=0)
                    {
                        if (v[i]>vmx)
                        {
                            vmx = v[i];
                            imx = i;
                        }
                    }
                    else
                    {
                        b[i] = temp_r[i * m + n] * ye[n] / ye[i];
                        s[i] = Math.Sqrt(temp_r[i * m + i]) * sd / ye[i];
                        if (Math.Abs(v[i]) < vmi)
                        {
                            vmi = Math.Abs(v[i]);
                            imi = i;
                        }
                    }
                }
                if (phi!=n-1)
                {
                    z = 0;
                    for (int i = 0; i < n; i++)
                    {
                        z = z + b[i] * xx[i];
                    }
                    b[n] = xx[n] - z;
                    s[n] = sd;
                    v[n] = q;
                }
                else
                {
                    b[n] = xx[n];
                    s[n] = sd;
                }
                fmi = vmi * phi / temp_r[n * m + n];
                fmx = (phi - 1) * vmx / (temp_r[n * m + n] - vmx);
                if ((fmi<f2)||(fmx>=f1))
                {
                    if (fmi<f2)
                    {
                        phi++;
                        l = imi;
                    }
                    else
                    {
                        phi--;
                        l = imx;
                    }
                    for (int i = 0; i <= n; i++)
                    {
                        if (i!=1)
                        {
                            for (int j = 0; j <= n; j++)
                            {
                                if (j!=l)
                                {
                                    temp_r[i * m + j] -= (temp_r[l * m + j] / temp_r[l * m + l]) * temp_r[i * m + l];
                                }
                            }
                        }
                    }
                    for (int j = 0; j <= n; j++)
                    {
                        if (j!=l)
                        {
                            temp_r[l * m + j] /= temp_r[l * m + l];
                        }
                    }
                    for (int i = 0; i <= n; i++)
                    {
                        if (i!=l)
                        {
                            temp_r[i * m + l] /= -temp_r[l * m + l];
                        }
                    }
                    temp_r[l * m + l] = 1 / temp_r[l * m + l];
                    q = temp_r[n * m + n] * ye[n] * ye[n];
                    sd = Math.Sqrt(temp_r[n * m + n] / phi) * ye[n];
                    dt = new double[]
                    {
                        Math.Sqrt(1-temp_r[n*m+n]),
                        (phi*(1-temp_r[n*m+n]))/((k-phi-1)*temp_r[n*m+n])
                    };
                    it = 1;
                }
            }
            for (int i = 0; i < k; i++)
            {
                z = 0;
                for (int j = 0; j < n; j++)
                {
                    z += b[j] * temp_x[i * m + j];
                }
                ye[i] = b[n] + z;
                yr[i] = temp_x[i * m + n] - ye[i];
            }

            r = Utility.C.Convert(temp_r, n + 1, n + 1);
        }

        /// <summary>
        /// 对于给定的n个数据点(Xi,Yi)(i=0,1,...,n-1)，用 y=bt^ax, t>0 作拟合。
        /// </summary>
        /// <param name="n">数据点数</param>
        /// <param name="x">x[n]: 存放n个数据点</param>
        /// <param name="y">y[n]: 存放n个数据点，要求所有的y>0</param>
        /// <param name="t">指数函数的底，要求t>0</param>
        /// <param name="a">a[7]: 返回拟合函数的参数以及各种统计量。
        /// <para>a(0): 拟合函数y=bt^ax中的 b</para>
        /// <para>a(1): 拟合函数y=bt^ax中的 a</para>
        /// <para>a(2): 偏差平方和q，即q=Σ(Yi-bt^aXi)^2</para>
        /// <para>a(3): 平均标准偏差s，即s=√(q/n)</para>
        /// <para>a(4): 最大偏差Umax，即Umax=max|Yi-bt^aXi|</para>
        /// <para>a(5): 最小偏差Umin，即Umin=min|Yi-bt^aXi|</para>
        /// <para>a(6): 偏差平均值u，即u=(1/n)Σ|Yi-bt^aXi|</para>
        /// </param>
        public static void LOG1(int n, double[] x, double[] y, double t, out double[] a)
        {
            double xx=0, yy=0, dx=0, dxy=0;
            a = new double[7];
            for (int i = 0; i < n; i++)
            {
                xx += x[i] / n;
                yy += Math.Log(y[i]) / Math.Log(t) / n;
            }
            for (int i = 0; i < n; i++)
            {
                a[2] = x[i] - xx;
                dx += a[2] * a[2];
                dxy += a[2] * (Math.Log(y[i]) / Math.Log(t) - yy);
            }
            a[1] = dxy / dx;
            a[0] = yy - a[1] * xx;
            a[0] = a[0] * Math.Log(t);
            a[0] = Math.Exp(a[0]);
            a[2] = 0;
            a[6] = 0;
            a[4] = 0;
            a[5] = 1e30;
            for (int i = 0; i < n; i++)
            {
                a[3] = a[1] * x[i] * Math.Log(t);
                a[3] = a[0] * Math.Exp(a[3]);
                a[2] = a[2] + (y[i] - a[3]) * (y[i] - a[3]);
                dx = Math.Abs(y[i] - a[3]);
                if (dx > a[4]) a[4] = dx;
                if (dx < a[5]) a[5] = dx;
                a[6] += dx / n;
            }
            a[3] = Math.Sqrt(a[2] / n);
        }

        /// <summary>
        /// 对于给定的n个数据点(Xk,Yk)(k=0,1,...,n-1)，用 y=bx^a, 作拟合。
        /// </summary>
        /// <param name="n">数据点数</param>
        /// <param name="x">x[n]: 存放n个数据点</param>
        /// <param name="y">y[n]: 存放n个数据点，要求所有的y>0</param>
        /// <param name="a">a[7]: 返回拟合函数的参数以及各种统计量。
        /// <para>a(0): 拟合函数y=bt^ax中的 b</para>
        /// <para>a(1): 拟合函数y=bt^ax中的 a</para>
        /// <para>a(2): 偏差平方和q，即q=Σ(Yi-bt^aXi)^2</para>
        /// <para>a(3): 平均标准偏差s，即s=√(q/n)</para>
        /// <para>a(4): 最大偏差Umax，即Umax=max|Yi-bt^aXi|</para>
        /// <para>a(5): 最小偏差Umin，即Umin=min|Yi-bt^aXi|</para>
        /// <para>a(6): 偏差平均值u，即u=(1/n)Σ|Yi-bt^aXi|</para>
        /// </param>
        public static void LOG2(int n, double[] x, double[] y, out double[] a)
        {
            double xx=0, yy=0, dx=0, dxy=0;
            a = new double[7];
            for (int i = 0; i < n; i++)
            {
                xx += Math.Log(x[i]) / n;
                yy += Math.Log(y[i]) / n;
            }
            for (int i = 0; i < n; i++)
            {
                a[2] = Math.Log(x[i]) - xx;
                dx += a[2] * a[2];
                dxy += a[2] * (Math.Log(y[i]) - yy);
            }
            a[1] = dxy / dx;
            a[0] = yy - a[1] * xx;
            a[0] = Math.Exp(a[0]);
            a[2] = 0;
            a[6] = 0;
            a[4] = 0;
            a[5] = 1e30;
            for (int i = 0; i < n; i++)
            {
                a[3] = a[1] * Math.Log(x[i]);
                a[3] = a[0] * Math.Exp(a[3]);
                a[2] += (y[i] - a[3]) * (y[i] - a[3]);
                dx = Math.Abs(y[i] - a[3]);
                if (dx > a[4]) a[4] = dx;
                if (dx < a[5]) a[5] = dx;
                a[6] += dx / n;
            }
            a[3] = Math.Sqrt(a[2] / n);
        }
    }
}
