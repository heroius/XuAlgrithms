using System;
using System.Collections.Generic;
using System.Text;

namespace Heroius.XuAlgrithms
{
    /// <summary>
    /// 常微分方程组求解
    /// </summary>
    public static class ODEs
    {
        /// <summary>
        /// 用改进的欧拉(Euler) 公式对一阶微分方程组进行定步长全区间积分
        /// </summary>
        /// <param name="t">对微分方程进行积分的起始点t0</param>
        /// <param name="y">存放 n 个未知函数在起始点 t 处的函数值 Yj(t) (j=0,1,...,n—1)</param>
        /// <param name="n">微分方程组中方程个数，也是未知函数的个数</param>
        /// <param name="h">积分的步长</param>
        /// <param name="k">积分步数（包括起始点这一步）</param>
        /// <param name="z">返回 k 个积分点（包括起始点）上的未知函数值</param>
        /// <param name="func">指向计算微分方程组中各方程右端函数值的函数名
        /// <para>第1参数：积分起点</para>
        /// <para>第2参数：未知函数在起始点处的函数值</para>
        /// <para>第3参数：方程数</para>
        /// <para>返回值：右端函数值</para>
        /// </param>
        public static void ELR1(double t, double[] y, int n, double h, int k, out double[,] z, Func<double, double[], int, double[]> func)
        {
            int i, j;
            double x;
            double[] d;
            double[] zz = new double[n * k];
            for (i = 0; i <= n - 1; i++) zz[i * k] = y[i];
            for (j = 1; j <= k - 1; j++)
            {
                x = t + (j - 1) * h;
                d = func(x, y, n);
                for (i = 0; i <= n - 1; i++)
                    y[i] = zz[i * k + j - 1] + h * d[i];
                x = t + j * h;
                d = func(x, y, n);
                for (i = 0; i <= n - 1; i++)
                    d[i] = zz[i * k + j - 1] + h * d[i];
                for (i = 0; i <= n - 1; i++)
                {
                    y[i] = (y[i] + d[i]) / 2;
                    zz[i * k + j] = y[i];
                }
            }
            z = Utility.C.Convert(zz, n, k);
        }

        /// <summary>
        /// 用变步长欧拉(Euler)方法对一阶微分方程组积分一步。
        /// </summary>
        /// <param name="t">对微分方程进行积分的起始点t0</param>
        /// <param name="h">积分的步长</param>
        /// <param name="y">存放 n 个未知函数在起始点 t 处的函数值 Yj(t) (j=0,1,...,n—1)
        /// <para>返回t+h点处的n个未知函数值Yj(t+h) (j=0,1,...,n—1)</para>
        /// </param>
        /// <param name="n">微分方程组中方程个数，也是未知函数的个数</param>
        /// <param name="eps">积分的精度要求</param>
        /// <param name="func">指向计算微分方程组中各方程右端函数值的函数名
        /// <para>第1参数：积分起点</para>
        /// <para>第2参数：未知函数在起始点处的函数值</para>
        /// <para>第3参数：方程数</para>
        /// <para>返回值：右端函数值</para>
        /// </param>
        public static void ELR2(double t, double h, ref double[] y, int n, double eps, Func<double, double[], int, double[]> func)
        {
            int i, j, m;
            double hh, p, x, q;
            double[]
                a = new double[n],
                b = new double[n],
                c = new double[n],
                d;
            hh = h; m = 1; p = 1.0 + eps;
            for (i = 0; i <= n - 1; i++) a[i] = y[i];
            while (p >= eps)
            {
                for (i = 0; i <= n - 1; i++)
                { b[i] = y[i]; y[i] = a[i]; }
                for (j = 0; j <= m - 1; j++)
                {
                    for (i = 0; i <= n - 1; i++) c[i] = y[i];
                    x = t + j * hh;
                    d = func(x, y, n);
                    for (i = 0; i <= n - 1; i++)
                        y[i] = c[i] + hh * d[i];
                    x = t + (j + 1) * hh;
                    d = func(x, y, n);
                    for (i = 0; i <= n - 1; i++)
                        d[i] = c[i] + hh * d[i];
                    for (i = 0; i <= n - 1; i++)
                        y[i] = (y[i] + d[i]) / 2.0;
                }
                p = 0.0;
                for (i = 0; i <= n - 1; i++)
                {
                    q = Math.Abs(y[i] - b[i]);
                    if (q > p) p = q;
                }
                hh = hh / 2.0; m = m + m;
            }
        }

        /// <summary>
        /// 用定步长维梯(Witty) 方法对一阶微分方程组进行全区间积分。
        /// </summary>
        /// <param name="t">对微分方程进行积分的起始点t0</param>
        /// <param name="y">存放 n 个未知函数在起始点 t 处的函数值 Yj(t) (j=0,1,...,n—1)</param>
        /// <param name="n">微分方程组中方程个数，也是未知函数的个数</param>
        /// <param name="h">积分的步长</param>
        /// <param name="k">积分步数（包括起始点这一步）</param>
        /// <param name="z">z[n,k] 返回k个积分点（包括起始点）上的未知函数值</param>
        /// <param name="func">指向计算微分方程组中各方程右端函数值的函数名
        /// <para>第1参数：积分起点</para>
        /// <para>第2参数：未知函数在起始点处的函数值</para>
        /// <para>第3参数：方程数</para>
        /// <para>返回值：右端函数值</para>
        /// </param>
        public static void WITY(double t, double[] y, int n, double h, int k, out double[,] z, Func<double, double[], int, double[]> func)
        {
            int i, j;
            double x;
            double[]
                a = new double[n],
                d;
            double[] zz = new double[n * k];
            for (i = 0; i <= n - 1; i++) zz[i * k] = y[i];
            d = func(t, y, n);
            for (j = 1; j <= k - 1; j++)
            {
                for (i = 0; i <= n - 1; i++)
                    a[i] = zz[i * k + j - 1] + h * d[i] / 2.0;
                x = t + (j - 0.5) * h;
                y = func(x, a, n);
                for (i = 0; i <= n - 1; i++)
                {
                    d[i] = 2.0 * y[i] - d[i];
                    zz[i * k + j] = zz[i * k + j - 1] + h * y[i];
                }
            }
            z = Utility.C.Convert(zz, n, k);
        }

        /// <summary>
        /// 用定步长四阶龙格－库塔(Runge- Kutta)法对一阶微分方程组进行全区间积分。
        /// </summary>
        /// <param name="t">对微分方程进行积分的起始点t0</param>
        /// <param name="y">存放 n 个未知函数在起始点 t 处的函数值 Yj(t) (j=0,1,...,n—1)</param>
        /// <param name="n">微分方程组中方程个数，也是未知函数的个数</param>
        /// <param name="h">积分的步长</param>
        /// <param name="k">积分步数（包括起始点这一步）</param>
        /// <param name="z">z[n,k] 返回k个积分点（包括起始点）上的未知函数值</param>
        /// <param name="func">指向计算微分方程组中各方程右端函数值的函数名
        /// <para>第1参数：积分起点</para>
        /// <para>第2参数：未知函数在起始点处的函数值</para>
        /// <para>第3参数：方程数</para>
        /// <para>返回值：右端函数值</para>
        /// </param>
        public static void RKT1(double t, double[] y, int n, double h, int k, out double[,] z, Func<double, double[], int, double[]> func)
        {
            int i, j, l;
            double[]
                a = new double[4] { h / 2, h / 2, h, h },
                b = new double[n],
                d = new double[n];
            double tt;
            double[] zz = new double[n * k];
            for (i = 0; i <= n - 1; i++) zz[i * k] = y[i];
            for (l = 1; l <= k - 1; l++)
            {
                d = func(t, y, n);
                for (i = 0; i <= n - 1; i++) b[i] = y[i];
                for (j = 0; j <= 2; j++)
                {
                    for (i = 0; i <= n - 1; i++)
                    {
                        y[i] = zz[i * k + l - 1] + a[j] * d[i];
                        b[i] = b[i] + a[j + 1] * d[i] / 3;
                    }
                    tt = t + a[j];
                    d = func(tt, y, n);
                }
                for (i = 0; i <= n - 1; i++)
                    y[i] = b[i] + h * d[i] / 6;
                for (i = 0; i <= n - 1; i++)
                    zz[i * k + l] = y[i];
                t += h;
            }
            z = Utility.C.Convert(zz, n, k);
        }

        /// <summary>
        /// 用变步长四阶龙格－库塔(Runge- Kutta) 法对一阶微分方程组积分一步
        /// </summary>
        /// <param name="t">对微分方程进行积分的起始点t0</param>
        /// <param name="h">积分的步长</param>
        /// <param name="y">存放 n 个未知函数在起始点 t 处的函数值 Yj(t) (j=0,1,...,n—1)
        /// <para>返回t+h点处的n个未知函数值Yj(t+h) (j=0,1,...,n—1)</param>
        /// <param name="n">微分方程组中方程个数，也是未知函数的个数</param>
        /// <param name="eps">积分的精度要求</param>
        /// <param name="func">指向计算微分方程组中各方程右端函数值的函数名
        /// <para>第1参数：积分起点</para>
        /// <para>第2参数：未知函数在起始点处的函数值</para>
        /// <para>第3参数：方程数</para>
        /// <para>返回值：右端函数值</para>
        /// </param>
        public static void RKT2(double t, double h, ref double[] y, int n, double eps, Func<double, double[], int, double[]> func)
        {
            int m = 1, i, j, k;
            double hh = h, p = 1 + eps, dt, x = t, tt, q;
            double[]
                a = new double[4],
                g = new double[n],
                b = new double[n],
                c = new double[n],
                d = new double[n],
                e = new double[n];
            for (i = 0; i < n; i++) c[i] = y[i];
            while (p >= eps)
            {
                a[0] = hh / 2; a[1] = a[0]; a[2] = hh; a[3] = hh;
                for (i = 0; i < n; i++)
                {
                    g[i] = y[i];
                    y[i] = c[i];
                }
                dt = h / m;
                t = x;
                for (j = 0; j < m; j++)
                {
                    d = func(t, y, n);
                    for (i = 0; i < n; i++)
                    {
                        b[i] = y[i];
                        e[i] = y[i];
                    }
                    for (k = 0; k < 3; k++)
                    {
                        for (i = 0; i < n; i++)
                        {
                            y[i] = e[i] + a[k] * d[i];
                            b[i] = b[i] + a[k + 1] * d[i] / 3;
                        }
                        tt = t + a[k];
                        d = func(tt, y, n);
                    }
                    for (i = 0; i < n; i++)
                    {
                        y[i] = b[i] + hh * d[i] / 6;
                    }
                    t += dt;
                }
                p = 0;
                for (i = 0; i < n; i++)
                {
                    q = Math.Abs(y[i] - g[i]);
                    if (q > p) p = q;
                }
                hh = hh / 2;
                m += m;
            }
        }

        /// <summary>
        /// 用变步长基尔(Gill) 公式对一阶微分方程组积分一步
        /// </summary>
        /// <param name="t">对微分方程进行积分的起始点t0</param>
        /// <param name="h">积分的步长</param>
        /// <param name="y">存放 n 个未知函数在起始点 t 处的函数值 Yj(t) (j=0,1,...,n—1)
        /// <para>返回t+h点处的n个未知函数值Yj(t+h) (j=0,1,...,n—1)</para></param>
        /// <param name="n">微分方程组中方程个数，也是未知函数的个数</param>
        /// <param name="eps">积分的精度要求</param>
        /// <paramref name="q">在主函数第一次调用本函数时，应赋值以0, 
        /// 即q[i]=0(i=0,1,…,n-1),以后每调用一次本函数（即每积分一步），将由本函数的返回值以便循环使用</paramref>
        /// <param name="func">指向计算微分方程组中各方程右端函数值的函数名
        /// <para>第1参数：积分起点</para>
        /// <para>第2参数：未知函数在起始点处的函数值</para>
        /// <para>第3参数：方程数</para>
        /// <para>返回值：右端函数值</para>
        /// </param>
        public static void GIL(double t, double h, ref double[] y, int n, double eps, ref double[] q, Func<double, double[], int, double[]> func)
        {
            int j, k, m = 1, ii;
            double x = t, p = 1 + eps, hh = h, r, s, t0, dt, qq;
            double[]
                a = new double[4] { 0.5, 0.29289321881, 1.7071067812, 0.166666667 },
                b = new double[4] { 2.0, 1.0, 1.0, 2.0 },
                c = new double[4] { a[0], a[1], a[2], 0.5 },
                e = new double[4] { 0.5, 0.5, 1.0, 1.0 },
                d = new double[n],
                u = new double[n],
                v = new double[n],
                g = new double[n];
            for (j = 0; j <= n - 1; j++) u[j] = y[j];
            while (p >= eps)
            {
                for (j = 0; j <= n - 1; j++)
                { v[j] = y[j]; y[j] = u[j]; g[j] = q[j]; }
                dt = h / m; t = x;
                for (k = 0; k <= m - 1; k++)
                {
                    d = func(t, y, n);
                    for (ii = 0; ii <= 3; ii++)
                    {
                        for (j = 0; j <= n - 1; j++) d[j] = d[j] * hh;
                        for (j = 0; j <= n - 1; j++)
                        {
                            r = (a[ii] * (d[j] - b[ii] * g[j]) + y[j]) - y[j];
                            y[j] = y[j] + r;
                            s = g[j] + 3.0 * r;
                            g[j] = s - c[ii] * d[j];
                        }
                        t0 = t + e[ii] * hh;
                        d = func(t0, y, n);
                    }
                    t = t + dt;
                }
                p = 0.0;
                for (j = 0; j <= n - 1; j++)
                {
                    qq = Math.Abs(y[j] - v[j]);
                    if (qq > p) p = qq;
                }
                hh = hh / 2.0; m = m + m;
            }
            for (j = 0; j <= n - 1; j++)
            {
                q[j] = g[j];
            }
        }

        /// <summary>
        /// 用变步长默森（Merson）方法对一阶微分方程组进行全区间积分
        /// </summary>
        /// <param name="t">对微分方程进行积分的起始点t0</param>
        /// <param name="h">积分的步长</param>
        /// <param name="n">微分方程组中方程个数，也是未知函数的个数</param>
        /// <param name="y">存放 n 个未知函数在起始点 t 处的函数值 Yj(t) (j=0,1,...,n—1)</param>
        /// <param name="eps">积分的精度要求</param>
        /// <param name="k">积分步数（包括起始点这一步）</param>
        /// <param name="z">z[n,k]返回k个积分点（包括起始点）的未知函数值</param>
        /// <param name="func">指向计算微分方程组中各方程右端函数值的函数名
        /// <para>第1参数：积分起点</para>
        /// <para>第2参数：未知函数在起始点处的函数值</para>
        /// <para>第3参数：方程数</para>
        /// <para>返回值：右端函数值</para></param>
        public static void MRSN(double t, double h, int n, double[] y, double eps, int k, out double[,] z, Func<double, double[], int, double[]> func)
        {
            int i, j, m, nn;
            double aa = t, bb, x, hh, p, dt, t0, qq;
            double[] d,
                a = new double[n],
                b = new double[n],
                c = new double[n],
                u = new double[n],
                v = new double[n],
                zz = new double[n * k];
            for (i = 0; i < n; i++) zz[i * k] = y[i];
            for (i = 1; i < k; i++)
            {
                x = aa + (i - 1) * h;
                nn = 1;
                hh = h;
                for (j = 0; j < n; j++) u[j] = y[j];
                p = 1.0 + eps;
                while (p >= eps)
                {
                    for (j = 0; j < n; j++)
                    {
                        v[j] = y[j];
                        y[j] = u[j];
                    }
                    dt = h / nn;
                    t = x;
                    for (m = 0; m < nn; m++)
                    {
                        d = func(t, y, n);
                        for (j = 0; j < n; j++)
                        {
                            a[j] = d[j];
                            y[j] = y[j] + hh * d[j] / 3.0;
                        }
                        t0 = t + hh / 3.0;
                        d = func(t0, y, n);
                        for (j = 0; j < n; j++)
                        {
                            b[j] = d[j];
                            y[j] = y[j] + hh * (d[j] - a[j]) / 6.0;
                        }
                        d = func(t0, y, n);
                        for (j = 0; j < n; j++)
                        {
                            b[j] = d[j];
                            bb = (d[j] - 4.0 * (b[j] + a[j] / 4.0) / 9.0) / 8.0;
                            y[j] = y[j] + 3.0 * hh * bb;
                        }
                        t0 = t + hh / 2.0;
                        d = func(t0, y, n);
                        for (j = 0; j < n; j++)
                        {
                            c[j] = d[j];
                            qq = d[j] - 15.0 * (b[j] - a[j] / 5.0) / 16.0;
                            y[j] = y[j] + 2.0 * hh * qq;
                        }
                        t0 = t + hh;
                        d = func(t0, y, n);
                        for (j = 0; j < n; j++)
                        {
                            qq = c[j] - 9.0 * (b[j] - 2.0 * a[j] / 9.0) / 8.0;
                            qq = d[j] - 8.0 * qq;
                            y[j] = y[j] + hh * qq / 6.0;
                        }
                        t = t + dt;
                    }
                    p = 0.0;
                    for (j = 0; j < n; j++)
                    {
                        qq = Math.Abs(y[j] - v[j]);
                        if (qq > p) p = qq;
                    }
                    hh = hh / 2.0;
                    nn = nn + nn;
                }
                for (j = 0; j < n; j++) zz[j * k + i] = y[j];
            }
            z = Utility.C.Convert(zz, n, k);
        }

        /// <summary>
        /// 用连分式法对一阶微分方程组积分一步。
        /// </summary>
        /// <param name="t">对微分方程进行积分的起始点t0</param>
        /// <param name="h">积分的步长</param>
        /// <param name="n">微分方程组中方程个数，也是未知函数的个数</param>
        /// <param name="y">存放 n 个未知函数在起始点 t 处的函数值 Yj(t) (j=0,1,...,n—1)
        /// <para>返回t+h点处的n个未知函数值Yj(t+h) (j=0,1,...,n—1)</para></param>
        /// <param name="eps">积分的精度要求</param>
        /// <param name="func">指向计算微分方程组中各方程右端函数值的函数名
        /// <para>第1参数：积分起点</para>
        /// <para>第2参数：未知函数在起始点处的函数值</para>
        /// <para>第3参数：方程数</para>
        /// <para>返回值：右端函数值</para></param>
        public static void PBS(double t, double h, int n, ref double[] y, double eps, Func<double, double[], int, double[]> func)
        {
            int k = 1, m, nn = 1, it = 1;
            double x = t, hh = h, dd, q, p;
            double[]
                g = new double[10],
                b = new double[10 * n],
                d = new double[n],
                u = new double[n],
                v = new double[n],
                w = new double[n],
                e = new double[n];

            for (int j = 0; j < n; j++)
            {
                v[j] = y[j];
            }
            g[0] = hh;
            PBS_RKT(x, hh, n, ref y, w, ref d, e, func);
            for (int i = 0; i < n; i++)
            {
                b[i] = y[i]; u[i] = y[i];
            }
            while (it == 1)
            {
                nn += nn;
                hh /= 2;
                it = 0;
                g[k] = hh;
                for (int i = 0; i < n; i++)
                {
                    y[i] = v[i];
                }
                t = x;
                for (int i = 0; i < nn; i++)
                {
                    PBS_RKT(t, hh, n, ref y, w, ref d, e, func);
                    t += hh;
                }
                for (int i = 0; i < n; i++)
                {
                    dd = y[i];
                    m = 0;
                    for (int j = 0; j < k; j++)
                    {
                        if (m == 0)
                        {
                            q = dd - b[j * n + i];
                            if (q == 0)
                            {
                                m = 1;
                            }
                            else
                            {
                                dd = (g[k] - g[j]) / q;
                            }
                        }
                    }
                    b[k * n + i] = dd;
                    if (m != 0) b[k * n + i] = 1e35;
                }
                for (int j = 0; j < n; j++)
                {
                    dd = 0;
                    for (int i = k - 1; i >= 0; i--)
                    {
                        dd = -g[i] / b[(i + 1) * n + j] + dd;
                    }
                    y[j] = dd + b[j];
                }
                p = 0;
                for (int j = 0; j < n; j++)
                {
                    q = Math.Abs(y[j] - u[j]);
                    if (q > p) p = q;
                }
                if ((p >= eps) && (k < 7))
                {
                    for (int j = 0; j < n; j++)
                    {
                        u[j] = y[j];
                    }
                    k++;
                    it = 1;
                }
            }
        }

        static void PBS_RKT(double t, double h, int n, ref double[] y, double[] b, ref double[] d, double[] e, Func<double, double[], int, double[]> func)
        {
            int i, k;
            double[] a = new double[] { h / 2, h / 2, h, h };
            double tt;
            for (i = 0; i < n; i++)
            {
                b[i] = y[i];
                e[i] = y[i];
            }
            for (k = 0; k < 3; k++)
            {
                for (i = 0; i < n; i++)
                {
                    y[i] = e[i] + a[k] * d[i];
                    b[i] = b[i] + a[k + 1] * d[i] / 3;
                }
                tt = t + a[k];
                d = func(tt, y, n);
            }
            for (i = 0; i < n; i++)
                y[i] = b[i] + h * d[i] / 6;
        }

        /// <summary>
        /// 用双边法对一阶微分方程组进行全区间积分。
        /// </summary>
        /// <param name="t">对微分方程进行积分的起始点t0</param>
        /// <param name="h">积分的步长</param>
        /// <param name="n">微分方程组中方程个数，也是未知函数的个数</param>
        /// <param name="y">存放 n 个未知函数在起始点 t 处的函数值 Yj(t) (j=0,1,...,n—1)</param>
        /// <param name="eps">积分的精度要求</param>
        /// <param name="k">积分步数（包括起始点这一步）</param>
        /// <param name="z">z[n,k]返回k个积分点（包括起始点）的未知函数值</param>
        /// <param name="func">指向计算微分方程组中各方程右端函数值的函数名
        /// <para>第1参数：积分起点</para>
        /// <para>第2参数：未知函数在起始点处的函数值</para>
        /// <para>第3参数：方程数</para>
        /// <para>返回值：右端函数值</para></param>
        public static void GJFQ(double t, double h, int n, double[] y, double eps, int k, out double[,] z, Func<double, double[], int, double[]> func)
        {
            double a = t, qq;
            double[]
                d = new double[n],
                p = new double[n],
                u = new double[n],
                v = new double[n],
                w = new double[n],
                zz = new double[n * k];
            for (int i = 0; i < n; i++)
            {
                p[i] = 0;
                zz[i * k] = y[i];
            }
            d = func(t, y, n);
            for (int j = 0; j < n; j++)
            {
                u[j] = d[j];
            }
            RKT2(t, h, ref y, n, eps, func);
            t = a + h;
            d = func(t, y, n);
            for (int j = 0; j < n; j++)
            {
                zz[j * k + 1] = y[j];
                v[j] = d[j];
            }
            for (int j = 0; j < n; j++)
            {
                p[j] = -4 * zz[j * k + 1] + 5 * zz[j * k] + 2 * h * (2 * v[j] + u[j]);
                y[j] = p[j];
            }
            t = a + 2 * h;
            d = func(t, y, n);
            for (int j = 0; j < n; j++)
            {
                qq = 2 * h * (d[j] - 2 * v[j] - 2 * u[j]) / 3;
                qq += 4 * zz[j * k + 1] - 3 * zz[j * k];
                zz[j * k + 2] = (p[j] + qq) / 2;
                y[j] = zz[j * k + 2];
            }
            for (int i = 3; i < k; i++)
            {
                t = a + (i - 1) * h;
                d = func(t, y, n);
                for (int j = 0; j < n; j++)
                {
                    u[j] = v[j];
                    v[j] = d[j];
                }
                for (int j = 0; j < n; j++)
                {
                    qq = -4 * zz[j * k + i - 1] + 5 * zz[j * k + i - 2];
                    p[j] = qq + 2 * h * (2 * v[j] + u[j]);
                    y[j] = p[j];
                }
                t += h;
                d = func(t, y, n);
                for (int j = 0; j < n; j++)
                {
                    qq = 2 * h * (d[j] - 2 * v[j] - 2 * u[j]) / 3;
                    qq += 4 * zz[j * k + i - 1] - 3 * zz[j * k + i - 2];
                    y[j] = (p[j] + qq) / 2;
                    zz[j * k + i] = y[j];
                }
            }
            z = Utility.C.Convert(zz, n, k);
        }

        /// <summary>
        /// 用阿当姆斯(Adams) 预报校正公式对一阶微分方程组进行全区间积分
        /// </summary>
        /// <param name="t">对微分方程进行积分的起始点t0</param>
        /// <param name="h">积分的步长</param>
        /// <param name="n">微分方程组中方程个数，也是未知函数的个数</param>
        /// <param name="y">存放 n 个未知函数在起始点 t 处的函数值 Yj(t) (j=0,1,...,n—1)</param>
        /// <param name="eps">积分的精度要求</param>
        /// <param name="k">积分步数（包括起始点这一步）</param>
        /// <param name="z">z[n,k]返回k个积分点（包括起始点）的未知函数值</param>
        /// <param name="func">指向计算微分方程组中各方程右端函数值的函数名
        /// <para>第1参数：积分起点</para>
        /// <para>第2参数：未知函数在起始点处的函数值</para>
        /// <para>第3参数：方程数</para>
        /// <para>返回值：右端函数值</para></param>
        public static void ADMS(double t, double h, int n, double[] y, double eps, int k, out double[,] z, Func<double, double[], int, double[]> func)
        {
            double a = t, q;
            double[]
                b = new double[4 * n],
                e = new double[n],
                s = new double[n],
                g = new double[n],
                d = new double[n],
                zz = new double[n * k];
            for (int i = 0; i < n; i++)
            {
                zz[i * k] = y[i];
            }
            d = func(t, y, n);
            for (int i = 0; i < n; i++)
            {
                b[i] = d[i];
            }
            for (int i = 0; i < 4; i++)
            {
                if (i <= k - 1)
                {
                    t = a + i * h;
                    RKT2(t, h, ref y, n, eps, func);
                    for (int j = 0; j < n; j++)
                    {
                        zz[j * k + i] = y[j];
                    }
                    d = func(t, y, n);
                    for (int j = 0; j < n; j++)
                    {
                        b[i * n + j] = d[j];
                    }
                }
            }
            for (int i = 4; i < k; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    q = 55 * b[3 * n + j] - 59 * b[2 * n + j];
                    q += 37 * b[n + j] - 9 * b[j];
                    y[j] = zz[j * k + i - 1] + h * q / 24;
                    b[j] = b[n + j];
                    b[n + j] = b[n + n + j];
                    b[n + n + j] = b[n + n + n + j];
                }
                t = a + i * h;
                d = func(t, y, n);
                for (int m = 0; m < n; m++)
                {
                    b[n + n + n + m] = d[m];
                }
                for (int j = 0; j < n; j++)
                {
                    q = 9 * b[3 * n + j] + 19 * b[n + n + j] - 5 * b[n + j] + b[j];
                    y[j] = zz[j * k + i - 1] + h * q / 24;
                    zz[j * k + i] = y[i];
                }
                d = func(t, y, n);
                for (int m = 0; m < n; m++)
                {
                    b[3 * n + m] = d[m];
                }
            }
            z = Utility.C.Convert(zz, n, k);
        }

        /// <summary>
        /// 用哈明(Hamming) 方法对一阶微分方程组进行全区间积分
        /// </summary>
        /// <param name="t">对微分方程进行积分的起始点t0</param>
        /// <param name="h">积分的步长</param>
        /// <param name="n">微分方程组中方程个数，也是未知函数的个数</param>
        /// <param name="y">存放 n 个未知函数在起始点 t 处的函数值 Yj(t) (j=0,1,...,n—1)</param>
        /// <param name="eps">积分的精度要求</param>
        /// <param name="k">积分步数（包括起始点这一步）</param>
        /// <param name="z">z[n,k]返回k个积分点（包括起始点）的未知函数值</param>
        /// <param name="func">指向计算微分方程组中各方程右端函数值的函数名
        /// <para>第1参数：积分起点</para>
        /// <para>第2参数：未知函数在起始点处的函数值</para>
        /// <para>第3参数：方程数</para>
        /// <para>返回值：右端函数值</para></param>
        public static void HAMG(double t, double h, int n, double[] y, double eps, int k, out double[,] z, Func<double, double[], int, double[]> func)
        {
            double a = t, q;
            double[]
                b = new double[4 * n],
                d = new double[n],
                u = new double[n],
                v = new double[n],
                w = new double[n],
                g = new double[n],
                zz = new double[n * k];
            for (int i = 0; i < n; i++)
            {
                zz[i * k] = y[i];
            }
            d = func(t, y, n);
            for (int i = 0; i < n; i++)
            {
                b[i] = d[i];
            }
            for (int i = 1; i < 4; i++)
            {
                if (i <= k - 1)
                {
                    t = a + i * h;
                    RKT2(t, h, ref y, n, eps, func);
                    for (int m = 0; m < n; m++)
                    {
                        zz[m * k + i] = y[m];
                    }
                    d = func(t, y, n);
                    for (int m = 0; m < n; m++)
                    {
                        b[i * n + m] = d[m];
                    }
                }
            }
            for (int i = 0; i < n; i++)
            {
                u[i] = 0;
            }
            for (int i = 4; i < k; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    q = 2 * b[3 * n + j] - b[n + n + j] + 2 * b[n + j];
                    y[j] = zz[j * k + i - 4] + 4 * h * q / 3;
                }
                for (int j = 0; j < n; j++)
                {
                    y[j] += 112 * u[j] / 121;
                }
                t = a + i * h;
                d = func(t, y, n);
                for (int j = 0; j < n; j++)
                {
                    q = 9 * zz[j * k + i - 1] - zz[j * k + i - 3];
                    q = (q + 3 * h * (d[j] + 2 * b[3 * n + j] - b[n + n + j])) / 8;
                    u[j] = q - y[j];
                    zz[j * k + i] = q - 9 * u[j] / 121;
                    y[j] = zz[j * k + i];
                    b[n + j] = b[n + n + j];
                    b[n + n + j] = b[n + n + n + j];
                }
                d = func(t, y, n);
                for (int m = 0; m < n; m++)
                {
                    b[3 * n + m] = d[m];
                }
            }
            z = Utility.C.Convert(zz, n, k);
        }

        /// <summary>
        /// 用特雷纳(Treanor) 方法对一阶刚性(stiff) 微分方程组积分一步
        /// </summary>
        /// <param name="t">对微分方程进行积分的起始点t0</param>
        /// <param name="h">积分的步长</param>
        /// <param name="n">微分方程组中方程个数，也是未知函数的个数</param>
        /// <param name="y">存放 n 个未知函数在起始点 t 处的函数值 Yj(t) (j=0,1,...,n—1)
        /// <para>返回t+h点处的n个未知函数值Yj(t+h) (j=0,1,...,n—1)</para></param>
        /// <param name="func">指向计算微分方程组中各方程右端函数值的函数名
        /// <para>第1参数：积分起点</para>
        /// <para>第2参数：未知函数在起始点处的函数值</para>
        /// <para>第3参数：方程数</para>
        /// <para>返回值：右端函数值</para></param>
        public static void TNR(double t, double h, int n, ref double[] y, Func<double, double[], int, double[]> func)
        {
            double s = t + h / 2, aa, bb, dd, g, dy, dy1;
            double[]
                w = new double[4 * n],
                q = new double[4 * n],
                r = new double[4 * n],
                d = new double[n],
                p = new double[n];
            for (int j = 0; j < n; j++)
            {
                w[j] = y[j];
            }
            d = func(t, y, n);
            for (int j = 0; j < n; j++)
            {
                q[j] = d[j];
                y[j] = w[j] + h * d[j] / 2;
                w[n + j] = y[j];
            }
            d = func(s, y, n);
            for (int j = 0; j < n; j++)
            {
                q[n + j] = d[j];
                y[j] = w[j] + h * d[j] / 2;
                w[n + n + j] = y[j];
            }
            d = func(s, y, n);
            for (int j = 0; j < n; j++)
            {
                q[n + n + j] = d[j];
            }
            for (int j = 0; j < n; j++)
            {
                aa = q[n + n + j] - q[n + j];
                bb = w[n + n + j] - w[n + j];
                if (-aa * bb * h > 0)
                {
                    p[j] = -aa / bb;
                    dd = -p[j] * h;
                    r[j] = Math.Exp(dd);
                    r[n + j] = (r[j] - 1) / dd;
                    r[n + n + j] = (r[n + j] - 1) / dd;
                    r[3 * n + j] = (r[n + n + j] - 1) / dd;
                }
                else
                {
                    p[j] = 0;
                }
                if (p[j] <= 0) g = q[n + n + j];
                else
                {
                    g = 2 * (q[n + n + j] - q[j]) * r[n + n + j];
                    g += (q[j] - q[n + j]) * r[n + j] + q[n + j];
                }
                w[3 * n + j] = w[j] + g * h;
                y[j] = w[3 * n + j];
            }
            s = t + h;
            d = func(s, y, n);
            for (int j = 0; j < n; j++)
            {
                q[3 * n + j] = d[j];
            }
            for (int j = 0; j < n; j++)
            {
                if (p[j] <= 0)
                {
                    dy = q[j] + 2 * (q[n + j] + q[n + n + j]);
                    dy = (dy + q[n + n + n + j]) * h / 6;
                }
                else
                {
                    dy = -3 * (q[j] + p[j] * w[j]) + 2 * (q[n + j] + p[j] * w[n + j]);
                    dy += 2 * (q[n + n + j] + p[j] * w[n + n + j]);
                    dy -= (q[n + n + n + j] + p[j] * w[n + n + n + j]);
                    dy = dy * r[n + n + j] + q[j] * r[n + j];
                    dy1 = q[j] - q[n + j] - q[n + n + j] + q[n + n + n + j];
                    dy1 += (w[j] - w[n + j] - w[n + n + j] + w[n + n + n + j]) * p[j];
                    dy = (dy + 4 * dy1 * r[n + n + n + j]) * h;
                }
                y[j] = w[j] + dy;
            }
        }

        /// <summary>
        /// 用吉尔(Gear) 方法积分一阶微分方程组的初值问题。
        /// </summary>
        /// <param name="a">积分区间起始点</param>
        /// <param name="b">积分区间终点</param>
        /// <param name="hmin">积分过程中所允许的最小步长</param>
        /// <param name="hmax">积分过程中所允许的最大步长</param>
        /// <param name="h">积分的拟定步长，在积分过程中将自动放大或缩小
        /// <para>hmin<<h<hmax</para></param>
        /// <param name="eps">误差检验常数</param>
        /// <param name="n">微分方程组中方程个数，也是未知函数的个数</param>
        /// <param name="y0">存放 n 个未知函数在起始点 t 处的函数值 Yj(t) (j=0,1,...,n—1)</param>
        /// <param name="k">拟定输出的积分点数</param>
        /// <param name="t">返回k 个输出点处的自变量值（包括起始点）</param>
        /// <param name="z">返回k 个输出点处的未知函数值</param>
        /// <param name="ss">指向计算雅可比矩阵的函数名</param>
        /// <param name="f">指向计算微分方程组中各方程右端函数值的函数名
        /// <para>第1参数：积分起点</para>
        /// <para>第2参数：未知函数在起始点处的函数值</para>
        /// <para>第3参数：方程数</para>
        /// <para>返回值：右端函数值</para></param>
        /// <returns>函数返回标志值。其意义如下：
        /// <para>返回值=+l 表示全区间积分成功</para>
        /// <para>返回值=-l 表示步长己小于最小步长hmin, 但精度还达不到，积分停止</para>
        /// <para>返回值=-2 表示阶数已大于6,积分停止</para>
        /// <para>返回值=-3 表示对于h≥hmin 校正迭代步收敛，积分停止</para>
        /// <para>返回值=-4 表示对所处理的问题，要求的精度太高，积分停止</para></returns>
        public static int GEAR(double a, double b, double hmin, double hmax, double h, double eps, int n, double[] y0, int k, out double[] t, out double[,] z, Func<double, double[], int, double[,]> ss, Func<double, double[], int, double[]> f)
        {
            int kf, jt, nn, nq, i, m, irt1, irt, j, nqd, idb=0, jjs = 0, flag;
            int iw=0, j1, j2, nt, nqw, l;
            double hw, hd, rm, t0, td, r = 0, dd, pr1, pr2, pr3, rr, h0;
            double enq1=0, enq2=0, enq3=0, eup=0, e=0, edwn=0, bnd=0, r1;
            double[]
                aa = new double[7],
                d = new double[n],
                s = new double[10 * n],
                s02 = new double[n],
                ym = new double[n],
                er = new double[n],
                yy = new double[n],
                y = new double[8 * n];
            t = new double[k];
            int[]
                _is = new int[n],
                js = new int[n];
            double[,]
                pp = new double[7, 3]{ {2.0,3.0,1.0},{4.5,6.0,1.0},
                       {7.333,9.167,0.5},{10.42,12.5,0.1667},
                       {13.7,15.98,0.04133},{17.15,1.0,0.008267},
                       {1.0,1.0,1.0}},
                p = new double[n, n];
            z = new double[n, k];
            aa[1] = -1.0; jt = 0; nn = 0; nq = 1; t0 = a; h0 = h;
            for (i = 0; i <= 8 * n - 1; i++) y[i] = 0.0;
            for (i = 0; i <= n - 1; i++)
            { y[i * 8] = y0[i]; yy[i] = y[i * 8]; }
            d = f(t0, yy, n);
            for (i = 0; i <= n - 1; i++) y[i * 8 + 1] = h * d[i];
            hw = h; m = 2;
            for (i = 0; i <= n - 1; i++) ym[i] = 1.0;
            l20:
            irt = 1; kf = 1; nn = nn + 1;
            t[nn - 1] = t0;
            for (i = 0; i <= n - 1; i++) z[i, nn - 1] = y[i * 8];
            if ((t0 >= b) || (nn == k))
            { flag = kf; h = h0; return flag; }
            for (i = 0; i <= n - 1; i++)
                for (j = 0; j <= m - 1; j++) s[i * 10 + j] = y[i * 8 + j];
            hd = hw;
            if (h != hd)
            {
                rm = h / hd; irt1 = 0;
                rr = Math.Abs(hmin / hd);
                if (rm < rr) rm = rr;
                rr = Math.Abs(hmax / hd);
                if (rm > rr) rm = rr;
                r = 1.0; irt1 = irt1 + 1;
                for (j = 1; j <= m - 1; j++)
                {
                    r = r * rm;
                    for (i = 0; i <= n - 1; i++) y[i * 8 + j] = s[i * 10 + j] * r;
                }
                h = hd * rm;
                for (i = 0; i <= n - 1; i++) y[i * 8] = s[i * 10];
                idb = m;
            }
            nqd = nq; td = t0; rm = 1.0;
            if (jt > 0) goto l80;
            l60:
            switch (nq)
            {
                case 1: aa[0] = -1.0; break;
                case 2: aa[0] = -2.0 / 3.0; aa[2] = -1.0 / 3.0; break;
                case 3:
                    aa[0] = -6.0 / 11.0; aa[2] = aa[0];
                    aa[3] = -1.0 / 11.0; break;
                case 4:
                    aa[0] = -0.48; aa[2] = -0.7; aa[3] = -0.2;
                    aa[4] = -0.02; break;
                case 5:
                    aa[0] = -120.0 / 274.0; aa[2] = -225.0 / 274.0;
                    aa[3] = -85.0 / 274.0; aa[4] = -15.0 / 274.0;
                    aa[5] = -1.0 / 274.0; break;
                case 6:
                    aa[0] = -720.0 / 1764.0; aa[2] = -1624.0 / 1764.0;
                    aa[3] = -735.0 / 1764.0; aa[4] = -175.0 / 1764.0;
                    aa[5] = -21.0 / 1764.0; aa[6] = -1.0 / 1764.0;
                    break;
                default: { flag = -2; h = h0; return flag; }
            }
            m = nq + 1; idb = m;
            enq2 = 0.5 / (nq + 1.0); enq3 = 0.5 / (nq + 2.0);
            enq1 = 0.5 / (nq + 0.0);
            eup = pp[nq - 1, 1] * eps; eup = eup * eup;
            e = pp[nq - 1, 0] * eps; e = e * e;
            edwn = pp[nq - 1, 2] * eps; edwn = edwn * edwn;
            if (edwn == 0.0)
            {
                for (i = 0; i <= n - 1; i++)
                    for (j = 0; j <= m - 1; j++) y[i * 8 + j] = s[i * 10 + j];
                h = hd; nq = nqd; jt = nq;
                flag = -4; h = h0; return flag;
            }
            bnd = eps * enq3 / (n + 0.0);
            iw = 1;
            if (irt == 2)
            {
                r1 = 1.0;
                for (j = 1; j <= m - 1; j++)
                {
                    r1 = r1 * r;
                    for (i = 0; i <= n - 1; i++) y[i * 8 + j] = y[i * 8 + j] * r1;
                }
                idb = m;
                for (i = 0; i <= n - 1; i++)
                    if (ym[i] < Math.Abs(y[i * 8])) ym[i] = Math.Abs(y[i * 8]);
                jt = nq;
                goto l20;
            }
        l80:
            t0 = t0 + h;
            for (j = 2; j <= m; j++)
                for (j1 = j; j1 <= m; j1++)
                {
                    j2 = m - j1 + j - 1;
                    for (i = 0; i <= n - 1; i++)
                        y[i * 8 + j2 - 1] = y[i * 8 + j2 - 1] + y[i * 8 + j2];
                }
            for (i = 0; i <= n - 1; i++) er[i] = 0.0;
            j1 = 1; nt = 1;
            for (l = 0; l <= 2; l++)
            {
                if ((j1 != 0) && (nt != 0))
                {
                    for (i = 0; i <= n - 1; i++) yy[i] = y[i * 8];
                    d = f(t0, yy, n);
                    if (iw >= 1)
                    {
                        for (i = 0; i <= n - 1; i++) yy[i] = y[i * 8];
                        p = ss(t0, yy, n);
                        r = aa[0] * h;
                        for (i = 0; i <= n - 1; i++)
                            for (j = 0; j <= n - 1; j++) p[i,j] = p[i,j] * r;
                        for (i = 0; i <= n - 1; i++) p[i,i] = 1.0 + p[i,i];
                        iw = -1;
                        jjs = Matrix.RINV(ref p, n);
                        j1 = jjs;
                    }
                    if (jjs != 0)
                    {
                        for (i = 0; i <= n - 1; i++) s02[i] = y[i * 8 + 1] - d[i] * h;
                        for (i = 0; i <= n - 1; i++)
                        {
                            dd = 0.0;
                            for (j = 0; j <= n - 1; j++) dd = dd + s02[j] * p[i,j];
                            s[i * 10 + 8] = dd;
                        }
                        nt = n;
                        for (i = 0; i <= n - 1; i++)
                        {
                            y[i * 8] = y[i * 8] + aa[0] * s[i * 10 + 8];
                            y[i * 8 + 1] = y[i * 8 + 1] - s[i * 10 + 8];
                            er[i] = er[i] + s[i * 10 + 8];
                            if (Math.Abs(s[i * 10 + 8]) <= (bnd * ym[i])) nt = nt - 1;
                        }
                    }
                }
            }
            if (nt > 0)
            {
                t0 = td;
                if ((h > (hmin * 1.00001)) || (iw >= 0))
                {
                    if (iw != 0) rm = 0.25 * rm;
                    iw = 1; irt1 = 2;
                    rr = Math.Abs(hmin / hd);
                    if (rm < rr) rm = rr;
                    rr = Math.Abs(hmax / hd);
                    if (rm > rr) rm = rr;
                    r = 1.0;
                    for (j = 1; j <= m - 1; j++)
                    {
                        r = r * rm;
                        for (i = 0; i <= n - 1; i++) y[i * 8 + j] = s[i * 10 + j] * r;
                    }
                    h = hd * rm;
                    for (i = 0; i <= n - 1; i++) y[i * 8] = s[i * 10];
                    idb = m;
                    goto l80;
                }
                for (i = 0; i <= n - 1; i++)
                    for (j = 0; j <= m - 1; j++) y[i * 8 + j] = s[i * 10 + j];
                h = hd; nq = nqd; jt = nq;
                flag = -3;
                h = h0; return flag;
            }
            dd = 0.0;
            for (i = 0; i <= n - 1; i++) dd = dd + (er[i] / ym[i]) * (er[i] / ym[i]);
            iw = 0;
            if (dd <= e)
            {
                if (m >= 3)
                    for (j = 2; j <= m - 1; j++)
                        for (i = 0; i <= n - 1; i++)
                            y[i * 8 + j] = y[i * 8 + j] + aa[j] * er[i];
                kf = 1; hw = h;
                if (idb > 1)
                {
                    idb = idb - 1;
                    if (idb <= 1)
                        for (i = 0; i <= n - 1; i++) s[i * 10 + 9] = er[i];
                    for (i = 0; i <= n - 1; i++)
                        if (ym[i] < Math.Abs(y[i * 8])) ym[i] = Math.Abs(y[i * 8]);
                    jt = nq;
                    goto l20;
                }
            }
            if (dd > e)
            {
                kf = kf - 2;
                if (h <= (hmin * 1.00001))
                {
                    hw = h; jt = nq;
                    flag = -1;
                    h = h0; return flag;
                }
                t0 = td;
                if (kf <= -5)
                {
                    if (nq == 1)
                    {
                        for (i = 0; i <= n - 1; i++)
                            for (j = 0; j <= m - 1; j++) y[i * 8 + j] = s[i * 10 + j];
                        h = hd; nq = nqd; jt = nq;
                        flag = -4;
                        h = h0; return flag;
                    }
                    for (i = 0; i <= n - 1; i++) yy[i] = y[i * 8];
                    d = f(t0, yy, n);
                    r = h / hd;
                    for (i = 0; i <= n - 1; i++)
                    {
                        y[i * 8] = s[i * 10];
                        s[i * 10 + 1] = hd * d[i];
                        y[i * 8 + 1] = s[i * 10 + 1] * r;
                    }
                    nq = 1; kf = 1; goto l60;
                }
            }
            pr2 = Math.Log(dd / e); pr2 = enq2 * pr2; pr2 = Math.Exp(pr2);
            pr2 = 1.2 * pr2;
            pr3 = 1.0e+20;
            if (nq < 7)
                if (kf > -1)
                {
                    dd = 0.0;
                    for (i = 0; i <= n - 1; i++)
                    {
                        pr3 = (er[i] - s[i * 10 + 9]) / ym[i];
                        dd = dd + pr3 * pr3;
                    }
                    pr3 = Math.Log(dd / eup); pr3 = enq3 * pr3;
                    pr3 = Math.Exp(pr3); pr3 = 1.4 * pr3;
                }
            pr1 = 1.0e+20;
            if (nq > 1)
            {
                dd = 0.0;
                for (i = 0; i <= n - 1; i++)
                {
                    pr1 = y[i * 8 + m - 1] / ym[i];
                    dd = dd + pr1 * pr1;
                }
                pr1 = Math.Log(dd / edwn); pr1 = enq1 * pr1;
                pr1 = Math.Exp(pr1); pr1 = 1.3 * pr1;
            }
            if (pr2 <= pr3)
            {
                if (pr2 > pr1)
                {
                    r = 1.0e+04;
                    if (pr1 > 1.0e-04) r = 1.0 / pr1;
                    nqw = nq - 1;
                }
                else
                {
                    nqw = nq; r = 1.0e+04;
                    if (pr2 > 1.0e-04) r = 1.0 / pr2;
                }
            }
            else
            {
                if (pr3 < pr1)
                {
                    r = 1.0e+04;
                    if (pr3 > 1.0e-04) r = 1.0 / pr3;
                    nqw = nq + 1;
                }
                else
                {
                    r = 1.0e+04;
                    if (pr1 > 1.0e-04) r = 1.0 / pr1;
                    nqw = nq - 1;
                }
            }
            idb = 10;
            if (kf == 1)
                if (r < 1.1)
                {
                    for (i = 0; i <= n - 1; i++)
                        if (ym[i] < Math.Abs(y[i * 8])) ym[i] = Math.Abs(y[i * 8]);
                    jt = nq; goto l20;
                }
            if (nqw > nq)
                for (i = 0; i <= n - 1; i++) y[i * 8 + nqw] = er[i] * aa[m - 1] / (m + 0.0);
            m = nqw + 1;
            if (kf == 1)
            {
                irt = 2; rr = hmax / Math.Abs(h);
                if (r > rr) r = rr;
                h = h * r; hw = h;
                if (nq == nqw)
                {
                    r1 = 1.0;
                    for (j = 1; j <= m - 1; j++)
                    {
                        r1 = r1 * r;
                        for (i = 0; i <= n - 1; i++) y[i * 8 + j] = y[i * 8 + j] * r1;
                    }
                    idb = m;
                    for (i = 0; i <= n - 1; i++)
                        if (ym[i] < Math.Abs(y[i * 8])) ym[i] = Math.Abs(y[i * 8]);
                    jt = nq; goto l20;
                }
                nq = nqw;
                goto l60;
            }
            rm = rm * r; irt1 = 3;
            rr = Math.Abs(hmin / hd);
            if (rm < rr) rm = rr;
            rr = Math.Abs(hmax / hd);
            if (rm > rr) rm = rr;
            r = 1.0;
            for (j = 1; j <= m - 1; j++)
            {
                r = r * rm;
                for (i = 0; i <= n - 1; i++) y[i * 8 + j] = s[i * 10 + j] * r;
            }
            h = hd * rm;
            for (i = 0; i <= n - 1; i++) y[i * 8] = s[i * 10];
            idb = m;
            if (nqw == nq) goto l80;
            nq = nqw; goto l60;
        }

        /// <summary>
        /// 用有限差分法求二阶线性微分方程边值问题的数值解。
        /// </summary>
        /// <param name="a">求解区间的左端点</param>
        /// <param name="b">求解区间的右端点</param>
        /// <param name="ya">未知函数在求解区间左端点处的函数值y(a)</param>
        /// <param name="yb">未知函数在求解区间右端点处的函数值y(b)</param>
        /// <param name="n">求解区间[a, b] 的等分点数（包括左右端点a 与b)</param>
        /// <param name="y">返回n 个等距离散点上的未知函数值</param>
        /// <param name="f">指向计算二阶微分方程中的函数u(x)、v(x)、w(x)、f(x) 值的函数
        /// <para>输入参数：自变量x值</para>
        /// <para>返回值：各函数值</para>
        /// </param>
        public static void DFTE(double a, double b, double ya, double yb, int n, out double[] y, Func<double, double[]> f) {
            int j, k, nn=2*n-1, m1;
            double h = (b-a)/(n-1), x;
            double[]
                z = new double[4],
                g = new double[6 * n],
                d = new double[2 * n];
            g[0] = 1;
            g[1] = 0;
            y = new double[n];
            y[0] = ya;
            y[n - 1] = yb;
            g[3 * n - 3] = 1;
            g[3 * n - 4] = 0;
            for (int i = 2; i < n; i++)
            {
                x = a + (i - 1) * h;
                z = f(x);
                k = 3 * (i - 1) - 1;
                g[k] = z[0] - h * z[1] / 2;
                g[k + 1] = h * h * z[2] - 2 * z[0];
                g[k + 2] = z[0] + h * z[1] / 2;
                y[i - 1] = h * h * z[3];
            }
            m1 = 3 * n - 2;
            LinearEquations.TRDE(g, n, m1, ref y);
            h = h / 2;
            g[0] = 1;
            g[1] = 0;
            d[0] = ya;
            d[nn - 1] = yb;
            g[3 * nn - 3] = 1;
            g[3 * nn - 4] = 0;
            for (int i = 2; i < nn; i++)
            {
                x = a + (i - 1) * h;
                z = f(x);
                k = 3 * (i - 1) - 1;
                g[k] = z[0] - h * z[1] / 2;
                g[k + 1] = h * h * z[2] - 2 * z[0];
                g[k + 2] = z[0] + h * z[1] / 2;
                d[i - 1] = h * h * z[3];
            }
            m1 = 3 * nn - 2;
            LinearEquations.TRDE(g, nn, m1, ref d);
            for (int i = 2; i < n; i++)
            {
                k = i + i - 1;
                y[i - 1] = (4 * d[k - 1] - y[i - 1]) / 3;
            }
        }
    }
}
