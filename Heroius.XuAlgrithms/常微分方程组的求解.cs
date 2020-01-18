using System;
using System.Collections.Generic;
using System.Text;

namespace Heroius.XuAlgrithms
{
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
                    for (k = 0; k < 3; i++)
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
                m++;
            }
        }

        /// <summary>
        /// 用变步长基尔(Gill) 公式对一阶微分方程组积分一步
        /// </summary>
        /// <param name="t">对微分方程进行积分的起始点t0</param>
        /// <param name="h">积分的步长</param>
        /// <param name="y">存放 n 个未知函数在起始点 t 处的函数值 Yj(t) (j=0,1,...,n—1)
        /// <para>返回t+h点处的n个未知函数值Yj(t+h) (j=0,1,...,n—1)</param>
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
                    nn++;
                }
                for (j = 0; j < n; j++) zz[j * k + i] = y[j];
            }
            z = Utility.C.Convert(zz, n, k);
        }
        //todo: MRSN 示例计算结果与书中差距较大


    }
}
