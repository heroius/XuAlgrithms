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
        /// <param name="func">指向计算微分方程组中各方程右端函数值的函数名</param>
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
        /// <param name="func">指向计算微分方程组中各方程右端函数值的函数名</param>
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
                    d= func(x, y,n);
                    for (i = 0; i <= n - 1; i++)
                        y[i] = c[i] + hh * d[i];
                    x = t + (j + 1) * hh;
                    d=func(x, y,n);
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
        /// <param name="func">指向计算微分方程组中各方程右端函数值的函数名</param>
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
                for (i = 0; i <= n-1; i++)
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
        /// <param name="func">指向计算微分方程组中各方程右端函数值的函数名</param>
        public static void RKT1(double t, double[] y, int n, double h, int k, out double[,] z, Func<double, double[], int, double[]> func)
        {

        }
    }
}
