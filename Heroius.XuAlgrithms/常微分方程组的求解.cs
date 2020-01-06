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


    }
}
