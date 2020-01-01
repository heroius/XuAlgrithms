using System;
using System.Collections.Generic;
using System.Text;

namespace Heroius.XuAlgrithms
{
    public static class InterpolationApproximation
    {
        /// <summary>
        /// 给定n个结点Xi(i=0,1,...,n—1)上的函数值Yi=f(Xi), 用拉格朗日(Lagrange)插值公式计算指定插值点t 处的函数近似值z = f(t) 。
        /// </summary>
        /// <param name="x">存放给定n个结点的值，要求X0<X1<...<Xn-1</param>
        /// <param name="y">存放n个给定结点上的函数值Y0,Y1,...,Yn-1</param>
        /// <param name="n">给定结点的个数</param>
        /// <param name="t">指定插值点</param>
        /// <returns>指定插值点t处的函数近似值z=f(t)</returns>
        public static double LGR(double[] x, double[] y, int n, double t)
        {
            int i, j, k, m;
            double z = 0, s;
            if (n < 1) return z;
            if (n == 1)
            {
                z = y[0];
                return z;
            }
            if (n == 2)
            {
                z = (y[0] * (t - x[1]) - y[1] * (t - x[0])) / (x[0] - x[1]);
                return z;
            }
            i = 0;
            while ((x[i] < t) && (i < n)) i++;
            k = i - 4;
            if (k < 0) k = 0;
            m = i + 3;
            if (m > n - 1) m = n - 1;
            for (i = k; i <= m; i++)
            {
                s = 1;
                for (j = k; j <= m; j++)
                    if (j != i) s = s * (t - x[j]) / (x[i] - x[j]);
                z = z + s * y[i];
            }
            return z;
        }

        /// <summary>
        /// 给定n个结点Xi(i=0,1,...,n—1)上的函数值Yi=f(Xi), 用抛物插值公式计算指定插值点t 处的函数近似值z = f(t) 。
        /// </summary>
        /// <param name="x">存放给定n个结点的值，要求X0<X1<...<Xn-1</param>
        /// <param name="y">存放n个给定结点上的函数值Y0,Y1,...,Yn-1</param>
        /// <param name="n">给定结点的个数</param>
        /// <param name="t">指定插值点</param>
        /// <returns>指定插值点t处的函数近似值z=f(t)</returns>
        public static double LG3(double[] x, double[] y, int n, double t)
        {
            int i, j, k, m;
            double z = 0, s;
            if (n < 1) return z;
            if (n == 1)
            {
                z = y[0];
                return z;
            }
            if (n == 2)
            {
                z = (y[0] * (t - x[1]) - y[1] * (t - x[0])) / (x[0] - x[1]);
                return z;
            }
            if (t <= x[1])
            {
                k = 0;
                m = 2;
            }
            else if (t >= x[n - 2])
            {
                k = n - 3;
                m = n - 1;
            }
            else
            {
                k = 1;
                m = n;
                while (m - k != 1)
                {
                    i = (k + m) / 2;
                    if (t < x[i - 1]) m = i;
                    else k = i;
                }
                k--;
                m--;
                if (Math.Abs(t - x[k]) < Math.Abs(t - x[m])) k--;
                else m++;
            }
            z = 0;
            for (i = k; i <= m; i++)
            {
                s = 1;
                for (j = k; j <= m; j++)
                    if (j != i) s = s * (t - x[j]) / (x[i] - x[j]);
                z = z + s * y[i];
            }
            return z;
        }

        /// <summary>
        /// 给定n个结点Xi(i=0,1,...,n—1)上的函数值Yi=f(Xi), 用连分式插值公式计算指定插值点t 处的函数近似值z = f(t) 。
        /// </summary>
        /// <param name="x">存放给定n个结点的值，要求X0<X1<...<Xn-1</param>
        /// <param name="y">存放n个给定结点上的函数值Y0,Y1,...,Yn-1</param>
        /// <param name="n">给定结点的个数</param>
        /// <param name="t">指定插值点</param>
        /// <returns>指定插值点t处的函数近似值z=f(t)</returns>
        public static double PQS(double[] x, double[] y, int n, double t)
        {
            int i, j, k, m, l;
            double z = 0, h;
            double[] b = new double[8];
            if (n < 1) return z;
            if (n == 1)
            {
                z = y[0];
                return z;
            }
            if (n <= 8)
            {
                k = 0;
                m = n;
            }
            else if (t < x[4])
            {
                k = 0;
                m = 8;
            }
            else if (t > x[n - 5])
            {
                k = n - 8;
                m = 8;
            }
            else
            {
                k = 1;
                j = n;
                while (j - k != 1)
                {
                    i = (k + j) / 2;
                    if (t < x[i - 1]) j = i;
                    else k = i;
                }
                k -= 4;
                m = 8;
            }
            b[0] = y[k];
            for (i = 2; i <= m; i++)
            {
                h = y[i + k - 1];
                l = 0;
                j = 1;
                while ((l == 0) & (j <= i - 1))
                {
                    if (Math.Abs(h - b[j - 1]) == 0) l = 1;
                    else h = (x[i + k - 1] - x[j + k - 1]) / (h - b[j - 1]);
                    j++;
                }
                b[i - 1] = h;
                if (l != 0) b[i - 1] = 1e35;
            }
            z = b[m - 1];
            for (i = m - 1; i >= 1; i--) z = b[i - 1] + (t - x[i + k - 1]) / z;
            return z;
        }

        /// <summary>
        /// 给定n个结点Xi(i=0,1,...,n—1)上的函数值Yi=f(Xi)以及一阶导数值Y'i=f'(Xi), 用埃尔米特(Hermite)插值公式计算指定插值点t 处的函数近似值z = f(t) 。
        /// </summary>
        /// <param name="x">存放给定n个结点的值，要求X0<X1<...<Xn-1</param>
        /// <param name="y">存放n个给定结点上的函数值Y0,Y1,...,Yn-1</param>
        /// <param name="dy">存放n个给定结点上的一阶导数值Y'0,Y'1,...,Y'n-1</param>
        /// <param name="n">给定结点的个数</param>
        /// <param name="t">指定插值点</param>
        /// <returns>指定插值点t处的函数近似值z=f(t)</returns>
        public static double HMT(double[] x, double[] y, double[] dy, int n, double t)
        {
            int i, j;
            double z, p, q, s;
            z = 0;
            for (i = 1; i <= n; i++)
            {
                s = 1;
                for (j = 1; j <= n; j++)
                    if (j != i) s = s * (t - x[j - 1]) / (x[i - 1] - x[j - 1]);
                s *= s;
                p = 0;
                for (j = 1; j <= n; j++)
                    if (j != i) p = p + 1 / (x[i - 1] - x[j - 1]);
                q = y[i - 1] + (t - x[i - 1]) * (dy[i - 1] - 2 * y[i - 1] * p);
                z = z + q * s;
            }
            return z;
        }

        /// <summary>
        /// 给定n个结点Xi(i=0,1,...,n—1)上的函数值Yi=f(Xi)以及一阶导数值Y'i=f'(Xi), 用埃特金(Aitken)插值公式计算指定插值点t 处的函数近似值z = f(t) 。
        /// </summary>
        /// <param name="x">存放给定n个结点的值，要求X0<X1<...<Xn-1</param>
        /// <param name="y">存放n个给定结点上的函数值Y0,Y1,...,Yn-1</param>
        /// <param name="n">给定结点的个数</param>
        /// <param name="t">指定插值点</param>
        /// <param name="eps">插值的精度要求</param>
        /// <returns>指定插值点t处的函数近似值z=f(t)</returns>
        public static double ATK(double[] x, double[] y, int n, double t, double eps)
        {
            int i, j, k, m, l = 0;
            double z = 0;
            double[] xx = new double[10];
            double[] yy = new double[10];
            if (n < 1) return z;
            if (n == 1)
            {
                z = y[0];
                return z;
            }
            m = 10;
            if (m > n) m = n;
            if (t <= x[0]) k = 1;
            else if (t >= x[n - 1]) k = n;
            else
            {
                k = 1;
                j = n;
                while ((k - j != 1) && (k - j != -1))
                {
                    l = (k + j) / 2;
                    if (t < x[l - 1]) j = l;
                    else k = l;
                }
                if (Math.Abs(t - x[l - 1]) > Math.Abs(t - x[j - 1])) k = j;
            }
            j = 1;
            l = 0;
            for (i = 1; i <= m; i++)
            {
                k = k + j * l;
                if ((k < 1) || (k > n))
                {
                    l++;
                    j = -j;
                    k = k + j * l;
                }
                xx[i - 1] = x[k - 1];
                yy[i - 1] = y[k - 1];
                l++;
                j = -j;
            }
            i = 0;
            do
            {
                i++;
                z = yy[i];
                for (j = 0; j <= i - 1; j++)
                    z = yy[j] + (t - xx[j]) * (yy[j] - z) / (xx[j] - xx[i]);
                yy[i] = z;
            }
            while ((i != m - 1) && (Math.Abs(yy[i] - yy[i - 1]) > eps));
            return z;
        }

        /// <summary>
        /// 给定n个结点Xi(i=0,1,...,n—1)上的函数值Yi=f(Xi), 用阿克玛(Akima)方法计算指定子区间上的三次插值多项式与指定插值点t 处的函数近似值z = f(t) 。
        /// </summary>
        /// <param name="x">存放给定n个结点的值，要求X0<X1<...<Xn-1</param>
        /// <param name="y">存放n个给定结点上的函数值Y0,Y1,...,Yn-1</param>
        /// <param name="n">给定结点的个数</param>
        /// <param name="k">当k≥0，则表示只计算第k个子区间[Xk,Xk+1]上的三次多项式的系数S0,S1,S2,S3(当k≥n-1时，按k=n-2处理)；
        /// 若k＜0，则表示需要计算指定插值点t处的函数近似值f(t)，并计算插值点t所在子区间上的三次多项式的系数S0,S1,S2,S3</param>
        /// <param name="t">指定插值点</param>
        /// <param name="s">S0,S1,S2,S3返回三次多项式的系数；S4返回指定插值点t处的函数近似值(当k<0时)或任意(当k≥0时)</param>
        public static void SPL(double[] x, double[] y, int n, int k, double t, out double[] s)
        {
            int kk, m, l;
            double p, q;
            double[] u = new double[5];
            s = new double[5] { 0, 0, 0, 0, 0 };
            if (n < 1) throw new Exception("n的点数不足");
            if (n == 1)
            {
                s[0] = y[0];
                s[4] = y[0];
                return;
            }
            if (n == 2)
            {
                s[0] = y[0];
                s[1] = (y[1] - y[0]) / (x[1] - x[0]);
                if (k < 0) s[4] = (y[0] * (t - x[1]) - y[1] * (t - x[0])) / (x[0] - x[1]);
                return;
            }
            if (k < 0)
            {
                if (t <= x[1]) kk = 0;
                else if (t >= x[n - 1]) kk = n - 2;
                else
                {
                    kk = 1;
                    m = n;
                    while (((kk - m) != 1) && ((kk - m) != -1))
                    {
                        l = (kk + m) / 2;
                        if (t < x[l - 1]) m = l;
                        else kk = l;
                    }
                    kk--;
                }
            }
            else kk = k;
            if (kk >= n - 1) kk = n - 2;
            u[2] = (y[kk + 1] - y[kk]) / (x[kk + 1] - x[kk]);
            if (n == 3)
            {
                if (kk == 0)
                {
                    u[3] = (y[2] - y[1]) / (x[2] - x[1]);
                    u[4] = 2 * u[3] - u[2];
                    u[1] = 2 * u[2] - u[3];
                    u[0] = 2 * u[1] - u[2];
                }
                else
                {
                    u[1] = (y[1] - y[0]) / (x[1] - x[0]);
                    u[0] = 2 * u[1] - u[2];
                    u[3] = 2 * u[2] - u[1];
                    u[4] = 2 * u[3] - u[2];
                }
            }
            else
            {
                if (kk <= 1)
                {
                    u[3] = (y[kk + 2] - y[kk + 1]) / (x[kk + 2] - x[kk + 1]);
                    if (kk == 1)
                    {
                        u[1] = (y[1] - y[0]) / (x[1] - x[0]);
                        u[0] = 2 * u[1] - u[2];
                        if (n == 4) u[4] = 2 * u[3] - u[2];
                        else u[4] = (y[4] - y[3]) / (x[4] - x[3]);
                    }
                    else
                    {
                        u[1] = 2 * u[2] - u[3];
                        u[0] = 2 * u[1] - u[2];
                        u[4] = (y[3] - y[2]) / (x[3] - x[2]);
                    }
                }
                else if (kk >= (n - 3))
                {
                    u[1] = (y[kk] - y[kk - 1]) / (x[kk] - x[kk - 1]);
                    if (kk == (n - 3))
                    {
                        u[3] = (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);
                        u[4] = 2 * u[3] - u[2];
                        if (n == 4) u[0] = 2 * u[1] - u[2];
                        else u[0] = (y[kk - 1] - y[kk - 2]) / (x[kk - 1] - x[kk - 2]);
                    }
                    else
                    {
                        u[3] = 2 * u[2] - u[1];
                        u[4] = 2 * u[3] - u[2];
                        u[0] = (y[kk - 1] - y[kk - 2]) / (x[kk - 1] - x[kk - 2]);
                    }
                }
                else
                {
                    u[1] = (y[kk] - y[kk - 1]) / (x[kk] - x[kk - 1]);
                    u[0] = (y[kk - 1] - y[kk - 2]) / (x[kk - 1] - x[kk - 2]);
                    u[3] = (y[kk + 2] - y[kk + 1]) / (x[kk + 2] - x[kk + 1]);
                    u[4] = (y[kk + 3] - y[kk + 2]) / (x[kk + 3] - x[kk + 2]);
                }
            }
            s[0] = Math.Abs(u[3] - u[2]);
            s[1] = Math.Abs(u[0] - u[1]);
            if ((s[0] == 0) && (s[1] == 0)) p = (u[1] + u[2]) / 2;
            else p = (s[0] * u[1] + s[1] * u[2]) / (s[0] + s[1]);
            s[0] = Math.Abs(u[3] - u[4]);
            s[1] = Math.Abs(u[2] - u[1]);
            if ((s[0] == 0) && (s[1] == 0)) q = (u[2] + u[3]) / 2;
            else q = (s[0] * u[2] + s[1] * u[3]) / (s[0] + s[1]);
            s[0] = y[kk];
            s[1] = p;
            s[3] = x[kk + 1] - x[kk];
            s[2] = (3 * u[2] - 2 * p - q) / s[3];
            s[3] = (q + p - 2 * u[2]) / (s[3] * s[3]);
            if (k < 0)
            {
                p = t - x[kk];
                s[4] = s[0] + s[1] * p + s[2] * p * p + s[3] * p * p * p;
            }
        }
    }
}
