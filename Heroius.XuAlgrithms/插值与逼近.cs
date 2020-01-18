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

        /// <summary>
        /// 给定n个结点Xi(i=0,1,...,n—1)上的函数值Yi=f(Xi)以及两端点上的一阶导数值Y'0=f'(X0)与Y'n-1=f'(Xn-1)，
        /// 利用三次样条函数计算各结点上的数值导数以及插值区间[X0,Xn-1]上的积分近似值s=∫f(x)dx，
        /// 并对函数f(x)进行成组插值与成组微商。
        /// </summary>
        /// <param name="x">存放给定n个结点的值</param>
        /// <param name="y">存放n个给定结点上的函数值</param>
        /// <param name="n">给定结点的个数</param>
        /// <param name="dy">dy[0]存放Y'0, dy[n-1]存放Y'n-1。
        /// 返回n个给定节点处的一阶导数值Y'k(k=0,...,n-1)</param>
        /// <param name="ddy">返回n个给定结点处的二阶导数值Y''k(k=0,...,n-1)</param>
        /// <param name="t">存放m个指定插值点的值，要求X0<t(j)<Xn-1(j=0,...,m-1)</param>
        /// <param name="m">指定插值点的个数</param>
        /// <param name="z">返回m个指定插值点处的函数值</param>
        /// <param name="dz">返回m个指定插值点处的一阶导数值</param>
        /// <param name="ddz">返回m个指定插值点处的二阶导数值</param>
        /// <returns>返回积分值s=∫f(x)dx</returns>
        public static double SPL1(double[] x, double[] y, int n, ref double[] dy, out double[] ddy, double[] t, int m, out double[] z, out double[] dz, out double[] ddz)
        {
            ddy = new double[n];
            z = new double[m];
            dz = new double[m];
            ddz = new double[m];

            int i, j;
            double h0, h1, alpha, beta, integ;

            double[] s = new double[n];

            s[0] = dy[0]; dy[0] = 0.0;
            h0 = x[1] - x[0];
            for (j = 1; j <= n - 2; j++)
            {
                h1 = x[j + 1] - x[j];
                alpha = h0 / (h0 + h1);
                beta = (1.0 - alpha) * (y[j] - y[j - 1]) / h0;
                beta = 3.0 * (beta + alpha * (y[j + 1] - y[j]) / h1);
                dy[j] = -alpha / (2.0 + (1.0 - alpha) * dy[j - 1]);
                s[j] = (beta - (1.0 - alpha) * s[j - 1]);
                s[j] = s[j] / (2.0 + (1.0 - alpha) * dy[j - 1]);
                h0 = h1;
            }
            for (j = n - 2; j >= 0; j--)
                dy[j] = dy[j] * dy[j + 1] + s[j];
            for (j = 0; j <= n - 2; j++) s[j] = x[j + 1] - x[j];
            for (j = 0; j <= n - 2; j++)
            {
                h1 = s[j] * s[j];
                ddy[j] = 6.0 * (y[j + 1] - y[j]) / h1 - 2.0 * (2.0 * dy[j] + dy[j + 1]) / s[j];
            }
            h1 = s[n - 2] * s[n - 2];
            ddy[n - 1] = 6.0 * (y[n - 2] - y[n - 1]) / h1 + 2.0 * (2.0 * dy[n - 1] + dy[n - 2]) / s[n - 2]; //origin code .* err
            integ = 0.0;
            for (i = 0; i <= n - 2; i++)
            {
                h1 = 0.5 * s[i] * (y[i] + y[i + 1]);
                h1 = h1 - s[i] * s[i] * s[i] * (ddy[i] + ddy[i + 1]) / 24.0;
                integ = integ + h1;
            }
            for (j = 0; j <= m - 1; j++)
            {
                if (t[j] >= x[n - 1]) i = n - 2;
                else
                {
                    i = 0;
                    while (t[j] > x[i + 1]) i = i + 1;
                }
                h1 = (x[i + 1] - t[j]) / s[i];
                h0 = h1 * h1;
                z[j] = (3.0 * h0 - 2.0 * h0 * h1) * y[i];
                z[j] = z[j] + s[i] * (h0 - h0 * h1) * dy[i];
                dz[j] = 6.0 * (h0 - h1) * y[i] / s[i];
                dz[j] = dz[j] + (3.0 * h0 - 2.0 * h1) * dy[i];
                ddz[j] = (6.0 - 12.0 * h1) * y[i] / (s[i] * s[i]);
                ddz[j] = ddz[j] + (2.0 - 6.0 * h1) * dy[i] / s[i];
                h1 = (t[j] - x[i]) / s[i];
                h0 = h1 * h1;
                z[j] = z[j] + (3.0 * h0 - 2.0 * h0 * h1) * y[i + 1];
                z[j] = z[j] - s[i] * (h0 - h0 * h1) * dy[i + 1];
                dz[j] = dz[j] - 6.0 * (h0 - h1) * y[i + 1] / s[i];
                dz[j] = dz[j] + (3.0 * h0 - 2.0 * h1) * dy[i + 1];
                ddz[j] = ddz[j] + (6.0 - 12.0 * h1) * y[i + 1] / (s[i] * s[i]);
                ddz[j] = ddz[j] - (2.0 - 6.0 * h1) * dy[i + 1] / s[i];
            }
            return integ;
        }

        /// <summary>
        /// 给定n个结点Xi(i=0,1,...,n—1)上的函数值Yi=f(Xi)以及两端点上的二阶导数值Y''0=f''(X0)与Y''n-1=f''(Xn-1)，
        /// 利用三次样条函数计算各结点上的数值导数以及插值区间[X0,Xn-1]上的积分近似值s=∫f(x)dx，
        /// 并对函数f(x)进行成组插值与成组微商。
        /// </summary>
        /// <param name="x">存放给定n个结点的值</param>
        /// <param name="y">存放n个给定结点上的函数值</param>
        /// <param name="n">给定结点的个数</param>
        /// <param name="dy">返回n个给定结点处的一阶导数值Y'k(k=0,...,n-1)</param>
        /// <param name="ddy">ddy[0]存放Y''0, ddy[n-1]存放Y''n-1。
        /// 返回n个给定节点处的二阶导数值Y''k(k=0,...,n-1)</param>
        /// <param name="t">存放m个指定插值点的值，要求X0<t(j)<Xn-1(j=0,...,m-1)</param>
        /// <param name="m">指定插值点的个数</param>
        /// <param name="z">返回m个指定插值点处的函数值</param>
        /// <param name="dz">返回m个指定插值点处的一阶导数值</param>
        /// <param name="ddz">返回m个指定插值点处的二阶导数值</param>
        /// <returns>返回积分值s=∫f(x)dx</returns>
        public static double SPL2(double[] x, double[] y, int n, out double[] dy, ref double[] ddy, double[] t, int m, out double[] z, out double[] dz, out double[] ddz)
        {
            dy = new double[n];
            z = new double[m];
            dz = new double[m];
            ddz = new double[m];

            int i, j;
            double h0, h1 = 0, alpha, beta, integ;

            double[] s = new double[n];

            dy[0] = -0.5;
            h0 = x[1] - x[0];
            s[0] = 3.0 * (y[1] - y[0]) / (2.0 * h0) - ddy[0] * h0 / 4.0;
            for (j = 1; j <= n - 2; j++)
            {
                h1 = x[j + 1] - x[j];
                alpha = h0 / (h0 + h1);
                beta = (1.0 - alpha) * (y[j] - y[j - 1]) / h0;
                beta = 3.0 * (beta + alpha * (y[j + 1] - y[j]) / h1);
                dy[j] = -alpha / (2.0 + (1.0 - alpha) * dy[j - 1]);
                s[j] = (beta - (1.0 - alpha) * s[j - 1]);
                s[j] = s[j] / (2.0 + (1.0 - alpha) * dy[j - 1]);
                h0 = h1;
            }
            dy[n - 1] = (3.0 * (y[n - 1] - y[n - 2]) / h1 + ddy[n - 1] * h1 /
                  2.0 - s[n - 2]) / (2.0 + dy[n - 2]);
            for (j = n - 2; j >= 0; j--)
                dy[j] = dy[j] * dy[j + 1] + s[j];
            for (j = 0; j <= n - 2; j++) s[j] = x[j + 1] - x[j];
            for (j = 0; j <= n - 2; j++)
            {
                h1 = s[j] * s[j];
                ddy[j] = 6.0 * (y[j + 1] - y[j]) / h1 - 2.0 * (2.0 * dy[j] + dy[j + 1]) / s[j];
            }
            h1 = s[n - 2] * s[n - 2];
            ddy[n - 1] = 6.0 * (y[n - 2] - y[n - 1]) / h1 + 2.0 * (2.0 * dy[n - 1] + dy[n - 2]) / s[n - 2];
            integ = 0.0;
            for (i = 0; i <= n - 2; i++)
            {
                h1 = 0.5 * s[i] * (y[i] + y[i + 1]);
                h1 = h1 - s[i] * s[i] * s[i] * (ddy[i] + ddy[i + 1]) / 24.0;
                integ = integ + h1;
            }
            for (j = 0; j <= m - 1; j++)
            {
                if (t[j] >= x[n - 1]) i = n - 2;
                else
                {
                    i = 0;
                    while (t[j] > x[i + 1]) i = i + 1;
                }
                h1 = (x[i + 1] - t[j]) / s[i];
                h0 = h1 * h1;
                z[j] = (3.0 * h0 - 2.0 * h0 * h1) * y[i];
                z[j] = z[j] + s[i] * (h0 - h0 * h1) * dy[i];
                dz[j] = 6.0 * (h0 - h1) * y[i] / s[i];
                dz[j] = dz[j] + (3.0 * h0 - 2.0 * h1) * dy[i];
                ddz[j] = (6.0 - 12.0 * h1) * y[i] / (s[i] * s[i]);
                ddz[j] = ddz[j] + (2.0 - 6.0 * h1) * dy[i] / s[i];
                h1 = (t[j] - x[i]) / s[i];
                h0 = h1 * h1;
                z[j] = z[j] + (3.0 * h0 - 2.0 * h0 * h1) * y[i + 1];
                z[j] = z[j] - s[i] * (h0 - h0 * h1) * dy[i + 1];
                dz[j] = dz[j] - 6.0 * (h0 - h1) * y[i + 1] / s[i];
                dz[j] = dz[j] + (3.0 * h0 - 2.0 * h1) * dy[i + 1];
                ddz[j] = ddz[j] + (6.0 - 12.0 * h1) * y[i + 1] / (s[i] * s[i]);
                ddz[j] = ddz[j] - (2.0 - 6.0 * h1) * dy[i + 1] / s[i];
            }
            return integ;
        }

        /// <summary>
        /// 给定n个结点Xi(i=0,1,...,n—1)上的函数值Yi=f(Xi)以及
        /// 第三种边界条件，
        /// 利用三次样条函数计算各结点上的数值导数以及插值区间[X0,Xn-1]上的积分近似值s=∫f(x)dx，
        /// 并对函数f(x)进行成组插值与成组微商。
        /// </summary>
        /// <param name="x">存放给定n个结点的值</param>
        /// <param name="y">存放n个给定结点上的函数值，要求Y0=Yn-1</param>
        /// <param name="n">给定结点的个数</param>
        /// <param name="dy">返回n个给定结点处的一阶导数值Y'k(k=0,...,n-1)</param>
        /// <param name="ddy">返回n个给定结点处的二阶导数值Y''k(k=0,...,n-1)</param>
        /// <param name="t">存放m个指定插值点的值，要求X0<t(j)<Xn-1(j=0,...,m-1)</param>
        /// <param name="m">指定插值点的个数</param>
        /// <param name="z">返回m个指定插值点处的函数值</param>
        /// <param name="dz">返回m个指定插值点处的一阶导数值</param>
        /// <param name="ddz">返回m个指定插值点处的二阶导数值</param>
        /// <returns>返回积分值s=∫f(x)dx</returns>
        public static double SPL3(double[] x, double[] y, int n, out double[] dy, out double[] ddy, double[] t, int m, out double[] z, out double[] dz, out double[] ddz)
        {
            dy = new double[n];
            ddy = new double[n];
            z = new double[m];
            dz = new double[m];
            ddz = new double[m];

            int i, j;
            double
                h0 = x[n - 1] - x[n - 2],
                y0 = y[n - 1] - y[n - 2],
                h1,
                y1, alpha = 0, beta = 0, u, g;

            double[] s = new double[n];

            dy[0] = 0;
            ddy[0] = 0.0;
            ddy[n - 1] = 0.0;
            s[0] = 1.0;
            s[n - 1] = 1.0;
            for (j = 1; j <= n - 1; j++)
            {
                h1 = h0;
                y1 = y0;
                h0 = x[j] - x[j - 1];
                y0 = y[j] - y[j - 1];
                alpha = h1 / (h1 + h0);
                beta = 3.0 * ((1.0 - alpha) * y1 / h1 + alpha * y0 / h0);
                if (j < n - 1)
                {
                    u = 2.0 + (1.0 - alpha) * dy[j - 1];
                    dy[j] = -alpha / u;
                    s[j] = (alpha - 1.0) * s[j - 1] / u;
                    ddy[j] = (beta - (1.0 - alpha) * ddy[j - 1]) / u;
                }
            }
            for (j = n - 2; j >= 1; j--)
            {
                s[j] = dy[j] * s[j + 1] + s[j];
                ddy[j] = dy[j] * ddy[j + 1] + ddy[j];
            }
            dy[n - 2] = (beta - alpha * ddy[1] - (1.0 - alpha) * ddy[n - 2]) /
                (alpha * s[1] + (1.0 - alpha) * s[n - 2] + 2.0);
            for (j = 2; j <= n - 1; j++)
                dy[j - 2] = s[j - 1] * dy[n - 2] + ddy[j - 1];
            dy[n - 1] = dy[0];
            for (j = 0; j <= n - 2; j++) s[j] = x[j + 1] - x[j];
            for (j = 0; j <= n - 2; j++)
            {
                h1 = s[j] * s[j];
                ddy[j] = 6.0 * (y[j + 1] - y[j]) / h1 - 2.0 * (2.0 * dy[j] + dy[j + 1]) / s[j];
            }
            h1 = s[n - 2] * s[n - 2];
            ddy[n - 1] = 6.0 * (y[n - 2] - y[n - 1]) / h1 + 2.0 * (2.0 * dy[n - 1] + dy[n - 2]) / s[n - 2];
            g = 0.0;
            for (i = 0; i <= n - 2; i++)
            {
                h1 = 0.5 * s[i] * (y[i] + y[i + 1]);
                h1 = h1 - s[i] * s[i] * s[i] * (ddy[i] + ddy[i + 1]) / 24.0;
                g = g + h1;
            } //todo: SPL3: g suppose to be 5.07754e-15
            for (j = 0; j <= m - 1; j++)
            {
                h0 = t[j];
                while (h0 >= x[n - 1]) h0 = h0 - (x[n - 1] - x[0]);
                while (h0 < x[0]) h0 = h0 + (x[n - 1] - x[0]);
                i = 0;
                while (h0 > x[i + 1]) i++;
                u = h0;
                h1 = (x[i + 1] - u) / s[i];
                h0 = h1 * h1;
                z[j] = (3.0 * h0 - 2.0 * h0 * h1) * y[i];
                z[j] = z[j] + s[i] * (h0 - h0 * h1) * dy[i];
                dz[j] = 6.0 * (h0 - h1) * y[i] / s[i];
                dz[j] = dz[j] + (3.0 * h0 - 2.0 * h1) * dy[i];
                ddz[j] = (6.0 - 12.0 * h1) * y[i] / (s[i] * s[i]);
                ddz[j] = ddz[j] + (2.0 - 6.0 * h1) * dy[i] / s[i];
                h1 = (u - x[i]) / s[i];
                h0 = h1 * h1;
                z[j] = z[j] + (3.0 * h0 - 2.0 * h0 * h1) * y[i + 1];
                z[j] = z[j] - s[i] * (h0 - h0 * h1) * dy[i + 1];
                dz[j] = dz[j] - 6.0 * (h0 - h1) * y[i + 1] / s[i];
                dz[j] = dz[j] + (3.0 * h0 - 2.0 * h1) * dy[i + 1];
                ddz[j] = ddz[j] + (6.0 - 12.0 * h1) * y[i + 1] / (s[i] * s[i]);
                ddz[j] = ddz[j] - (2.0 - 6.0 * h1) * dy[i + 1] / s[i];
            }
            return g;
        }

        /// <summary>
        /// 给定矩形域上 n×m 个结点 (Xk,Yj)(k=0,l,…,n-1; j=0,1,…,m-1) 上的
        /// 函数值 Zkj = z(Xk, Yj), 利用二元三点插值公式计算指定插值点 (u, v) 处的函数值 w =z(u, v) 。
        /// </summary>
        /// <param name="x">存放给定 n×m 个结点 X 方向的 n 个坐标</param>
        /// <param name="y">存放给定 n×m 个结点 Y 方向的 m 个坐标</param>
        /// <param name="z">存放给定 n×m 个结点上的函数值</param>
        /// <param name="n">给定结点在 X 方向上的坐标个数</param>
        /// <param name="m">给定结点在 Y 方向上的坐标个数</param>
        /// <param name="u">指定插值点的 X 坐标</param>
        /// <param name="v">指定插值点的 Y 坐标</param>
        /// <returns>函数返回指定插值点 (u, v) 处的函数近似值</returns>
        public static double SLQ3(double[] x, double[] y, double[,] z, int n, int m, double u, double v)
        {
            int nn = 3, mm = 3, ip, iq, i, j, k, l;
            double h, w;
            double[] b = new double[3];
            double[] zz = Utility.C.Convert(z);
            if (n <= 3) { ip = 0; nn = n; }
            else if (u <= x[1]) ip = 0;
            else if (u >= x[n - 2]) ip = n - 3;
            else
            {
                i = 1;
                j = n;
                while (((i - j) != 1) && ((i - j) != -1))
                {
                    l = (i + j) / 2;
                    if (u < x[l - 1]) j = l;
                    else i = l;
                }
                if (Math.Abs(u - x[i - 1]) < Math.Abs(u - x[j - 1])) ip = i - 2;
                else ip = i - 1;
            }
            if (m <= 3) { iq = 0; mm = m; }
            else if (v <= y[1]) iq = 0;
            else if (v >= y[m - 2]) iq = m - 3;
            else
            {
                i = 1;
                j = m;
                while (((i - j) != 1) && ((i - j) != -1))
                {
                    l = (i + j) / 2;
                    if (v < y[l - 1]) j = l;
                    else i = l;
                }
                if (Math.Abs(v - y[i - 1]) < Math.Abs(v - y[j - 1])) iq = i - 2;
                else iq = i - 1;
            }
            for (i = 0; i <= nn - 1; i++)
            {
                b[i] = 0;
                for (j = 0; j <= mm - 1; j++)
                {
                    k = m * (ip + i) + (iq + j);
                    h = zz[k];
                    for (k = 0; k <= mm - 1; k++)
                        if (k != j)
                            h = h * (v - y[iq + k]) / (y[iq + j] - y[iq + k]);
                    b[i] = b[i] + h;
                }
            }
            w = 0;
            for (i = 0; i <= nn - 1; i++)
            {
                h = b[i];
                for (j = 0; j <= nn - 1; j++)
                    if (j != i)
                        h = h * (u - x[ip + j]) / (x[ip + i] - x[ip + j]);
                w += h;
            }
            return w;
        }

        /// <summary>
        /// 给定矩形域上 n×m 个结点 (Xk,Yj)(k=0,l,…,n-1; j=0,1,…,m-1) 上的
        /// 函数值 Zkj = z(Xk, Yj), 利用二元插值公式计算指定插值点 (u, v) 处的函数值 w =z(u, v) 。
        /// </summary>
        /// <param name="x">存放给定 n×m 个结点 X 方向的 n 个坐标</param>
        /// <param name="y">存放给定 n×m 个结点 Y 方向的 m 个坐标</param>
        /// <param name="z">存放给定 n×m 个结点上的函数值</param>
        /// <param name="n">给定结点在 X 方向上的坐标个数</param>
        /// <param name="m">给定结点在 Y 方向上的坐标个数</param>
        /// <param name="u">指定插值点的 X 坐标</param>
        /// <param name="v">指定插值点的 Y 坐标</param>
        /// <returns>函数返回指定插值点 (u, v) 处的函数近似值</returns>
        public static double SLGQ(double[] x, double[] y, double[,] z, int n, int m, double u, double v)
        {
            int ip, ipp, i, j, kk, iq, iqq, k;
            double h; double[] b = new double[10];
            double f;

            if (u <= x[0]) { ip = 1; ipp = 4; }
            else if (u >= x[n - 1]) { ip = n - 3; ipp = n; }
            else
            {
                i = 1; j = n;
                while (((i - j) != 1) && ((i - j) != -1))
                {
                    kk = (i + j) / 2;
                    if (u < x[kk - 1]) j = kk;
                    else i = kk;
                }
                ip = i - 3; ipp = i + 4;
            }
            if (ip < 1) ip = 1;
            if (ipp > n) ipp = n;
            if (v <= y[0]) { iq = 1; iqq = 4; }
            else if (v >= y[m - 1]) { iq = m - 3; iqq = m; }
            else
            {
                i = 1; j = m;
                while (((i - j) != 1) && ((i - j) != -1))
                {
                    kk = (i + j) / 2;
                    if (v < y[kk - 1]) j = kk;
                    else i = kk;
                }
                iq = i - 3; iqq = i + 4;
            }
            if (iq < 1) iq = 1;
            if (iqq > m) iqq = m;
            for (i = ip - 1; i <= ipp - 1; i++)
            {
                b[i - ip + 1] = 0.0;
                for (j = iq - 1; j <= iqq - 1; j++)
                {
                    h = z[i, j];
                    for (k = iq - 1; k <= iqq - 1; k++)
                        if (k != j) h = h * (v - y[k]) / (y[j] - y[k]);
                    b[i - ip + 1] = b[i - ip + 1] + h;
                }
            }
            f = 0.0;
            for (i = ip - 1; i <= ipp - 1; i++)
            {
                h = b[i - ip + 1];
                for (j = ip - 1; j <= ipp - 1; j++)
                    if (j != i) h = h * (u - x[j]) / (x[i] - x[j]);
                f = f + h;
            }
            return f;
        }

        /// <summary>
        /// 用最小二乘法求给定数据点的拟合多项式
        /// </summary>
        /// <param name="x">存放给定 n 个数据点的 X 坐标</param>
        /// <param name="y">存放给定 n 个数据点的 Y 坐标</param>
        /// <param name="n">给定数据点的个数</param>
        /// <param name="a">返回 m-1 次拟合多项式的 m 个系数</param>
        /// <param name="m">拟合多项式的项数，即拟合多项式的最高次为 m-1 。要求 m≤n 且 m≤20 。若
        /// m>n 或 m>20, 则本函数自动按 m = min{n,20} 处理</param>
        /// <param name="dt">dt[0] 返回拟合多项式与各数据点误差的平方和，dt[1]返回拟合多项式与各数据
        /// 点误差的绝对值之和，dt[2]返回拟合多项式与各数据点误差绝对值的最大值</param>
        public static void PIR1(double[] x, double[] y, int n, out double[] a, int m, out double[] dt)
        {
            int i, j, k;
            double z, p, c, g, q = 0, d1, d2;
            double[] s = new double[20], t = new double[20], b = new double[20];

            a = new double[m];
            dt = new double[3];

            for (i = 0; i <= m; i++) a[i] = 0.0;
            if (m + 1 > n) m = n - 1;
            if (m > 19) m = 19;
            z = 0.0;
            for (i = 0; i <= n - 1; i++) z = z + x[i] / (1.0 * n);
            b[0] = 1.0; d1 = 1.0 * n; p = 0.0; c = 0.0;
            for (i = 0; i <= n - 1; i++)
            { p = p + (x[i] - z); c = c + y[i]; }
            c = c / d1; p = p / d1;
            a[0] = c * b[0];
            if (m > 0)
            {
                t[1] = 1.0; t[0] = -p;
                d2 = 0.0; c = 0.0; g = 0.0;
                for (i = 0; i <= n - 1; i++)
                {
                    q = x[i] - z - p; d2 = d2 + q * q;
                    c = c + y[i] * q;
                    g = g + (x[i] - z) * q * q;
                }
                c = c / d2; p = g / d2; q = d2 / d1;
                d1 = d2;
                a[1] = c * t[1]; a[0] = c * t[0] + a[0];
            }
            for (j = 2; j <= m; j++)
            {
                s[j] = t[j - 1];
                s[j - 1] = -p * t[j - 1] + t[j - 2];
                if (j >= 3)
                    for (k = j - 2; k >= 1; k--)
                        s[k] = -p * t[k] + t[k - 1] - q * b[k];
                s[0] = -p * t[0] - q * b[0];
                d2 = 0.0; c = 0.0; g = 0.0;
                for (i = 0; i <= n - 1; i++)
                {
                    q = s[j];
                    for (k = j - 1; k >= 0; k--)
                        q = q * (x[i] - z) + s[k];
                    d2 = d2 + q * q; c = c + y[i] * q;
                    g = g + (x[i] - z) * q * q;
                }
                c = c / d2; p = g / d2; q = d2 / d1;
                d1 = d2;
                a[j] = c * s[j]; t[j] = s[j];
                for (k = j - 1; k >= 0; k--)
                {
                    a[k] = c * s[k] + a[k];
                    b[k] = t[k]; t[k] = s[k];
                }
            }
            dt[0] = 0.0; dt[1] = 0.0; dt[2] = 0.0;
            for (i = 0; i <= n - 1; i++)
            {
                q = a[m];
                for (k = m - 1; k >= 0; k--)
                    q = a[k] + q * (x[i] - z);
                p = q - y[i];
                if (Math.Abs(p) > dt[2]) dt[2] = Math.Abs(p);
                dt[0] = dt[0] + p * p;
                dt[1] = dt[1] + Math.Abs(p);
            }
        }

        /// <summary>
        /// 给定 n 个数据点，求切比雪夫 (Chebyshev) 意义上的最佳拟合多项式 。
        /// </summary>
        /// <param name="x">存放给定 n 个数据点的 X 坐标</param>
        /// <param name="y">存放给定 n 个数据点的 Y 坐标</param>
        /// <param name="n">给定数据点的个数</param>
        /// <param name="a">前 m 个 元素返回 m-1 次拟合多项式的 m 个系数；
        /// 最后一个元素 a[m] 返回拟合多项式 Pm-1(x) 的偏差最大值。
        /// 若 a[m] 为负值，则说明在迭代过程中参考偏差不再增大，其绝对值为当前选择的参考偏差</param>
        /// <param name="m">拟合多项式的项数，即拟合多项式的最高次为 m-1 。要求 m≤n 且 m≤20 。若
        /// m>n 或 m>20, 则本函数自动按 m = min{n,20} 处理</param>
        public static void CHIR(double[] x, double[] y, int n, out double[] a, int m)
        {
            int m1, i, j, l, ii, k, im;
            int[] ix = new int[21];
            double ha, hh, y1, y2, h1, h2, d, hm;
            double[] h = new double[21];

            a = new double[m + 1];

            for (i = 0; i <= m + 1; i++) a[i] = 0.0;
            if (m >= n) m = n - 2;
            if (m >= 20) m = 18;
            m1 = m + 2;
            ha = 0.0;
            ix[0] = 0; ix[m + 1] = n - 1;
            l = (n - 1) / (m + 1); j = l;
            for (i = 1; i <= m; i++)
            { ix[i] = j; j = j + l; }
            while (1 == 1)
            {
                hh = 1.0;
                for (i = 0; i <= m + 1; i++)
                { a[i] = y[ix[i]]; h[i] = -hh; hh = -hh; }
                for (j = 1; j <= m + 1; j++)
                {
                    ii = m1; y2 = a[ii - 1]; h2 = h[ii - 1];
                    for (i = j; i <= m + 1; i++)
                    {
                        d = x[ix[ii - 1]] - x[ix[m1 - i - 1]];
                        y1 = a[m - i + j];
                        h1 = h[m - i + j];
                        a[ii - 1] = (y2 - y1) / d;
                        h[ii - 1] = (h2 - h1) / d;
                        ii = m - i + j + 1; y2 = y1; h2 = h1;
                    }
                }
                hh = -a[m + 1] / h[m + 1];
                for (i = 0; i <= m + 1; i++)
                    a[i] = a[i] + h[i] * hh;
                for (j = 1; j <= m; j++)
                {
                    ii = m - j + 1; d = x[ix[ii - 1]];
                    y2 = a[ii - 1];
                    for (k = m1 - j; k <= m + 1; k++)
                    {
                        y1 = a[k - 1]; a[ii - 1] = y2 - d * y1;
                        y2 = y1; ii = k;
                    }
                }
                hm = Math.Abs(hh);
                if (hm <= ha) { a[m + 1] = -hm; return; }
                a[m + 1] = hm; ha = hm; im = ix[0]; h1 = hh;
                j = 0;
                for (i = 0; i <= n - 1; i++)
                {
                    if (i == ix[j])
                    { if (j < m + 1) j = j + 1; }
                    else
                    {
                        h2 = a[m];
                        for (k = m - 1; k >= 0; k--)
                            h2 = h2 * x[i] + a[k];
                        h2 = h2 - y[i];
                        if (Math.Abs(h2) > hm)
                        { hm = Math.Abs(h2); h1 = h2; im = i; }
                    }
                }
                if (im == ix[0]) return;
                i = 0; l = 1;
                while (l == 1)
                {
                    l = 0;
                    if (im >= ix[i])
                    {
                        i = i + 1;
                        if (i <= m + 1) l = 1;
                    }
                }
                if (i > m + 1) i = m + 1;
                if (i == (i / 2) * 2) h2 = -hh;
                else h2 = hh;
                if (h1 * h2 >= 0.0) ix[i] = im;
                else
                {
                    if (im < ix[0])
                    {
                        for (j = m; j >= 0; j--)
                            ix[j + 1] = ix[j];
                        ix[0] = im;
                    }
                    else
                    {
                        if (im > ix[m + 1])
                        {
                            for (j = 1; j <= m + 1; j++)
                                ix[j - 1] = ix[j];
                            ix[m + 1] = im;
                        }
                        else ix[i - 1] = im;
                    }
                }
            }
        }

        /// <summary>
        /// 用里米兹 (Remez) 方法求给定函数的最佳一致逼近多项式。
        /// </summary>
        /// <param name="a">区间左端点值</param>
        /// <param name="b">区间右端点值</param>
        /// <param name="p">前 n 个元素返回 n-1 次最佳一致逼近多项式 Pn-1(x) 的 n 个系数；
        /// 最后一个元素 p[n] 返回 Pn-1(x) 的偏差绝对值µ</param>
        /// <param name="n">n-1 次最佳一致逼近多项式的项数，即最佳一致逼近多项式的最高次为 n-1 。
        /// 要求 n≤20 。若 n>20 , 则本函数自动取 n = 20</param>
        /// <param name="eps">控制精度要求，一般在 10^(-10) ~ 10^(-35)之间</param>
        /// <param name="func">指向计算函数值 f(x) 的函数</param>
        public static void REMZ(double a, double b, out double[] p, int n, double eps, Func<double, double> func)
        {
            int i, j, k, m, nn;
            double[] x = new double[21], g = new double[21];
            double d, t, u, s, xx, x0, h, yy;
            p = new double[n + 1];

            if (n > 19) n = 19;
            nn = n + 1;
            m = nn + 1; d = 1.0e+35;
            for (k = 0; k <= nn; k++)
            {
                t = Math.Cos((nn - k) * 3.1415926 / (1.0 * nn));
                x[k] = (b + a + (b - a) * t) / 2.0;
            }
            while (1 == 1)
            {
                u = 1.0;
                for (i = 0; i <= m - 1; i++)
                {
                    p[i] = func(x[i]);
                    g[i] = -u; u = -u;
                }
                for (j = 0; j <= nn - 1; j++)
                {
                    k = m; s = p[k - 1]; xx = g[k - 1];
                    for (i = j; i <= nn - 1; i++)
                    {
                        t = p[nn - i + j - 1]; x0 = g[nn - i + j - 1];
                        p[k - 1] = (s - t) / (x[k - 1] - x[m - i - 2]);
                        g[k - 1] = (xx - x0) / (x[k - 1] - x[m - i - 2]);
                        k = nn - i + j; s = t; xx = x0;
                    }
                }
                u = -p[m - 1] / g[m - 1];
                for (i = 0; i <= m - 1; i++)
                    p[i] = p[i] + g[i] * u;
                for (j = 1; j <= nn - 1; j++)
                {
                    k = nn - j; h = x[k - 1]; s = p[k - 1];
                    for (i = m - j; i <= nn; i++)
                    {
                        t = p[i - 1]; p[k - 1] = s - h * t;
                        s = t; k = i;
                    }
                }
                p[m - 1] = Math.Abs(u); u = p[m - 1];
                if (Math.Abs(u - d) <= eps) return;
                d = u; h = 0.1 * (b - a) / (1.0 * nn);
                xx = a; x0 = a;
                while (x0 <= b)
                {
                    s = func(x0); t = p[nn - 1];
                    for (i = nn - 2; i >= 0; i--)
                        t = t * x0 + p[i];
                    s = Math.Abs(s - t);
                    if (s > u) { u = s; xx = x0; }
                    x0 = x0 + h;
                }
                s = func(xx); t = p[nn - 1];
                for (i = nn - 2; i >= 0; i--) t = t * xx + p[i];
                yy = s - t; i = 1; j = nn + 1;
                while ((j - i) != 1)
                {
                    k = (i + j) / 2;
                    if (xx < x[k - 1]) j = k;
                    else i = k;
                }
                if (xx < x[0])
                {
                    s = func(x[0]); t = p[nn - 1];
                    for (k = nn - 2; k >= 0; k--) t = t * x[0] + p[k];
                    s = s - t;
                    if (s * yy > 0.0) x[0] = xx;
                    else
                    {
                        for (k = nn - 1; k >= 0; k--) x[k + 1] = x[k];
                        x[0] = xx;
                    }
                }
                else
                {
                    if (xx > x[nn])
                    {
                        s = func(x[nn]); t = p[nn - 1];
                        for (k = nn - 2; k >= 0; k--) t = t * x[nn] + p[k];
                        s = s - t;
                        if (s * yy > 0.0) x[nn] = xx;
                        else
                        {
                            for (k = 0; k <= nn - 1; k++) x[k] = x[k + 1];
                            x[nn] = xx;
                        }
                    }
                    else
                    {
                        i = i - 1; j = j - 1;
                        s = func(x[i]); t = p[nn - 1];
                        for (k = nn - 2; k >= 0; k--) t = t * x[i] + p[k];
                        s = s - t;
                        if (s * yy > 0.0) x[i] = xx;
                        else x[j] = xx;
                    }
                }
            }
        }

        /// <summary>
        /// 用最小二乘法求矩形域上 n×m 个数据点的拟合曲面
        /// </summary>
        /// <param name="x">存放给定数据点的n个X坐标</param>
        /// <param name="y">存放给定数据点的m个Y坐标</param>
        /// <param name="z">存放矩形区域内 n×m 个网点上的函数值</param>
        /// <param name="n">X 坐标个数</param>
        /// <param name="m">Y 坐标个数</param>
        /// <param name="a">返回二元拟合多项式 f(x,y)=ΣΣa(x-x)^i(y-y)^j 的各系数</param>
        /// <param name="p">拟合多项式中 x 的最高次数加1。要求 p≤n 且 p≤20 ，若不满足这个条件，本函数自动取 p=min{n, 20}</param>
        /// <param name="q">拟合多项式中 y 的最高次数加1。要求 q≤m 且 q≤20 ，若不满足这个条件，本函数自动取 q=min{m, 20}</param>
        /// <param name="dt">dt[0] 返回拟合多项式与数据点误差的平方和，dt[1]返回拟合多项式与数据
        /// 点误差的绝对值之和，dt[2]返回拟合多项式与数据点误差绝对值的最大值</param>
        public static void PIR2(double[] x, double[] y, double[,] z, int n, int m, out double[,] a, int p, int q, out double[] dt)
        {
            int i, j, k, l, kk;
            double[] apx = new double[20], apy = new double[20], bx = new double[20], by = new double[20];
            double[,] u = new double[20, 20];
            double[] t = new double[20], t1 = new double[20], t2 = new double[20]; double xx, yy, d1, d2, g = 0, g1, g2;
            double x2, dd, y1, x1;

            double[,] v = new double[20, m];

            a = new double[p, q];
            dt = new double[3];

            for (i = 0; i < p; i++)
                for (j = 0; j < q; j++) a[i, j] = 0.0;
            if (p > n) p = n;
            if (p > 20) p = 20;
            if (q > m) q = m ;
            if (q > 20) q = 20;
            xx = 0.0;
            for (i = 0; i <= n - 1; i++)
                xx = xx + x[i] / (1.0 * n);
            yy = 0.0;
            for (i = 0; i <= m - 1; i++)
                yy = yy + y[i] / (1.0 * m);
            d1 = 1.0 * n; apx[0] = 0.0;
            for (i = 0; i <= n - 1; i++)
                apx[0] = apx[0] + x[i] - xx;
            apx[0] = apx[0] / d1;
            for (j = 0; j <= m - 1; j++)
            {
                v[0, j] = 0.0;
                for (i = 0; i <= n - 1; i++)
                    v[0, j] = v[0, j] + z[i, j];
                v[0, j] = v[0, j] / d1;
            }
            if (p > 1)
            {
                d2 = 0.0; apx[1] = 0.0;
                for (i = 0; i <= n - 1; i++)
                {
                    g = x[i] - xx - apx[0];
                    d2 = d2 + g * g;
                    apx[1] = apx[1] + (x[i] - xx) * g * g;
                }
                apx[1] = apx[1] / d2;
                bx[1] = d2 / d1;
                for (j = 0; j <= m - 1; j++)
                {
                    v[1, j] = 0.0;
                    for (i = 0; i <= n - 1; i++)
                    {
                        g = x[i] - xx - apx[0];
                        v[1, j] = v[1, j] + z[i, j] * g;
                    }
                    v[1, j] = v[1, j] / d2;
                }
                d1 = d2;
            }
            for (k = 2; k <= p-1; k++)
            {
                d2 = 0.0; apx[k] = 0.0;
                for (j = 0; j <= m - 1; j++) v[k, j] = 0.0;
                for (i = 0; i <= n - 1; i++)
                {
                    g1 = 1.0; g2 = x[i] - xx - apx[0];
                    for (j = 2; j <= k; j++)
                    {
                        g = (x[i] - xx - apx[j - 1]) * g2 - bx[j - 1] * g1;
                        g1 = g2; g2 = g;
                    }
                    d2 = d2 + g * g;
                    apx[k] = apx[k] + (x[i] - xx) * g * g;
                    for (j = 0; j <= m - 1; j++)
                        v[k, j] = v[k, j] + z[i, j] * g;
                }
                for (j = 0; j <= m - 1; j++)
                    v[k, j] = v[k, j] / d2;
                apx[k] = apx[k] / d2;
                bx[k] = d2 / d1;
                d1 = d2;
            }
            d1 = m; apy[0] = 0.0;
            for (i = 0; i <= m - 1; i++)
                apy[0] = apy[0] + y[i] - yy;
            apy[0] = apy[0] / d1;
            for (j = 0; j <= p-1; j++)
            {
                u[j, 0] = 0.0;
                for (i = 0; i <= m - 1; i++)
                    u[j, 0] = u[j, 0] + v[j, i];
                u[j, 0] = u[j, 0] / d1;
            }
            if (q > 1)
            {
                d2 = 0.0; apy[1] = 0.0;
                for (i = 0; i <= m - 1; i++)
                {
                    g = y[i] - yy - apy[0];
                    d2 = d2 + g * g;
                    apy[1] = apy[1] + (y[i] - yy) * g * g;
                }
                apy[1] = apy[1] / d2;
                by[1] = d2 / d1;
                for (j = 0; j <= p-1; j++)
                {
                    u[j, 1] = 0.0;
                    for (i = 0; i <= m - 1; i++)
                    {
                        g = y[i] - yy - apy[0];
                        u[j, 1] = u[j, 1] + v[j, i] * g;
                    }
                    u[j, 1] = u[j, 1] / d2;
                }
                d1 = d2;
            }
            for (k = 2; k <= q-1; k++)
            {
                d2 = 0.0; apy[k] = 0.0;
                for (j = 0; j <= p-1; j++) u[j, k] = 0.0;
                for (i = 0; i <= m - 1; i++)
                {
                    g1 = 1.0;
                    g2 = y[i] - yy - apy[0];
                    for (j = 2; j <= k; j++)
                    {
                        g = (y[i] - yy - apy[j - 1]) * g2 - by[j - 1] * g1;
                        g1 = g2; g2 = g;
                    }
                    d2 = d2 + g * g;
                    apy[k] = apy[k] + (y[i] - yy) * g * g;
                    for (j = 0; j <= p-1; j++)
                        u[j, k] = u[j, k] + v[j, i] * g;
                }
                for (j = 0; j <= p-1; j++)
                    u[j, k] = u[j, k] / d2;
                apy[k] = apy[k] / d2;
                by[k] = d2 / d1;
                d1 = d2;
            }
            v[0, 0] = 1.0; v[1, 0] = -apy[0]; v[1, 1] = 1.0;
            for (i = 0; i <= p-1; i++)
                for (j = 0; j <= q-1; j++)
                    a[i, j] = 0.0;
            for (i = 2; i <= q-1; i++)
            {
                v[i, i] = v[i - 1, i - 1];
                v[i, i - 1] = -apy[i - 1] * v[i - 1, i - 1] + v[i - 1, i - 2];
                if (i >= 3)
                    for (k = i - 2; k >= 1; k--)
                        v[i, k] = -apy[i - 1] * v[i - 1, k] +
                                v[i - 1, k - 1] - by[i - 1] * v[i - 2, k];
                v[i, 0] = -apy[i - 1] * v[i - 1, 0] - by[i - 1] * v[i - 2, 0];
            }
            for (i = 0; i <= p-1; i++)
            {
                if (i == 0) { t[0] = 1.0; t1[0] = 1.0; }
                else
                {
                    if (i == 1)
                    {
                        t[0] = -apx[0]; t[1] = 1.0;
                        t2[0] = t[0]; t2[1] = t[1];
                    }
                    else
                    {
                        t[i] = t2[i - 1];
                        t[i - 1] = -apx[i - 1] * t2[i - 1] + t2[i - 2];
                        if (i >= 3)
                            for (k = i - 2; k >= 1; k--)
                                t[k] = -apx[i - 1] * t2[k] + t2[k - 1]
                                     - bx[i - 1] * t1[k];
                        t[0] = -apx[i - 1] * t2[0] - bx[i - 1] * t1[0];
                        t2[i] = t[i];
                        for (k = i - 1; k >= 0; k--)
                        { t1[k] = t2[k]; t2[k] = t[k]; }
                    }
                }
                for (j = 0; j <= q-1; j++)
                    for (k = i; k >= 0; k--)
                        for (l = j; l >= 0; l--)
                            a[k, l] = a[k, l] + u[i, j] * t[k] * v[j, l];
            }
            dt[0] = 0.0; dt[1] = 0.0; dt[2] = 0.0;
            for (i = 0; i <= n - 1; i++)
            {
                x1 = x[i] - xx;
                for (j = 0; j <= m - 1; j++)
                {
                    y1 = y[j] - yy;
                    x2 = 1.0; dd = 0.0;
                    for (k = 0; k <= p-1; k++)
                    {
                        g = a[k, q-1];
                        for (kk = q - 2; kk >= 0; kk--)
                            g = g * y1 + a[k, kk];
                        g = g * x2; dd = dd + g; x2 = x2 * x1;
                    }
                    dd = dd - z[i, j];
                    if (Math.Abs(dd) > dt[2]) dt[2] = Math.Abs(dd);
                    dt[0] = dt[0] + dd * dd;
                    dt[1] = dt[1] + Math.Abs(dd);
                }
            }
        }
    }
}
