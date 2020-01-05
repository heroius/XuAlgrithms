using System;
using System.Collections.Generic;
using System.Text;

namespace Heroius.XuAlgrithms
{
    /// <summary>
    /// 数值积分
    /// </summary>
    public static class NumericalIntegration
    {
        /// <summary>
        /// 用变步长梯形求积法计算定积分T=∫f(x)dx
        /// </summary>
        /// <param name="a">积分下限</param>
        /// <param name="b">积分上限，要求b>a</param>
        /// <param name="eps">积分精度要求</param>
        /// <param name="func">指向计算被积函数f(x)值的函数</param>
        /// <returns>返回一个积分值</returns>
        public static double FFTS(double a, double b, double eps, Func<double, double> func)
        {
            int n, k;
            double fa, fb, h, t1, p, s, x, t = 0;
            fa = func(a); fb = func(b);
            n = 1; h = b - a;
            t1 = h * (fa + fb) / 2.0;
            p = eps + 1.0;
            while (p >= eps)
            {
                s = 0.0;
                for (k = 0; k <= n - 1; k++)
                {
                    x = a + (k + 0.5) * h;
                    s = s + func(x);
                }
                t = (t1 + h * s) / 2.0;
                p = Math.Abs(t1 - t);
                t1 = t; n = n + n; h = h / 2.0;
            }
            return t;
        }

        /// <summary>
        /// 用变步长辛卜生(Simpson)求积法计算定积分S=∫f(x)dx
        /// </summary>
        /// <param name="a">积分下限</param>
        /// <param name="b">积分上限，要求b>a</param>
        /// <param name="eps">积分精度要求</param>
        /// <param name="func">指向计算被积函数f(x)值的函数</param>
        /// <returns>返回一个积分值</returns>
        public static double SIMP(double a, double b, double eps, Func<double, double> func)
        {
            int n, k;
            double h, t1, t2, s1, s2 = 0, ep, p, x;
            n = 1; h = b - a;
            t1 = h * (func(a) + func(b)) / 2.0;
            s1 = t1;
            ep = eps + 1.0;
            while (ep >= eps)
            {
                p = 0.0;
                for (k = 0; k <= n - 1; k++)
                {
                    x = a + (k + 0.5) * h;
                    p = p + func(x);
                }
                t2 = (t1 + h * p) / 2.0;
                s2 = (4.0 * t2 - t1) / 3.0;
                ep = Math.Abs(s2 - s1);
                t1 = t2; s1 = s2; n = n + n; h = h / 2.0;
            }
            return s2;
        }

        /// <summary>
        /// 用自适应梯形求积法计算被积函数f(x)在积分区间内有强峰的定积分T=∫f(x)dx
        /// </summary>
        /// <param name="a">积分下限</param>
        /// <param name="b">积分上限，要求b>a</param>
        /// <param name="eps">积分精度要求</param>
        /// <param name="d">对积分区间进行分割的最小步长，当子区间的宽度小于d时，即使没有满足精度
        /// 要求，也不再往下进行分割</param>
        /// <param name="func">指向计算被积函数f(x)值的函数</param>
        /// <returns>返回一个积分值</returns>
        public static double FPTS(double a, double b, double eps, double d, Func<double, double> func)
        {
            double h, f0, f1, t0;
            double[] t = new double[2];
            h = b - a; t[0] = 0.0;
            f0 = func(a); f1 = func(b);
            t0 = h * (f0 + f1) / 2.0;
            FPTS_PPP(a, b, h, f0, f1, t0, eps, d, t, func);
            return t[0];
        }

        static void FPTS_PPP(double x0, double x1, double h, double f0, double f1, double t0, double eps, double d, double[] t, Func<double, double> func)
        {
            double x, f, t1, t2, p, g, eps1;
            x = x0 + h / 2.0; f = func(x);
            t1 = h * (f0 + f) / 4.0; t2 = h * (f + f1) / 4.0;
            p = Math.Abs(t0 - (t1 + t2));
            if ((p < eps) || (h / 2.0 < d))
            { t[0] = t[0] + (t1 + t2); return; }
            else
            {
                g = h / 2.0; eps1 = eps / 1.4;
                FPTS_PPP(x0, x, g, f0, f, t1, eps1, d, t, func);
                FPTS_PPP(x, x1, g, f, f1, t2, eps1, d, t, func);
                return;
            }
        }

        /// <summary>
        /// 用龙贝格(Romberg)求积法计算定积分T=∫f(x)dx
        /// </summary>
        /// <param name="a">积分下限</param>
        /// <param name="b">积分上限，要求b>a</param>
        /// <param name="eps">积分精度要求</param>
        /// <param name="func">指向计算被积函数f(x)值的函数</param>
        /// <returns>返回一个积分值</returns>
        public static double ROMB(double a, double b, double eps, Func<double, double> func)
        {
            int m, n, i, k;
            double h, ep, p, x, s, q = 0;
            double[] y = new double[10];
            h = b - a;
            y[0] = h * (func(a) + func(b)) / 2.0;
            m = 1; n = 1; ep = eps + 1.0;
            while ((ep >= eps) && (m <= 9))
            {
                p = 0.0;
                for (i = 0; i <= n - 1; i++)
                {
                    x = a + (i + 0.5) * h;
                    p = p + func(x);
                }
                p = (y[0] + h * p) / 2.0;
                s = 1.0;
                for (k = 1; k <= m; k++)
                {
                    s = 4.0 * s;
                    q = (s * p - y[k - 1]) / (s - 1.0);
                    y[k - 1] = p; p = q;
                }
                ep = Math.Abs(q - y[m - 1]);
                m = m + 1; y[m - 1] = q; n = n + n; h = h / 2.0;
            }
            return q;
        }

        /// <summary>
        /// 用连分式法计算定积分S=∫f(x)dx
        /// </summary>
        /// <param name="a">积分下限</param>
        /// <param name="b">积分上限，要求b>a</param>
        /// <param name="eps">积分精度要求</param>
        /// <param name="func">指向计算被积函数f(x)值的函数</param>
        /// <returns>返回一个积分值</returns>
        public static double FPQG(double a, double b, double eps, Func<double, double> func)
        {
            int m, n, k, l, j;
            double hh, t1, s1, ep, s, x, t2, g = 0;
            double[] h = new double[10], bb = new double[10];
            m = 1; n = 1;
            hh = b - a; h[0] = hh;
            t1 = hh * (func(a) + func(b)) / 2.0;
            s1 = t1; bb[0] = s1; ep = 1.0 + eps;
            while ((ep >= eps) && (m <= 9))
            {
                s = 0.0;
                for (k = 0; k <= n - 1; k++)
                {
                    x = a + (k + 0.5) * hh;
                    s = s + func(x);
                }
                t2 = (t1 + hh * s) / 2.0;
                m = m + 1;
                h[m - 1] = h[m - 2] / 2.0;
                g = t2;
                l = 0; j = 2;
                while ((l == 0) && (j <= m))
                {
                    s = g - bb[j - 2];
                    if (Math.Abs(s) + 1.0 == 1.0) l = 1;
                    else g = (h[m - 1] - h[j - 2]) / s;
                    j = j + 1;
                }
                bb[m - 1] = g;
                if (l != 0) bb[m - 1] = 1.0e+35;
                g = bb[m - 1];
                for (j = m; j >= 2; j--)
                    g = bb[j - 2] - h[j - 2] / g;
                ep = Math.Abs(g - s1);
                s1 = g; t1 = t2; hh = hh / 2.0; n = n + n;
            }
            return g;
        }

        /// <summary>
        /// 用分部积分法计算高振荡函数的积分∫f(x)sinmxdx与∫f(x)cosmxdx
        /// </summary>
        /// <param name="a">积分下限</param>
        /// <param name="b">积分上限，要求b>a</param>
        /// <param name="m">被积函数中振荡函数的角频率</param>
        /// <param name="n">给定积分区间两端点上f(x)的导数最高阶数+1</param>
        /// <param name="fa">存放f(x)在积分区间端点x=a处的各阶导数值。即fa(k)=f(a), k=0,...,n-1</param>
        /// <param name="fb">存放f(x)在积分区间端点x=b处的各阶导数值。即fb(k)=f(b), k=0,...,n-1</param>
        /// <param name="s">s[0]返回积分值∫f(x)cosmxdx；
        /// s[1]返回积分值∫f(x)sinmxdx</param>
        public static void PART(double a, double b, int m, int n, double[] fa, double[] fb, out double[] s)
        {
            int mm, k, j;
            s = new double[2];
            double[] sa = new double[4], sb = new double[4], ca = new double[4], cb = new double[4];
            double sma, smb, cma, cmb;
            sma = Math.Sin(m * a); smb = Math.Sin(m * b);
            cma = Math.Cos(m * a); cmb = Math.Cos(m * b);
            sa[0] = sma; sa[1] = cma; sa[2] = -sma; sa[3] = -cma;
            sb[0] = smb; sb[1] = cmb; sb[2] = -smb; sb[3] = -cmb;
            ca[0] = cma; ca[1] = -sma; ca[2] = -cma; ca[3] = sma;
            cb[0] = cmb; cb[1] = -smb; cb[2] = -cmb; cb[3] = smb;
            s[0] = 0.0; s[1] = 0.0;
            mm = 1;
            for (k = 0; k <= n - 1; k++)
            {
                j = k;
                while (j >= 4) j = j - 4;
                mm = mm * m;
                s[0] = s[0] + (fb[k] * sb[j] - fa[k] * sa[j]) / (1.0 * mm);
                s[1] = s[1] + (fb[k] * cb[j] - fa[k] * ca[j]) / (1.0 * mm);
            }
            s[1] = -s[1];
        }

        /// <summary>
        /// 用变步长勒让德-高斯(Legendre-Gauss)求积法计算定积分G=∫f(x)dx
        /// </summary>
        /// <param name="a">积分下限</param>
        /// <param name="b">积分上限，要求b>a</param>
        /// <param name="eps">积分精度要求</param>
        /// <param name="func">指向计算被积函数f(x)值的函数</param>
        /// <returns>返回一个积分值</returns>
        public static double LRGS(double a, double b, double eps, Func<double, double> func)
        {
            int m, i, j;
            double s, p, ep, h, aa, bb, w, x, g = 0;
            double[] t = new double[5] {-0.9061798459,-0.5384693101,0.0,
                         0.5384693101,0.9061798459};
            double[] c = new double[5] {0.2369268851,0.4786286705,0.5688888889,
                        0.4786286705,0.2369268851};
            m = 1;
            h = b - a; s = Math.Abs(0.001 * h);
            p = 1.0e+35; ep = eps + 1.0;
            while ((ep >= eps) && (Math.Abs(h) > s))
            {
                g = 0.0;
                for (i = 1; i <= m; i++)
                {
                    aa = a + (i - 1.0) * h; bb = a + i * h;
                    w = 0.0;
                    for (j = 0; j <= 4; j++)
                    {
                        x = ((bb - aa) * t[j] + (bb + aa)) / 2.0;
                        w = w + func(x) * c[j];
                    }
                    g = g + w;
                }
                g = g * h / 2.0;
                ep = Math.Abs(g - p) / (1.0 + Math.Abs(g));
                p = g; m = m + 1; h = (b - a) / m;
            }
            return g;
        }

        /// <summary>
        /// 用拉盖尔－高斯(Laguerre-Gauss) 求积公式计算半无限区间[0,∞]上的积分G=∫f(x)dx
        /// </summary>
        /// <param name="func">指向计算被积函数f(x)值的函数</param>
        /// <returns>返回一个积分值</returns>
        public static double LAGS(Func<double, double> func)
        {
            int i;
            double[] t = new double[5]{0.26355990,1.41340290,
                3.59642600,7.08580990,12.64080000};
            double[] c = new double[5] {0.6790941054,1.638487956,
                 2.769426772,4.315944000,7.104896230};
            double g = 0.0;
            for (i = 0; i <= 4; i++)
            {
                g = g + c[i] * func(t[i]);
            }
            return g;
        }

        /// <summary>
        /// 用埃尔米特－高斯(Hermite-Gauss)求积公式计算无限区间(-∞,∞]上的积分G=∫f(x)dx
        /// </summary>
        /// <param name="func">指向计算被积函数f(x)值的函数</param>
        /// <returns>返回一个积分值</returns>
        public static double HMGS(Func<double, double> func)
        {
            int i;
            double[] t = new double[5]{-2.02018200,-0.95857190,
                           0.0,0.95857190,2.02018200};
            double[] c = new double[5]{1.181469599,0.9865791417,
               0.9453089237,0.9865791417,1.181469599};
            double g = 0.0;
            for (i = 0; i <= 4; i++)
            {
                g = g + c[i] * func(t[i]);
            }
            return g;
        }

        /// <summary>
        /// 用变步长切比雪夫(Chebyshev)求积公式计算定积分S=∫f(x)dx
        /// </summary>
        /// <param name="a">积分下限</param>
        /// <param name="b">积分上限，要求b>a</param>
        /// <param name="eps">积分精度要求</param>
        /// <param name="func">指向计算被积函数f(x)值的函数</param>
        /// <returns>返回一个积分值</returns>
        public static double CBSV(double a, double b, double eps, Func<double, double> func)
        {
            int m, i, j;
            double h, d, p, ep, g = 0, aa, bb, s, x;
            double[] t = new double[5]{-0.8324975,-0.3745414,0.0,
                                  0.3745414,0.8324975};
            m = 1;
            h = b - a; d = Math.Abs(0.001 * h);
            p = 1.0e+35; ep = 1.0 + eps;
            while ((ep >= eps) && (Math.Abs(h) > d))
            {
                g = 0.0;
                for (i = 1; i <= m; i++)
                {
                    aa = a + (i - 1.0) * h; bb = a + i * h;
                    s = 0.0;
                    for (j = 0; j <= 4; j++)
                    {
                        x = ((bb - aa) * t[j] + (bb + aa)) / 2.0;
                        s = s + func(x);
                    }
                    g = g + s;
                }
                g = g * h / 5.0;
                ep = Math.Abs(g - p) / (1.0 + Math.Abs(g));
                p = g; m = m + 1; h = (b - a) / m;
            }
            return g;
        }

        /// <summary>
        /// 用蒙特卡罗(Monte Carlo) 法计算定积分S=∫f(x)dx
        /// </summary>
        /// <param name="a">积分下限</param>
        /// <param name="b">积分上限，要求b>a</param>
        /// <param name="func">指向计算被积函数f(x)值的函数</param>
        /// <returns>返回一个积分值</returns>
        public static double MTCL(double a, double b, Func<double, double> func)
        {
            int m;
            double x, r = 1.0, s = 0.0, d = 10000.0;
            for (m = 0; m <= 9999; m++)
            {
                x = a + (b - a) * Random.RND1(r);
                s = s + func(x) / d;
            }
            s = s * (b - a);
            return s;
        }

        /// <summary>
        /// 用变步长辛卜生(Simpson) 方法计算二重积分S=∫dx∫f(x,y)dy
        /// </summary>
        /// <param name="a">积分下限</param>
        /// <param name="b">积分上限，要求b>a</param>
        /// <param name="eps">积分精度要求</param>
        /// <param name="ss">指向计算上下限 y1(x) 与 y0(x) ( 要求y1(x)>y0(x) ）的函数</param>
        /// <param name="func">指向计算被积函数f(x,y)值的函数</param>
        /// <returns>返回一个积分值</returns>
        public static double SIM2(double a, double b, double eps, Func<double, double[]> ss, Func<double, double, double> func)
        {
            int n, j;
            double h, d, s1, s2, t1, x, t2, g, s = 0, s0, ep;
            n = 1; h = 0.5 * (b - a);
            d = Math.Abs((b - a) * 1.0e-06);
            s1 = SIM2_SIM1(a, eps, ss, func); s2 = SIM2_SIM1(b, eps, ss, func);
            t1 = h * (s1 + s2);
            s0 = 1.0e+35; ep = 1.0 + eps;
            while (((ep >= eps) && (Math.Abs(h) > d)) || (n < 16))
            {
                x = a - h; t2 = 0.5 * t1;
                for (j = 1; j <= n; j++)
                {
                    x = x + 2.0 * h;
                    g = SIM2_SIM1(x, eps, ss, func);
                    t2 = t2 + h * g;
                }
                s = (4.0 * t2 - t1) / 3.0;
                ep = Math.Abs(s - s0) / (1.0 + Math.Abs(s));
                n = n + n; s0 = s; t1 = t2; h = h * 0.5;
            }
            return s;
        }

        static double SIM2_SIM1(double x, double eps, Func<double, double[]> s, Func<double, double, double> func)   //一维积分函数
        {
            int n, i;
            double h, d, t1, yy, t2, g = 0, ep, g0;
            n = 1;
            double[] y = s(x);

            h = 0.5 * (y[1] - y[0]);
            d = Math.Abs(h * 2.0e-06);
            t1 = h * (func(x, y[0]) + func(x, y[1]));
            ep = 1.0 + eps; g0 = 1.0e+35;
            while (((ep >= eps) && (Math.Abs(h) > d)) || (n < 16))
            {
                yy = y[0] - h;
                t2 = 0.5 * t1;
                for (i = 1; i <= n; i++)
                {
                    yy = yy + 2.0 * h;
                    t2 = t2 + h * func(x, yy);
                }
                g = (4.0 * t2 - t1) / 3.0;
                ep = Math.Abs(g - g0) / (1.0 + Math.Abs(g));
                n = n + n; g0 = g; t1 = t2; h = 0.5 * h;
            }
            return g;
        }

        /// <summary>
        /// 用高斯(Gauss)方法计算n重积分 S=∫dX0∫dX1...∫f(X0,X1,...,Xn-1)dXn-1
        /// </summary>
        /// <param name="n">积分重数</param>
        /// <param name="js">js[k]表示第k层积分区间所划分的子区间个数</param>
        /// <param name="ss">指向计算各层积分上下限（要求所有的上限＞下限）的函数
        /// <para>第一参数：当前积分层 （0-based）</para>
        /// <para>第二参数：积分重数</para>
        /// <para>第三参数：点</para>
        /// <para>返回当前层的上限y[1]和下限y[0]</para>
        /// </param>
        /// <param name="func">指向计算被积函数f(X0,X1,...Xn-1)值的函数
        /// <para>第一参数：积分重数</para>
        /// <para>第二参数：点</para>
        /// <para>返回：值</para>
        /// </param>
        /// <returns>返回积分值</returns>
        public static double GAUS(int n, int[] js, Func<int, int, double[], double[]> ss, Func<int, double[], double> func)
        {
            int mm, j, k, q, l;
            double s, p;
            double[] t = new double[5] {-0.9061798459,-0.5384693101,0.0,
                                 0.5384693101,0.9061798459};
            double[] c = new double[5] {0.2369268851,0.4786286705,0.5688888889,
                                         0.4786286705,0.2369268851};
            double[] y, x = new double[n], a = new double[2 * (n + 1)], b = new double[n + 1];
            int[] _is = new int[2 * (n + 1)];
            mm = 1; l = 1;
            a[n] = 1.0; a[2 * n + 1] = 1.0;
            while (l == 1)
            {
                for (j = mm; j <= n; j++)
                {
                    y = ss(j - 1, n, x);
                    a[j - 1] = 0.5 * (y[1] - y[0]) / js[j - 1];
                    b[j - 1] = a[j - 1] + y[0];
                    x[j - 1] = a[j - 1] * t[0] + b[j - 1];
                    a[n + j] = 0.0;
                    _is[j - 1] = 1; _is[n + j] = 1;
                }
                j = n; q = 1;
                while (q == 1)
                {
                    k = _is[j - 1];
                    if (j == n) p = func(n, x);
                    else p = 1.0;
                    a[n + j] = a[n + j + 1] * a[j] * p * c[k - 1] + a[n + j];
                    _is[j - 1] = _is[j - 1] + 1;
                    if (_is[j - 1] > 5)
                        if (_is[n + j] >= js[j - 1])
                        {
                            j = j - 1; q = 1;
                            if (j == 0)
                            {
                                s = a[n + 1] * a[0];
                                return s;
                            }
                        }
                        else
                        {

                            _is[n + j] = _is[n + j] + 1;
                            b[j - 1] = b[j - 1] + a[j - 1] * 2.0;
                            _is[j - 1] = 1; k = _is[j - 1];
                            x[j - 1] = a[j - 1] * t[k - 1] + b[j - 1];
                            if (j == n) q = 1;
                            else q = 0;
                        }
                    else
                    {
                        k = _is[j - 1];
                        x[j - 1] = a[j - 1] * t[k - 1] + b[j - 1];
                        if (j == n) q = 1;
                        else q = 0;
                    }
                }
                mm = j + 1;
            }
            return 0.0;
        }

        /// <summary>
        /// 用连分式计算二重积分S=∫dx∫f(x,y)dy
        /// </summary>
        /// <param name="a">积分下限</param>
        /// <param name="b">积分上限，要求b>a</param>
        /// <param name="eps">积分精度要求</param>
        /// <param name="ss">指向计算上下限 y1(x) 与 y0(x) ( 要求y1(x)>y0(x) ）的函数
        /// <para>第一参数：点</para>
        /// <para>返回上限y[1]和下限y[0]</para>
        /// </param>
        /// <param name="func">指向计算被积函数f(x,y)值的函数</param>
        /// <returns>返回一个积分值</returns>
        public static double PQG2(double a, double b, double eps, Func<double, double[]> ss, Func<double, double, double> func)
        {
            int m, n, k, l, j;
            double[] bb = new double[10], h = new double[10];
            double hh, s1, s2, t1, t2, x, g, s0, ep, s = 0;
            m = 1; n = 1;
            hh = b - a; h[0] = hh;
            s1 = PQG2_PQG1(a, eps, ss, func); s2 = PQG2_PQG1(b, eps, ss, func);
            t1 = hh * (s1 + s2) / 2.0;
            s0 = t1; bb[0] = t1; ep = 1.0 + eps;
            while ((ep >= eps) && (m <= 9))
            {
                t2 = 0.5 * t1;
                for (k = 0; k <= n - 1; k++)
                {
                    x = a + (k + 0.5) * hh;
                    s1 = PQG2_PQG1(x, eps, ss, func);
                    t2 = t2 + 0.5 * s1 * hh;
                }
                m = m + 1;
                h[m - 1] = h[m - 2] / 2.0;
                g = t2; l = 0; j = 2;
                while ((l == 0) && (j <= m))
                {
                    s = g - bb[j - 2];
                    if (Math.Abs(s) + 1.0 == 1.0) l = 1;
                    else g = (h[m - 1] - h[j - 2]) / s;
                    j = j + 1;
                }
                bb[m - 1] = g;
                if (l != 0) bb[m - 1] = 1.0e+35;
                s = bb[m - 1];
                for (j = m; j >= 2; j--) s = bb[j - 2] - h[j - 2] / s;
                ep = Math.Abs(s - s0) / (1.0 + Math.Abs(s));
                n = n + n; t1 = t2; s0 = s; hh = hh / 2.0;
            }
            return s;
        }

        static double PQG2_PQG1(double x, double eps, Func<double, double[]> ss, Func<double, double, double> func)
        {
            int m, n, k, l, j;
            double[] b = new double[10], h = new double[10];
            double hh, t1, t2, s0, yy, g, ep, s = 0;
            m = 1; n = 1;
            double[] y = ss(x);
            hh = y[1] - y[0]; h[0] = hh;
            t1 = 0.5 * hh * (func(x, y[0]) + func(x, y[1]));
            s0 = t1; b[0] = t1; ep = 1.0 + eps;
            while ((ep >= eps) && (m <= 9))
            {
                t2 = 0.5 * t1;
                for (k = 0; k <= n - 1; k++)
                {
                    yy = y[0] + (k + 0.5) * hh;
                    t2 = t2 + 0.5 * hh * func(x, yy);
                }
                m = m + 1;
                h[m - 1] = h[m - 2] / 2.0;
                g = t2; l = 0; j = 2;
                while ((l == 0) && (j <= m))
                {
                    s = g - b[j - 2];
                    if (Math.Abs(s) + 1.0 == 1.0) l = 1;
                    else g = (h[m - 1] - h[j - 2]) / s;
                    j = j + 1;
                }
                b[m - 1] = g;
                if (l != 0) b[m - 1] = 1.0e+35;
                s = b[m - 1];
                for (j = m; j >= 2; j--) s = b[j - 2] - h[j - 2] / s;
                ep = Math.Abs(s - s0) / (1.0 + Math.Abs(s));
                n = n + n; t1 = t2; s0 = s; hh = 0.5 * hh;
            }
            return s;
        }

        /// <summary>
        /// 用蒙特卡罗(Monte Carlo) 法计算多重积分 S=∫∫...∫f(X0,X1,...,Xn-1)dX0dX1...dXn-1
        /// </summary>
        /// <param name="n">积分的重数</param>
        /// <param name="a">存放各层积分的下限值</param>
        /// <param name="b">存放各层积分的上限值</param>
        /// <param name="func">指向计算被积函数值f(X0,X1,...,Xn-1)的函数
        /// <para>第一参数：积分重数</para>
        /// <para>第二参数：点</para>
        /// <para>返回：值</para>
        /// </param>
        /// <returns>函数返回积分值</returns>
        public static double MTML(int n, double[] a, double[] b, Func<int, double[], double> func)
        {
            int m, i;
            double r = 1.0, d = 10000.0, s = 0.0;
            double[] x = new double[n];

            for (m = 0; m <= 9999; m++)
            {
                for (i = 0; i <= n - 1; i++)
                    x[i] = a[i] + (b[i] - a[i]) * Random.RND1(r);
                s = s + func(n, x) / d;
            }
            for (i = 0; i <= n - 1; i++) s = s * (b[i] - a[i]);
            return s;
        }
    }
}
