using System;
using System.Collections.Generic;
using System.Text;

namespace Heroius.XuAlgrithms
{
    /// <summary>
    /// 非线性方程与方程组的求解
    /// </summary>
    public class NonlinearEquations
    {
        /// <summary>
        /// 用对分法搜索方程f(x)在区间[a,b]内的实根，返回实际搜索到的实根个数
        /// </summary>
        /// <param name="a">求根区间的左端点</param>
        /// <param name="b">求根区间的右端点</param>
        /// <param name="h">搜索求根时采用的步长</param>
        /// <param name="eps">控制精度要求</param>
        /// <param name="m">在区间[a,b]内实根个数的预估值</param>
        /// <param name="func">计算方程左端函数值f(x)的函数</param>
        /// <returns>在区间[a,b]内实际搜索到的的实根个数</returns>
        public static int DHRT(double a, double b, double h, double eps, int m, Func<double, double> func, out double[] x)
        {
            int js;
            double z, y, z1, y1, z0, y0;
            int n = 0; z = a; y = func(z);
            x = new double[m];
            while ((z <= b + h / 2.0) && (n != m))
            {
                if (Math.Abs(y) < eps)
                {
                    n = n + 1; x[n - 1] = z;
                    z = z + h / 2.0; y = func(z);
                }
                else
                {
                    z1 = z + h; y1 = func(z1);
                    if (Math.Abs(y1) < eps)
                    {
                        n = n + 1; x[n - 1] = z1;
                        z = z1 + h / 2.0; y = func(z);
                    }
                    else if (y * y1 > 0.0)
                    { y = y1; z = z1; }
                    else
                    {
                        js = 0;
                        while (js == 0)
                        {
                            if (Math.Abs(z1 - z) < eps)
                            {
                                n = n + 1; x[n - 1] = (z1 + z) / 2.0;
                                z = z1 + h / 2.0; y = func(z);
                                js = 1;
                            }
                            else
                            {
                                z0 = (z1 + z) / 2.0; y0 = func(z0);
                                if (Math.Abs(y0) < eps)
                                {
                                    x[n] = z0; n = n + 1; js = 1;
                                    z = z0 + h / 2.0; y = func(z);
                                }
                                else if ((y * y0) < 0.0)
                                { z1 = z0; y1 = y0; }
                                else { z = z0; y = y0; }
                            }
                        }
                    }
                }
            }
            return n;
        }

        /// <summary>
        /// 用牛顿 (Newton) 迭代法求方程 f(x) =O 的一个实根。
        /// </summary>
        /// <param name="x">迭代初值。返回时指向迭代终值</param>
        /// <param name="eps">控制精度要求</param>
        /// <param name="js">最大迭代次数</param>
        /// <param name="func">指向计算f(x)与f'(x)值的函数：输入x，返回y和dy</param>
        /// <returns>函数返回标志值。若返回的标志值小于0, 则表示出现了f'(x)=O的情况；
        /// 若返回的标志值等于最大迭代次数js, 则说明返回的终值有可能不满足精度要求；
        /// 若返回的标志值大于等于0且小于最大迭代次数，则表示正常返回</returns>
        public static int NEWT(ref double x, double eps, int js, Func<double, Tuple<double, double>> func)
        {
            int k, l;
            double d, p, x0, x1 = 0;
            l = js; k = 1; x0 = x;
            var y = func(x0);
            d = eps + 1.0;
            while ((d >= eps) && (l != 0))
            {
                if (Math.Abs(y.Item1) + 1.0 == 1.0)
                {
                    throw new Exception("f(x)的导数为0！");
                }
                x1 = x0 - y.Item1 / y.Item2;
                y = func(x1);
                d = Math.Abs(x1 - x0); p = Math.Abs(y.Item1);
                if (p > d) d = p;
                x0 = x1; l = l - 1;
            }
            x = x1;
            k = js - l;
            return k;
        }

        /// <summary>
        /// 用埃特金 (Aitken) 迭代法求非线性方程 x=φ(x) 的一个实根。
        /// </summary>
        /// <param name="x">指向迭代初值，返回时指向迭代终值</param>
        /// <param name="eps">控制精度要求</param>
        /// <param name="js">最大迭代次数</param>
        /// <param name="func">指向计算φ(x)值的函数</param>
        /// <returns>函数返回标志值。
        /// 若返回的标志值为0，则表示实际迭代次数已等于最大迭代次数js, 说明返回的终值有可能不满足精度要求；
        /// 若返回的标志值大于0, 则表示正常返回</returns>
        public static int ATKN(ref double x, double eps, int js, Func<double, double> func)
        {
            int flag, l;
            double u, v, x0;
            l = 0; x0 = x; flag = 0;
            while ((flag == 0) && (l != js))
            {
                l++;
                u = func(x0);
                v = func(u);
                if (Math.Abs(u - v) < eps)
                {
                    x0 = v; flag = 1;
                }
                else
                {
                    x0 = v - (v - u) * (v - u) / (v - 2.0 * u + x0);
                }
            }
            x = x0; l = js - l;
            return l;
        }

        /// <summary>
        /// 利用试位法求非线性方程 f(x) = O 在给定区间 [ a, b ] 内的一个实根
        /// </summary>
        /// <param name="a">求根区间的左端点</param>
        /// <param name="b">求根区间的右端点</param>
        /// <param name="eps">控制精度要求</param>
        /// <param name="func">指向计算方程左端函数 f( x) 值的函数</param>
        /// <param name="x">指向方程实根</param>
        /// <returns>函数返回标志值。若返回的标志值为-1, 则 表示 f(a)f(b)>O;
        /// 若返回的标志值大于 0, 则表示迭代次数</returns>
        public static int FALSE(double a, double b, double eps, Func<double, double> func, out double x)
        {
            int m = 0;
            double fa = func(a), fb = func(b), y;
            if (fa * fb > 0)
            {
                x = 0;
                return -1;
            }
            do
            {
                m++;
                x = (a * fb - b * fa) / (fb - fa);
                y = func(x);
                if (y * fa < 0)
                {
                    b = x; fb = y;
                }
                else
                {
                    a = x; fa = y;
                }
            } while (Math.Abs(y) >= eps);
            return m;
        }

        /// <summary>
        /// 利用连分式求非线性方程 f(x) =O 的一个实根。
        /// </summary>
        /// <param name="x">指向迭代初值，返回时指向迭代终值</param>
        /// <param name="eps">控制精度要求</param>
        /// <param name="func">指向计算 f(x) 值的函数</param>
        /// <returns>函数返回迭代次数。
        /// 若返回的迭代次数等于10, 返回的迭代终值有可能没有满足精度要求；
        /// 若返回的迭代次数小于10, 则表示正常返回</returns>
        public static int PQRT(ref double x, double eps, Func<double, double> func)
        {
            int i, j, m, it = 0, l = 10;
            double[] a = new double[10], y = new double[10];
            double z, h = 0.0, x0 = x, q = 1.0e+35;
            while (l != 0)
            {
                l--;
                j = 0; it = l;
                while (j <= 7)
                {
                    if (j <= 2)
                    {
                        z = x0 + 0.1 * j;
                    }
                    else
                    {
                        z = h;
                    }
                    y[j] = func(z);
                    h = z;
                    if (j == 0)
                    {
                        a[0] = z;
                    }
                    else
                    {
                        m = 0; i = 0;
                        while ((m == 0) && (i <= j - 1))
                        {
                            if (Math.Abs(h - a[i]) == 0)
                            {
                                m = 1;
                            }
                            else
                            {
                                h = (y[j] - y[i]) / (h - a[i]);
                            }
                            i++;
                        }
                        a[j] = h;
                        if (m != 0) a[j] = q;
                        h = 0.0;
                        for (i = j - 1; i >= 0; i--)
                        {
                            if (Math.Abs(a[i + 1] + h) == 0.0)
                            {
                                h = q;
                            }
                            else
                            {
                                h = -y[i] / (a[i + 1] + h);
                            }
                        }
                        h = h + a[0];
                    }
                    if (Math.Abs(y[j]) >= eps)
                    {
                        j++;
                    }
                    else
                    {
                        j = 10; l = 0;
                    }
                }
                x0 = h;
            }
            x = h;
            return 10 - it;
        }

        /// <summary>
        /// 用 QR 方法求实系数 n 次多项式方程的全部实根与复根
        /// </summary>
        /// <param name="a">存放 n 次多项式方程的 n+l 个系数</param>
        /// <param name="n">多项式方程的次数</param>
        /// <param name="eps">QR 方法中的控制精度要求</param>
        /// <param name="jt">控制 QR 方法中最大迭代次数</param>
        /// <param name="xr">返回 n 个根的实部</param>
        /// <param name="xi">返回 n 个根的虚部</param>
        /// <returns>函数返回标志值。
        /// 若返回的标志值小于0, 则表示在QR方法中没有满足精度要求；
        /// 若返回的标志值大于0, 则表示正常返回</returns>
        public static void QRRT(double[] a, int n, double eps, int jt, out double[] xr, out double[] xi)
        {
            if (a.Length != n + 1)
            {
                throw new Exception("系数长度应为n+1");
            }
            int i, j;
            double[] q = new double[n * n];
            for (j = 0; j < n - 1; j++)
            {
                q[j] = -a[n - j - 1] / a[n];
            }
            for (j = n; j <= n * n - 1; j++)
            {
                q[j] = 0.0;
            }
            for (i = 0; i <= n - 2; i++)
            {
                q[(i + 1) * n + i] = 1.0;
            }
            MatrixEigen.HHQR(Utility.C.Convert(q, n, n), eps, jt, out xr, out xi);
        }

        /// <summary>
        /// 用牛顿下山法求实系数代数方程的全部根
        /// </summary>
        /// <param name="a">存放 n 次多项式方程的实系数</param>
        /// <param name="n">多项式方程的次数</param>
        /// <param name="xr">返回 n 个根的实部</param>
        /// <param name="xi">返回 n 个根的虚部</param>
        /// <returns>函数返回标志值。
        /// 若返回的标志值小于0,则表示方程中所有系数为0; 
        /// 若返回的标志值大于0, 则表示正常返回</returns>
        public static int SRRT(double[] a, int n, out double[] xr, out double[] xi)
        {
            if (a.Length != n + 1)
            {
                throw new Exception("系数长度应为n+1");
            }
            int m = 0, i = 0, jt = 0, k = 0, _is = 0, it = 0;
            double t = 0, x = 0, y = 0, x1 = 0, y1 = 0, dx = 0, dy = 0, p = 0, q = 0, w = 0, dd = 0, dc = 0, c = 0, g = 0, u = 0, v = 0, pq = 0, g1 = 0, u1 = 0, v1 = 0;
            m = n;
            xr = new double[n];
            xi = new double[n];
            while ((m > 0) && (Math.Abs(a[m]) + 1.0 == 1.0))
            {
                m = m - 1;
            }
            if (m <= 0)
            {
                return -1; // fail
            }
            for (i = 0; i <= m / 2; i++)
            {
                w = a[i];
                a[i] = a[m - i];
                a[m - i] = w;
            }
            k = m;
            _is = 0;
            w = 1.0;
            jt = 1;
            while (jt == 1)
            {
                pq = Math.Abs(a[k]);
                while (pq < 1.0e-12)
                {
                    xr[k - 1] = 0.0;
                    xi[k - 1] = 0.0;
                    k = k - 1;
                    if (k == 1)
                    {
                        xr[0] = -a[1] * w / a[0];
                        xi[0] = 0.0;
                        return 1;
                    }
                    pq = Math.Abs(a[k]);
                }
                q = Math.Log(pq);
                q = q / (1.0 * k);
                q = Math.Exp(q);
                p = q;
                w = w * p;
                for (i = 1; i <= k; i++)
                {
                    a[i] = a[i] / q;
                    q = q * p;
                }
                x = 0.0001;
                x1 = x;
                y = 0.2;
                y1 = y;
                dx = 1.0;
                g = 1.0e37;
            l40:
                u = a[0];
                v = 0.0;
                for (i = 1; i <= k; i++)
                {
                    p = u * x1;
                    q = v * y1;
                    pq = (u + v) * (x1 + y1);
                    u = p - q + a[i];
                    v = pq - p - q;
                }
                g1 = u * u + v * v;
                if (g1 >= g)
                {
                    if (_is != 0)
                    {
                        it = 1;
                        SRRT_G65(ref x, ref y, ref x1, ref y1, ref dx, ref dy, ref dd, ref dc, ref c, ref k, ref _is, ref it);
                        if (it == 0) goto l40;
                    }
                    else
                    {
                        SRRT_G60(ref t, ref x, ref y, ref x1, ref y1, ref dx, ref dy, ref p, ref q, ref k, ref it);
                        if (t >= 1.0e-3) goto l40;
                        if (g > 1.0e-18)
                        {
                            it = 0;
                            SRRT_G65(ref x, ref y, ref x1, ref y1, ref dx, ref dy, ref dd, ref dc, ref c, ref k, ref _is, ref it);
                            if (it == 0) goto l40;
                        }
                    }
                    SRRT_G90(xr, xi, a, ref x, ref y, ref p, ref q, ref w, ref k);
                }
                else
                {
                    g = g1;
                    x = x1;
                    y = y1;
                    _is = 0;
                    if (g <= 1.0e-22)
                        SRRT_G90(xr, xi, a, ref x, ref y, ref p, ref q, ref w, ref k);
                    else
                    {
                        u1 = k * a[0];
                        v1 = 0.0;
                        for (i = 2; i <= k; i++)
                        {
                            p = u1 * x;
                            q = v1 * y;
                            pq = (u1 + v1) * (x + y);
                            u1 = p - q + (k - i + 1) * a[i - 1];
                            v1 = pq - p - q;
                        }
                        p = u1 * u1 + v1 * v1;
                        if (p <= 1.0e-20)
                        {
                            it = 0;
                            SRRT_G65(ref x, ref y, ref x1, ref y1, ref dx, ref dy, ref dd, ref dc, ref c, ref k, ref _is, ref it);
                            if (it == 0) goto l40;
                            SRRT_G90(xr, xi, a, ref x, ref y, ref p, ref q, ref w, ref k);
                        }
                        else
                        {
                            dx = (u * u1 + v * v1) / p;
                            dy = (u1 * v - v1 * u) / p;
                            t = 1.0 + 4.0 / k;
                            SRRT_G60(ref t, ref x, ref y, ref x1, ref y1, ref dx, ref dy, ref p, ref q, ref k, ref it);
                            if (t >= 1.0e-3) goto l40;
                            if (g > 1.0e-18)
                            {
                                it = 0;
                                SRRT_G65(ref x, ref y, ref x1, ref y1, ref dx, ref dy, ref dd, ref dc, ref c, ref k, ref _is, ref it);
                                if (it == 0) goto l40;
                            }
                            SRRT_G90(xr, xi, a, ref x, ref y, ref p, ref q, ref w, ref k);
                        }
                    }
                }
                if (k == 1) jt = 0;
                else jt = 1;
            }
            return 1;
        }

        #region SRRT sub func

        /// <summary>
        /// SRRT sub func
        /// </summary>
        /// <param name="t"></param>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="x1"></param>
        /// <param name="y1"></param>
        /// <param name="dx"></param>
        /// <param name="dy"></param>
        /// <param name="p"></param>
        /// <param name="q"></param>
        /// <param name="k"></param>
        /// <param name="it"></param>
        static void SRRT_G60(ref double t, ref double x, ref double y, ref double x1, ref double y1, ref double dx, ref double dy, ref double p, ref double q, ref int k, ref int it)
        {
            it = 1;
            while (it == 1)
            {
                t = t / 1.67;
                it = 0;
                x1 = x - t * dx;
                y1 = y - t * dy;
                if (k >= 50)
                {
                    p = Math.Sqrt(x1 * x1 + y1 * y1);
                    q = Math.Exp(85.0 / k);
                    if (p >= q) it = 1;
                }
            }
        }
        /// <summary>
        /// SRRT sub func
        /// </summary>
        /// <param name="xr"></param>
        /// <param name="xi"></param>
        /// <param name="a"></param>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="p"></param>
        /// <param name="q"></param>
        /// <param name="w"></param>
        /// <param name="k"></param>
        static void SRRT_G90(double[] xr, double[] xi, double[] a, ref double x, ref double y, ref double p, ref double q, ref double w, ref int k)
        {
            int i;
            if (Math.Abs(y) <= 1.0e-6)
            {
                p = -x;
                y = 0.0;
                q = 0.0;
            }
            else
            {
                p = -2.0 * x;
                q = x * x + y * y;
                xr[k - 1] = x * w;
                xi[k - 1] = -y * w;
                k = k - 1;
            }
            for (i = 1; i <= k; i++)
            {
                a[i] = a[i] - a[i - 1] * p;
                a[i + 1] = a[i + 1] - a[i - 1] * q;
            }
            xr[k - 1] = x * w;
            xi[k - 1] = y * w;
            k = k - 1;
            if (k == 1)
            {
                xr[0] = -a[1] * w / a[0];
                xi[0] = 0.0;
            }
        }
        /// <summary>
        /// SRRT sub func
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="x1"></param>
        /// <param name="y1"></param>
        /// <param name="dx"></param>
        /// <param name="dy"></param>
        /// <param name="dd"></param>
        /// <param name="dc"></param>
        /// <param name="c"></param>
        /// <param name="k"></param>
        /// <param name="_is"></param>
        /// <param name="it"></param>
        static void SRRT_G65(ref double x, ref double y, ref double x1, ref double y1, ref double dx, ref double dy, ref double dd, ref double dc, ref double c, ref int k, ref int _is, ref int it)
        {
            if (it == 0)
            {
                _is = 1;
                dd = Math.Sqrt(dx * dx + dy * dy);
                if (dd > 1.0) dd = 1.0;
                dc = 6.28 / (4.5 * k);
                c = 0.0;
            }
            while (true)
            {
                c += dc;
                dx = dd * Math.Cos(c);
                dy = dd * Math.Sin(c);
                x1 = x + dx;
                y1 = y + dy;
                if (c <= 6.29)
                {
                    it = 0;
                    return;
                }
                dd = dd / 1.67;
                if (dd <= 1.0e-7)
                {
                    it = 1;
                    return;
                }
                c = 0.0;
            }
        }

        #endregion

        /// <summary>
        /// 用牛顿下山法求复系数代数方程的全部根
        /// </summary>
        /// <param name="ar">存放 n 次多项式方程系数的实部</param>
        /// <param name="ai">存放 n 次多项式方程系数的虚部</param>
        /// <param name="n">多项式方程的次数</param>
        /// <param name="xr">返回 n 个根的实部</param>
        /// <param name="xi">返回 n 个根的虚部</param>
        /// <returns>函数返回标志值。
        /// 若返回的标志值小于0,则表示方程中所有系数为0; 
        /// 若返回的标志值大于0, 则表示正常返回</returns>
        public static int CSRT(double[] ar, double[] ai, int n, out double[] xr, out double[] xi)
        {
            if (ar.Length != n + 1 || ai.Length != n + 1)
            {
                throw new Exception("系数长度应为n+1");
            }
            int m = 0, i = 0, jt = 0, k = 0, _is = 0, it = 0;
            double t = 0, x = 0, y = 0, x1 = 0, y1 = 0, dx = 0, dy = 0, p = 0, q = 0, w = 0, dd = 0, dc = 0, c = 0, g = 0, u = 0, v = 0, pq = 0, g1 = 0, u1 = 0, v1 = 0;
            m = n;
            xr = new double[n];
            xi = new double[n];
            p = Math.Sqrt(ar[m] * ar[m] + ai[m] * ai[m]);
            while ((m > 0) && (p + 1.0 == 1.0))
            {
                m = m - 1;
                p = Math.Sqrt(ar[m] * ar[m] + ai[m] * ai[m]);
            }
            if (m <= 0)
            {
                return -1; // fail
            }
            for (i = 0; i <= m; i++)
            {
                ar[i] = ar[i] / p;
                ai[i] = ai[i] / p;
            }
            for (i = 0; i <= m / 2; i++)
            {
                w = ar[i];
                ar[i] = ar[m - i];
                ar[m - i] = w;
                w = ai[i];
                ai[i] = ai[m - i];
                ai[m - i] = w;
            }
            k = m;
            _is = 0;
            w = 1.0;
            jt = 1;
            while (jt == 1)
            {
                pq = Math.Sqrt(ar[k] * ar[k] + ai[k] * ai[k]);
                while (pq < 1.0e-12)
                {
                    xr[k - 1] = 0.0;
                    xi[k - 1] = 0.0;
                    k = k - 1;
                    if (k == 1)
                    {
                        p = ar[0] * ar[0] + ai[0] * ai[0];
                        xr[0] = -w * (ar[0] * ar[1] + ai[0] * ai[1]) / p;
                        xi[0] = w * (ar[1] * ai[0] - ar[0] * ai[1]) / p;
                        return 1;
                    }
                    pq = Math.Sqrt(ar[k] * ar[k] + ai[k] * ai[k]);
                }
                q = Math.Log(pq);
                q = q / (1.0 * k);
                q = Math.Exp(q);
                p = q;
                w = w * p;
                for (i = 1; i <= k; i++)
                {
                    ar[i] = ar[i] / q;
                    ai[i] = ai[i] / q;
                    q = q * p;
                }
                x = 0.0001;
                x1 = x;
                y = 0.2;
                y1 = y;
                dx = 1.0;
                g = 1.0e37;
            l40:
                u = ar[0];
                v = ai[0];
                for (i = 1; i <= k; i++)
                {
                    p = u * x1;
                    q = v * y1;
                    pq = (u + v) * (x1 + y1);
                    u = p - q + ar[i];
                    v = pq - p - q + ai[i];
                }
                g1 = u * u + v * v;
                if (g1 >= g)
                {
                    if (_is != 0)
                    {
                        it = 1;
                        CSRT_G65(ref x, ref y, ref x1, ref y1, ref dx, ref dy, ref dd, ref dc, ref c, ref k, ref _is, ref it);
                        if (it == 0) goto l40;
                    }
                    else
                    {
                        CSRT_G60(ref t, ref x, ref y, ref x1, ref y1, ref dx, ref dy, ref p, ref q, ref k, ref it);
                        if (t >= 1.0e-3) goto l40;
                        if (g > 1.0e-18)
                        {
                            it = 0;
                            CSRT_G65(ref x, ref y, ref x1, ref y1, ref dx, ref dy, ref dd, ref dc, ref c, ref k, ref _is, ref it);
                            if (it == 0) goto l40;
                        }
                    }
                    CSRT_G90(xr, xi, ar, ai, ref x, ref y, ref p, ref q, ref w, ref k);
                }
                else
                {
                    g = g1;
                    x = x1;
                    y = y1;
                    _is = 0;
                    if (g <= 1.0e-22)
                        CSRT_G90(xr, xi, ar, ai, ref x, ref y, ref p, ref q, ref w, ref k);
                    else
                    {
                        u1 = k * ar[0];
                        v1 = ai[0];
                        for (i = 2; i <= k; i++)
                        {
                            p = u1 * x;
                            q = v1 * y;
                            pq = (u1 + v1) * (x + y);
                            u1 = p - q + (k - i + 1) * ar[i - 1];
                            v1 = pq - p - q + (k - i + 1) * ai[i - 1];
                        }
                        p = u1 * u1 + v1 * v1;
                        if (p <= 1.0e-20)
                        {
                            it = 0;
                            CSRT_G65(ref x, ref y, ref x1, ref y1, ref dx, ref dy, ref dd, ref dc, ref c, ref k, ref _is, ref it);
                            if (it == 0) goto l40;
                            CSRT_G90(xr, xi, ar, ai, ref x, ref y, ref p, ref q, ref w, ref k);
                        }
                        else
                        {
                            dx = (u * u1 + v * v1) / p;
                            dy = (u1 * v - v1 * u) / p;
                            t = 1.0 + 4.0 / k;
                            CSRT_G60(ref t, ref x, ref y, ref x1, ref y1, ref dx, ref dy, ref p, ref q, ref k, ref it);
                            if (t >= 1.0e-3) goto l40;
                            if (g > 1.0e-18)
                            {
                                it = 0;
                                CSRT_G65(ref x, ref y, ref x1, ref y1, ref dx, ref dy, ref dd, ref dc, ref c, ref k, ref _is, ref it);
                                if (it == 0) goto l40;
                            }
                            CSRT_G90(xr, xi, ar, ai, ref x, ref y, ref p, ref q, ref w, ref k);
                        }
                    }
                }
                if (k == 1) jt = 0;
                else jt = 1;
            }
            return 1;
        }

        #region CSRT sub func

        /// <summary>
        /// CSRT sub func
        /// </summary>
        /// <param name="t"></param>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="x1"></param>
        /// <param name="y1"></param>
        /// <param name="dx"></param>
        /// <param name="dy"></param>
        /// <param name="p"></param>
        /// <param name="q"></param>
        /// <param name="k"></param>
        /// <param name="it"></param>
        static void CSRT_G60(ref double t, ref double x, ref double y, ref double x1, ref double y1, ref double dx, ref double dy, ref double p, ref double q, ref int k, ref int it)
        {
            it = 1;
            while (it == 1)
            {
                t = t / 1.67;
                it = 0;
                x1 = x - t * dx;
                y1 = y - t * dy;
                if (k >= 30)
                {
                    p = Math.Sqrt(x1 * x1 + y1 * y1);
                    q = Math.Exp(75.0 / k);
                    if (p >= q) it = 1;
                }
            }
        }
        /// <summary>
        /// CSRT sub func
        /// </summary>
        /// <param name="xr"></param>
        /// <param name="xi"></param>
        /// <param name="ar"></param>
        /// <param name="ai"></param>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="p"></param>
        /// <param name="q"></param>
        /// <param name="w"></param>
        /// <param name="k"></param>
        static void CSRT_G90(double[] xr, double[] xi, double[] ar, double[] ai, ref double x, ref double y, ref double p, ref double q, ref double w, ref int k)
        {
            int i;
            for (i = 1; i <= k; i++)
            {
                ar[i] = ar[i] + ar[i - 1] * x - ai[i - 1] * y;
                ai[i] = ai[i] + ar[i - 1] * y + ai[i - 1] * x;
            }
            xr[k - 1] = x * w;
            xi[k - 1] = y * w;
            k--;
            if (k == 1)
            {
                p = ar[0] * ar[0] + ai[0] * ai[0];
                xr[0] = -w * (ar[0] * ar[1] + ai[0] * ai[1]) / p;
                xi[0] = w * (ar[1] * ai[0] - ar[0] * ai[1]) / p;
            }
        }
        /// <summary>
        /// CSRT sub func
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="x1"></param>
        /// <param name="y1"></param>
        /// <param name="dx"></param>
        /// <param name="dy"></param>
        /// <param name="dd"></param>
        /// <param name="dc"></param>
        /// <param name="c"></param>
        /// <param name="k"></param>
        /// <param name="_is"></param>
        /// <param name="it"></param>
        static void CSRT_G65(ref double x, ref double y, ref double x1, ref double y1, ref double dx, ref double dy, ref double dd, ref double dc, ref double c, ref int k, ref int _is, ref int it)
        {
            if (it == 0)
            {
                _is = 1;
                dd = Math.Sqrt(dx * dx + dy * dy);
                if (dd > 1.0) dd = 1.0;
                dc = 6.28 / (4.5 * k);
                c = 0.0;
            }
            while (true)
            {
                c += dc;
                dx = dd * Math.Cos(c);
                dy = dd * Math.Sin(c);
                x1 = x + dx;
                y1 = y + dy;
                if (c <= 6.29)
                {
                    it = 0;
                    return;
                }
                dd = dd / 1.67;
                if (dd <= 1.0e-7)
                {
                    it = 1;
                    return;
                }
                c = 0.0;
            }
        }

        #endregion

        /// <summary>
        /// 用梯度法（即最速下降法）求非线性方程组的一组实数解。
        /// </summary>
        /// <param name="n">方程个数，也是未知数个数</param>
        /// <param name="eps">控制精度要求</param>
        /// <param name="x">存放一组初值X0,X1,…,Xn-1，返回一组实根</param>
        /// <param name="js">允许最大迭代次数</param>
        /// <param name="func">指向计算目标函数(Σf²)值与偏导数值的函数.
        /// <para>参数1：给定的初始X点</para>
        /// <para>参数2：输出X处的目标函数偏导数值</para>
        /// <para>参数3：方程数</para>
        /// <para>返回：在给定的X处的目标函数值</para>
        /// </param>
        /// <returns>函数返回实际迭代次数。
        /// 当返回值小于0时，说明在迭代过程中遇到了D=O,此时可改变初值再试；
        /// 当返回值等千允许的最大迭代次数js, 说明迭代了js次还未满足精度要求， 
        /// 此时可以改变初值或放松精度要求或增大js值后再进行尝试</returns>
        public static int SNSE(int n, double eps, double[] x, int js, Func<double[], double[], int, double> func)
        {
            int l = js, j;
            double f, d, s;
            double[] y = new double[n];
            f = func(x, y, n);
            while (f >= eps)
            {
                l--;
                if (l == 0) return js;
                d = 0;
                for (j = 0; j <= n - 1; j++) d += y[j] * y[j];
                if (d == 0) return -1;
                s = f / d;
                for (j = 0; j <= n - 1; j++) x[j] -= s * y[j];
                f = func(x, y, n);
            }
            return js - 1;
        }

        /// <summary>
        /// 用拟牛顿法求非线性方程组的一组实数解。
        /// </summary>
        /// <param name="n">方程组中方程个数，也是未知数个数</param>
        /// <param name="eps">控制精度要求</param>
        /// <param name="t">控制 h 大小的变量， O<t<l</param>
        /// <param name="h">增量初值，在本函数中要被破坏</param>
        /// <param name="x">存放初值，返回方程组的一组实数解</param>
        /// <param name="k">允许的最大迭代次数</param>
        /// <param name="func">指向计算方程组左端函数值的函数
        /// <para>参数1：点x</para>
        /// <para>参数2：输出y</para>
        /// <para>参数3：方程数</para>
        /// </param>
        /// <returns>函数返回实际迭代次数。若返回值等于0, 则说明迭代了 K 次还未满足精度要求；
        /// 若返回值等于-1, 则说明线性代数方程组 AZ = B 奇异；
        /// 若返回值等于-2，则说明ΣZj=1</returns>
        public static int NETN(int n, double eps, double t, double h, double[] x, int k, Action<double[], double[], int> func)
        {
            if (n != x.Length)
            {
                throw new Exception("x长度与指定的n不一致");
            }
            double[] y = new double[n];
            double[] a = new double[n * n];
            double[] b = new double[n];
            int l = k, i, j;
            double am = 1 + eps, z, beta, d;
            while (am >= eps)
            {
                func(x, b, n);
                am = 0;
                for (i = 0; i <= n - 1; i++)
                {
                    z = Math.Abs(b[i]);
                    if (z > am) am = z;
                }
                if (am >= eps)
                {
                    l--;
                    if (l == 0)
                    {
                        return 0;//fail
                    }
                    for (j = 0; j <= n - 1; j++)
                    {
                        z = x[j];
                        x[j] += h;
                        func(x, y, n);
                        for (i = 0; i <= n - 1; i++) a[i * n + j] = y[i];
                        x[j] = z;
                    }
                    try
                    {
                        LinearEquations.GAUS(Utility.C.Convert(a, n, n), b);
                    }
                    catch
                    {
                        return -2;//fail
                    }
                    beta = 1;
                    for (i = 0; i <= n - 1; i++) beta -= b[i];
                    if (beta == 0)
                    {
                        return -2;//fail
                    }
                    d = h / beta;
                    for (i = 0; i <= n - 1; i++) x[i] = x[i] - d * b[i];
                    h = t * h;
                }
            }
            return k - l;
        }

        /// <summary>
        /// 利用广义逆法求解无约束条件下的优化问题
        /// <para>f(X0,X1,...,Xn-1)=0,i=0,1,...,m-1,m≥n</para>
        /// 当 m = n 时，即为求解非线性方程组
        /// </summary>
        /// <param name="m">非线性方程组中方程个数</param>
        /// <param name="n">非线性方程组中未知数个数</param>
        /// <param name="eps1">控制最小二乘解的精度要求</param>
        /// <param name="eps2">用于奇异值分解中的控制精度要求</param>
        /// <param name="x">存放非线性方程组解的初始近似值X(O),要求各分扯不全为0。返回最小二乘解，当 m = n 时，即为非线性方程组的一组解</param>
        /// <param name="ka">ka=max(m,n)+1</param>
        /// <param name="f">指向计算非线性方程组中各方程左端函数值的函数</param>
        /// <param name="s">指向计算雅可比矩阵的函数名</param>
        /// <returns>函数返回一个标志值 。 
        /// 若返回值小于 0, 则说明在 奇异值分解 中迭代超过了 60 次还未满足精度要求；
        /// 若返回值等千 0, 则说明本函数迭代次还未满足精度要求；
        /// 若返回值大于 o, 则表示正常返回</returns>
        public static int NGIN(int m, int n, double eps1, double eps2, double[] x, int ka, NGIN_Func f, NGIN_Func s)
        {
            int i, j, k, l = 60, kk, jt;
            double[] y = new double[10], b = new double[10], p = new double[m*n], d = new double[m], pp = new double[n*m], dx = new double[n], u=new double[m*m], v = new double[n*n], w = new double[ka];
            double alpha = 1.0, z=0, h2, y1, y2, y3, y0, h1;
            double[,] p_t, pp_t, u_t, v_t;
            while (l>0)
            {
                f(m, n, x, ref d);
                s(m, n, x, ref p);
                p_t = Utility.C.Convert(p, m, n);
                try
                {
                    dx = LinearEquations.GMIV(p_t, d, eps2, out pp_t, out u_t, out v_t);
                }
                catch
                {
                    return -1;
                }
                p = Utility.C.Convert(p_t);
                pp = Utility.C.Convert(pp_t);
                u = Utility.C.Convert(u_t);
                v = Utility.C.Convert(v_t);
                j = 0;
                jt = 1;
                h2 = 0;
                while (jt==1)
                {
                    jt = 0;
                    if (j <= 2) z = alpha + 0.01 * j;
                    else z = h2;
                    for (i = 0; i <= n - 1; i++) w[i] = x[i] - z * dx[i];
                    f(m, n, w, ref d);
                    y1 = 0;
                    for (i = 0; i <= m - 1; i++) y1 = y1 + d[i] * d[i];
                    for (i = 0; i <= n - 1; i++) w[i] = x[i] - (z + 0.00001) * dx[i];
                    f(m, n, w, ref d);
                    y2 = 0;
                    for (i = 0; i <= m - 1; i++) y2 = y2 + d[i] * d[i];
                    y0 = (y2 - y1) / 0.00001;
                    if (Math.Abs(y0)>1e-10)
                    {
                        h1 = y0;
                        h2 = z;
                        if (j == 0)
                        {
                            y[0] = h1;
                            b[0] = h2;
                        }
                        else
                        {
                            y[j] = h1;
                            kk = 0;
                            k = 0;
                            while ((kk==0)&&(k<=j-1))
                            {
                                y3 = h2 - b[k];
                                if (Math.Abs(y3) == 0) kk = 1;
                                else h2 = (h1 - y[k]) / y3;
                                k++;
                            }
                            b[j] = h2;
                            if (kk != 0) b[j] = 1e35;
                            h2 = 0;
                            for (k = j - 1; k >= 0; k--) h2 = -y[k] / (b[k + 1] + h2);
                            h2 += b[0];
                        }
                        j++;
                        if (j <= 7) jt = 1;
                        else z = h2;
                    }
                }
                alpha = z;
                y1 = 0;
                y2 = 0;
                for (i = 0; i <= n - 1; i++)
                {
                    dx[i] = -alpha * dx[i];
                    x[i] = x[i] + dx[i];
                    y1 = y1 + Math.Abs(dx[i]);
                    y2 = y2 + Math.Abs(x[i]);
                }
                if (y1 < eps1 * y2)
                {
                    return 1;
                }
                l--;
            }
            return 0;
        }
        /// <summary>
        /// 指向计算非线性方程组中各方程左端函数值的函数，或计算雅可比矩阵的函数名
        /// </summary>
        /// <param name="m">非线性方程组中方程个数</param>
        /// <param name="n">非线性方程组中未知数个数</param>
        /// <param name="x">存放非线性方程组解的初始近似值,要求各分扯不全为 0 </param>
        /// <param name="dp">左端函数值，或雅可比矩阵</param>
        public delegate void NGIN_Func(int m, int n, double[] x, ref double[] dp);
    }
}
