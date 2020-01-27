using System;
using System.Collections.Generic;
using System.Text;

namespace Heroius.XuAlgrithms
{
    /// <summary>
    /// 线性方程组的求解
    /// </summary>
    public static partial class LinearEquations
    {
        /// <summary>
        /// 用全选主元高斯（Gauss）消去法求解n阶线性代数方程组AX=B.
        /// </summary>
        /// <param name="a">存放方程组的系数矩阵，返回时将被破坏</param>
        /// <param name="b">存放方程组右端的常数向量，返回方程组的解向量</param>
        public static void GAUS(double[,] a, double[] b)
        {
            int n = a.GetLength(0);
            if (a.GetLength(1) != n)
            {
                throw new Exception("非系数矩阵");
            }
            if (b.Length != n)
            {
                throw new Exception("非常数向量");
            }
            int[] js;
            int l, k, i, j, s = 0;
            double d, t;
            js = new int[n];
            l = 1;
            for (k = 0; k <= n - 2; k++)
            {
                d = 0.0;
                for (i = k; i <= n - 1; i++)
                    for (j = k; j <= n - 1; j++)
                    {
                        t = Math.Abs(a[i, j]);
                        if (t > d) { d = t; js[k] = j; s = i; }
                    }
                if (d + 1.0 == 1.0) l = 0;
                else
                {
                    if (js[k] != k)
                        for (i = 0; i <= n - 1; i++)
                        {
                            t = a[i, k];
                            a[i, k] = a[i, js[k]];
                            a[i, js[k]] = t;
                        }
                    if (s != k)
                    {
                        for (j = k; j <= n - 1; j++)
                        {
                            t = a[k, j];
                            a[k, j] = a[s, j];
                            a[s, j] = t;
                        }
                        t = b[k]; b[k] = b[s]; b[s] = t;
                    }
                }
                if (l == 0)
                {
                    throw new Exception("系数矩阵奇异！无解.");
                }
                d = a[k, k];
                for (j = k + 1; j <= n - 1; j++)
                    a[k, j] = a[k, j] / d;
                b[k] = b[k] / d;
                for (i = k + 1; i <= n - 1; i++)
                {
                    for (j = k + 1; j <= n - 1; j++)
                        a[i, j] = a[i, j] - a[i, k] * a[k, j];
                    b[i] = b[i] - a[i, k] * b[k];
                }
            }
            d = a[n - 1, n - 1];
            if (Math.Abs(d) + 1.0 == 1.0)
            {
                throw new Exception("系数矩阵奇异！无解.");
            }
            b[n - 1] = b[n - 1] / d;
            for (i = n - 2; i >= 0; i--)
            {
                t = 0.0;
                for (j = i + 1; j <= n - 1; j++)
                    t = t + a[i, j] * b[j];
                b[i] = b[i] - t;
            }
            js[n - 1] = n - 1;
            for (k = n - 1; k >= 0; k--)
                if (js[k] != k)
                {
                    t = b[k]; b[k] = b[js[k]]; b[js[k]] = t;
                }
        }

        /// <summary>
        /// 用全选主元高斯-约当（Gauss-Jordan）消去法同时求解系数矩阵相同而右端具有m组常数向量的n阶线性代数方程组AX=B.
        /// </summary>
        /// <param name="a">存放诚成祖的系数矩阵，返回时将被破坏</param>
        /// <param name="b">存放方程组右端的m组常数向量，返回方程组的m组解向量</param>
        public static void GJDN(double[,] a, double[,] b)
        {
            int n = a.GetLength(0);
            if (a.GetLength(1) != n)
            {
                throw new Exception("非系数矩阵");
            }
            if (b.GetLength(0) != n)
            {
                throw new Exception("非常数向量");
            }
            int m = b.GetLength(1);
            int l, k, i, j, ss = 0;
            double d, t;
            int[] js = new int[n];
            l = 1;
            for (k = 0; k <= n - 1; k++)
            {
                d = 0.0;
                for (i = k; i <= n - 1; i++)
                    for (j = k; j <= n - 1; j++)
                    {
                        t = Math.Abs(a[i, j]);
                        if (t > d) { d = t; js[k] = j; ss = i; }
                    }
                if (d + 1.0 == 1.0) l = 0;
                else
                {
                    if (js[k] != k)
                        for (i = 0; i <= n - 1; i++)
                        {
                            t = a[i, k];
                            a[i, k] = a[i, js[k]];
                            a[i, js[k]] = t;
                        }
                    if (ss != k)
                    {
                        for (j = k; j <= n - 1; j++)
                        {
                            t = a[k, j];
                            a[k, j] = a[ss, j];
                            a[ss, j] = t;
                        }
                        for (j = 0; j <= m - 1; j++)
                        {
                            t = b[k, j];
                            b[k, j] = b[ss, j];
                            b[ss, j] = t;
                        }
                    }
                }
                if (l == 0)
                {
                    throw new Exception("系数矩阵奇异！无解. ");
                }
                d = a[k, k];
                for (j = k + 1; j <= n - 1; j++)
                    a[k, j] = a[k, j] / d;
                for (j = 0; j <= m - 1; j++)
                    b[k, j] = b[k, j] / d;
                for (j = k + 1; j <= n - 1; j++)
                    for (i = 0; i <= n - 1; i++)
                    {
                        if (i != k)
                            a[i, j] = a[i, j] - a[i, k] * a[k, j];
                    }
                for (j = 0; j <= m - 1; j++)
                    for (i = 0; i <= n - 1; i++)
                    {
                        if (i != k)
                            b[i, j] = b[i, j] - a[i, k] * b[k, j];
                    }
            }
            for (k = n - 1; k >= 0; k--)
                if (js[k] != k)
                    for (j = 0; j <= m - 1; j++)
                    {
                        t = b[k, j]; b[k, j] = b[js[k], j]; b[js[k], j] = t;
                    }
        }

        /// <summary>
        /// 用全选主元高斯(Gauss)消去法求解n 阶复系数线性代数方程组AX=B。
        /// </summary>
        /// <param name="ar">存放方程组复系数矩阵的实部，返回时将被破坏</param>
        /// <param name="ai">存放方程组复系数矩阵的虚部，返回时将被破坏</param>
        /// <param name="br">存放方程组右端复常数向量的实部，返回方程组解向址的实部</param>
        /// <param name="bi">存放方程组右端复常数向量的虚部，返回方程组解向址的虚部</param>
        public static void CGAS(double[,] ar, double[,] ai, double[] br, double[] bi)
        {
            int n = ar.GetLength(0);
            if (ar.GetLength(1) != n)
            {
                throw new Exception("系数实部非n阶矩阵");
            }
            if (ai.GetLength(0) != n || ai.GetLength(1) != n)
            {
                throw new Exception("系数虚部非n阶矩阵");
            }
            if (br.Length != n || bi.Length != n)
            {
                throw new Exception("常数实部或虚部n阶向量");
            }
            int k, i, j, ss = 0;
            double p, q, s, d;
            var js = new int[n];
            for (k = 0; k <= n - 2; k++)
            {
                d = 0.0;
                for (i = k; i <= n - 1; i++)
                    for (j = k; j <= n - 1; j++)
                    {
                        p = ar[i, j] * ar[i, j] + ai[i, j] * ai[i, j];
                        if (p > d) { d = p; js[k] = j; ss = i; }
                    }
                if (d + 1.0 == 1.0)
                {
                    throw new Exception("系数矩阵奇异！无解. ");
                }
                if (ss != k)
                {
                    for (j = k; j <= n - 1; j++)
                    {
                        p = ar[k, j]; ar[k, j] = ar[ss, j]; ar[ss, j] = p;
                        p = ai[k, j]; ai[k, j] = ai[ss, j]; ai[ss, j] = p;
                    }
                    p = br[k]; br[k] = br[ss]; br[ss] = p;
                    p = bi[k]; bi[k] = bi[ss]; bi[ss] = p;
                }
                if (js[k] != k)
                    for (i = 0; i <= n - 1; i++)
                    {
                        p = ar[i, k]; ar[i, k] = ar[i, js[k]]; ar[i, js[k]] = p;
                        p = ai[i, k]; ai[i, k] = ai[i, js[k]]; ai[i, js[k]] = p;
                    }
                for (j = k + 1; j <= n - 1; j++)
                {
                    p = ar[k, j] * ar[k, k]; q = -ai[k, j] * ai[k, k];
                    s = (ar[k, k] - ai[k, k]) * (ar[k, j] + ai[k, j]);
                    ar[k, j] = (p - q) / d; ai[k, j] = (s - p - q) / d;
                }
                p = br[k] * ar[k, k]; q = -bi[k] * ai[k, k];
                s = (ar[k, k] - ai[k, k]) * (br[k] + bi[k]);
                br[k] = (p - q) / d; bi[k] = (s - p - q) / d;
                for (i = k + 1; i <= n - 1; i++)
                {
                    for (j = k + 1; j <= n - 1; j++)
                    {
                        p = ar[i, k] * ar[k, j]; q = ai[i, k] * ai[k, j];
                        s = (ar[i, k] + ai[i, k]) * (ar[k, j] + ai[k, j]);
                        ar[i, j] = ar[i, j] - p + q;
                        ai[i, j] = ai[i, j] - s + p + q;
                    }
                    p = ar[i, k] * br[k]; q = ai[i, k] * bi[k];
                    s = (ar[i, k] + ai[i, k]) * (br[k] + bi[k]);
                    br[i] = br[i] - p + q; bi[i] = bi[i] - s + p + q;
                }
            }
            d = ar[n - 1, n - 1] * ar[n - 1, n - 1] + ai[n - 1, n - 1] * ai[n - 1, n - 1];
            if (d + 1.0 == 1.0)
            {
                throw new Exception("系数矩阵奇异！无解. ");
            }
            p = ar[n - 1, n - 1] * br[n - 1]; q = -ai[n - 1, n - 1] * bi[n - 1];
            s = (ar[n - 1, n - 1] - ai[n - 1, n - 1]) * (br[n - 1] + bi[n - 1]);
            br[n - 1] = (p - q) / d; bi[n - 1] = (s - p - q) / d;
            for (i = n - 2; i >= 0; i--)
                for (j = i + 1; j <= n - 1; j++)
                {
                    p = ar[i, j] * br[j]; q = ai[i, j] * bi[j];
                    s = (ar[i, j] + ai[i, j]) * (br[j] + bi[j]);
                    br[i] = br[i] - p + q;
                    bi[i] = bi[i] - s + p + q;
                }
            js[n - 1] = n - 1;
            for (k = n - 1; k >= 0; k--)
                if (js[k] != k)
                {
                    p = br[k]; br[k] = br[js[k]]; br[js[k]] = p;
                    p = bi[k]; bi[k] = bi[js[k]]; bi[js[k]] = p;
                }
        }

        /// <summary>
        /// 用全选主元高斯－约当(Gauss-Jordan) 消去法同时求解系数矩阵相同而右端具有m 组复常数向量的n 阶复系数线性代数方程组AX=B 。
        /// </summary>
        /// <param name="ar">存放方程组复系数矩阵的实部，返回时将被破坏</param>
        /// <param name="ai">存放方程组复系数矩阵的虚部，返回时将被破坏</param>
        /// <param name="br">存放方程组右端m 组复常数向最的实部, 返回方程组复解向最的实部</param>
        /// <param name="bi">存放方程组右端m 组复常数向最的虚部, 返回方程组复解向最的虚部</param>
        public static void CJDN(double[,] ar, double[,] ai, double[,] br, double[,] bi)
        {
            int n = ar.GetLength(0);
            int m = br.GetLength(1);
            if (ar.GetLength(1) != n)
            {
                throw new Exception("系数实部非n阶矩阵");
            }
            if (ai.GetLength(0) != n || ai.GetLength(1) != n)
            {
                throw new Exception("系数虚部非n阶矩阵");
            }
            if (br.GetLength(0) != n || bi.GetLength(0) != n)
            {
                throw new Exception("常数实部或虚部行数不正确");
            }
            if (bi.GetLength(1) != m)
            {
                throw new Exception("常数虚部行数不正确");
            }

            int k, i, j, ss = 0;
            double p, q, s, d;
            var js = new int[n];
            for (k = 0; k <= n - 1; k++)
            {
                d = 0.0;
                for (i = k; i <= n - 1; i++)
                    for (j = k; j <= n - 1; j++)
                    {
                        p = ar[i, j] * ar[i, j] + ai[i, j] * ai[i, j];
                        if (p > d) { d = p; js[k] = j; ss = i; }
                    }
                if (d + 1.0 == 1.0)
                {
                    throw new Exception("系数矩阵奇异！无解. ");
                }
                if (ss != k)
                {
                    for (j = k; j <= n - 1; j++)
                    {
                        p = ar[k, j]; ar[k, j] = ar[ss, j]; ar[ss, j] = p;
                        p = ai[k, j]; ai[k, j] = ai[ss, j]; ai[ss, j] = p;
                    }
                    for (j = 0; j <= m - 1; j++)
                    {
                        p = br[k, j]; br[k, j] = br[ss, j]; br[ss, j] = p;
                        p = bi[k, j]; bi[k, j] = bi[ss, j]; bi[ss, j] = p;
                    }
                }
                if (js[k] != k)
                    for (i = 0; i <= n - 1; i++)
                    {
                        p = ar[i, k]; ar[i, k] = ar[i, js[k]]; ar[i, js[k]] = p;
                        p = ai[i, k]; ai[i, k] = ai[i, js[k]]; ai[i, js[k]] = p;
                    }
                for (j = k + 1; j <= n - 1; j++)
                {
                    p = ar[k, j] * ar[k, k];
                    q = -ai[k, j] * ai[k, k];
                    s = (ar[k, k] - ai[k, k]) * (ar[k, j] + ai[k, j]);
                    ar[k, j] = (p - q) / d;
                    ai[k, j] = (s - p - q) / d;
                }
                for (j = 0; j <= m - 1; j++)
                {
                    p = br[k, j] * ar[k, k];
                    q = -bi[k, j] * ai[k, k];
                    s = (ar[k, k] - ai[k, k]) * (br[k, j] + bi[k, j]);
                    br[k, j] = (p - q) / d;
                    bi[k, j] = (s - p - q) / d;
                }
                for (i = 0; i <= n - 1; i++)
                    if (i != k)
                    {
                        for (j = k + 1; j <= n - 1; j++)
                        {
                            p = ar[i, k] * ar[k, j];
                            q = ai[i, k] * ai[k, j];
                            s = (ar[i, k] + ai[i, k]) * (ar[k, j] + ai[k, j]);
                            ar[i, j] = ar[i, j] - p + q;
                            ai[i, j] = ai[i, j] - s + p + q;
                        }
                        for (j = 0; j <= m - 1; j++)
                        {
                            p = ar[i, k] * br[k, j];
                            q = ai[i, k] * bi[k, j];
                            s = (ar[i, k] + ai[i, k]) * (br[k, j] + bi[k, j]);
                            br[i, j] = br[i, j] - p + q;
                            bi[i, j] = bi[i, j] - s + p + q;
                        }
                    }
            }
            for (k = n - 1; k >= 0; k--)
                if (js[k] != k)
                    for (j = 0; j <= m - 1; j++)
                    {
                        p = br[k, j]; br[k, j] = br[js[k], j]; br[js[k], j] = p;
                        p = bi[k, j]; bi[k, j] = bi[js[k], j]; bi[js[k], j] = p;
                    }
        }

        /// <summary>
        /// 用追赶法求解n 阶三对角线方程组AX=D 。
        /// </summary>
        /// <param name="b">以行为主存放三对角矩阵中3 条对角线上的元素</param>
        /// <param name="m">三对角矩阵3 条对角线上的元素个数, 其值应为m =3 n- 2</param>
        /// <param name="n">方程组的阶数</param>
        /// <param name="d">存放方程组右端的常数向量，返回方程组的解向拭</param>
        public static void TRDE(double[] b, int n, int m, ref double[] d)
        {
            //int m = b.Length;
            //int n = d.Length;
            if (m != 3 * n - 2)
            {
                throw new Exception("三对角矩阵3 条对角线上的元素个数应为m = 3 n - 2");
            }
            int k, j;
            double s;
            for (k = 0; k <= n - 2; k++)
            {
                j = 3 * k; s = b[j];
                if (Math.Abs(s) + 1.0 == 1.0)
                {
                    throw new Exception("程序工作失败！无解.");
                }
                b[j + 1] = b[j + 1] / s;
                d[k] = d[k] / s;
                b[j + 3] = b[j + 3] - b[j + 2] * b[j + 1];
                d[k + 1] = d[k + 1] - b[j + 2] * d[k];
            }
            s = b[3 * n - 3];
            if (Math.Abs(s) + 1.0 == 1.0)
            {
                throw new Exception("程序工作失败！无解.");
            }
            d[n - 1] = d[n - 1] / s;
            for (k = n - 2; k >= 0; k--)
                d[k] = d[k] - b[3 * k + 1] * d[k + 1];
        }

        /// <summary>
        /// 用列选主元高斯(Gauss) 消去法求解右端具有m 组常数向量的n 阶一般带型方程组AX=D 。其中A 为n 阶带型矩阵。
        /// </summary>
        /// <param name="b">存放带型矩阵A 中带区内的元素，返回时将被破坏</param>
        /// <param name="d">存放方程组右端的m 组常数向散，返回方程组的m 组解向量</param>
        public static void BAND(double[,] b, double[,] d)
        {
            int n = b.GetLength(0);
            int il = b.GetLength(1);
            int m = d.GetLength(1);
            if (d.GetLength(0) != n)
            {
                throw new Exception("b 与 d 不兼容");
            }
            int ls, k, i, j, ss = 0;
            double p, t;
            ls = (il - 1) / 2;
            for (k = 0; k <= n - 2; k++)
            {
                p = 0.0;
                for (i = k; i <= ls; i++)
                {
                    t = Math.Abs(b[i, 0]);
                    if (t > p) { p = t; ss = i; }
                }
                if (p + 1.0 == 1.0)
                {
                    throw new Exception("程序工作失败！无解.");
                }
                for (j = 0; j <= m - 1; j++)
                {
                    t = d[k, j];
                    d[k, j] = d[ss, j];
                    d[ss, j] = t;
                }
                for (j = 0; j <= il - 1; j++)
                {
                    t = b[k, j];
                    b[k, j] = b[ss, j];
                    b[ss, j] = t;
                }
                for (j = 0; j <= m - 1; j++)
                    d[k, j] = d[k, j] / b[k, 0];
                for (j = 1; j <= il - 1; j++)
                    b[k, j] = b[k, j] / b[k, 0];
                for (i = k + 1; i <= ls; i++)
                {
                    t = b[i, 0];
                    for (j = 0; j <= m - 1; j++)
                        d[i, j] = d[i, j] - t * d[k, j];
                    for (j = 1; j <= il - 1; j++)
                        b[i, j - 1] = b[i, j] - t * b[k, j];
                    b[i, il - 1] = 0.0;
                }
                if (ls != (n - 1)) ls = ls + 1;
            }
            p = b[n - 1, 0];
            if (Math.Abs(p) + 1.0 == 1.0)
            {
                throw new Exception("程序工作失败！无解.");
            }
            for (j = 0; j <= m - 1; j++)
                d[n - 1, j] = d[n - 1, j] / p;
            ls = 1;
            for (i = n - 2; i >= 0; i--)
            {
                for (k = 0; k <= m - 1; k++)
                {
                    for (j = 1; j <= ls; j++)
                    {
                        d[i, k] = d[i, k] - b[i, j] * d[i + j, k];
                    }
                }
                if (ls != (il - 1)) ls = ls + 1;
            }
        }

        /// <summary>
        /// 用分解法求解系数矩阵为对称、且右端具有m 组常数向量的线性代数方程组 AX = C 。
        /// </summary>
        /// <param name="a">存放方程组的系数矩阵（应为对称矩阵），返回时将被破坏</param>
        /// <param name="c">存放方程组右端m 组常数向量，返回方程组的m 组解向扯</param>
        public static void LDLE(double[,] a, double[,] c)
        {
            int n = a.GetLength(0);
            if (a.GetLength(1) != n)
            {
                throw new Exception("a 非方阵");
            }
            if (c.GetLength(0) != n)
            {
                throw new Exception("c 的行数不匹配");
            }
            int m = c.GetLength(1);
            int i, j, k, k1, k2, k3;
            double p;
            if (Math.Abs(a[0, 0]) + 1.0 == 1.0)
            {
                throw new Exception("系数矩阵奇异！无解.");
            }
            for (i = 1; i <= n - 1; i++)
                a[i, 0] = a[i, 0] / a[0, 0];
            for (i = 1; i <= n - 2; i++)
            {
                for (j = 1; j <= i; j++)
                {
                    a[i, i] = a[i, i] - a[i, j - 1] * a[i, j - 1] * a[j - 1, j - 1];
                }
                p = a[i, i];
                if (Math.Abs(p) + 1.0 == 1.0)
                {
                    throw new Exception("系数矩阵奇异！无解.");
                }
                for (k = i + 1; k <= n - 1; k++)
                {
                    for (j = 1; j <= i; j++)
                    {
                        a[k, i] = a[k, i] - a[k, j - 1] * a[i, j - 1] * a[j - 1, j - 1];
                    }
                    a[k, i] = a[k, i] / p;
                }
            }
            for (j = 1; j <= n - 1; j++)
            {
                a[n - 1, n - 1] = a[n - 1, n - 1] - a[n - 1, j - 1] * a[n - 1, j - 1] * a[j - 1, j - 1];
            }
            p = a[n - 1, n - 1];
            if (Math.Abs(p) + 1.0 == 1.0)
            {
                throw new Exception("系数矩阵奇异！无解.");
            }
            for (j = 0; j <= m - 1; j++)
                for (i = 1; i <= n - 1; i++)
                {
                    for (k = 1; k <= i; k++)
                    {
                        c[i, j] = c[i, j] - a[i, k - 1] * c[k - 1, j];
                    }
                }
            for (i = 1; i <= n - 1; i++)
            {
                for (j = i; j <= n - 1; j++)
                {
                    a[i - 1, j] = a[i - 1, i - 1] * a[j, i - 1];
                }
            }
            for (j = 0; j <= m - 1; j++)
            {
                c[n - 1, j] = c[n - 1, j] / p;
                for (k = 1; k <= n - 1; k++)
                {
                    k1 = n - k; k3 = k1 - 1;
                    for (k2 = k1; k2 <= n - 1; k2++)
                    {
                        c[k3, j] = c[k3, j] - a[k3, k2] * c[k2, j];
                    }
                    c[k3, j] = c[k3, j] / a[k3, k3];
                }
            }
        }

        /// <summary>
        /// 用乔里斯基(Cholesky) 分解法（即平方根法）求解系数矩阵为对称正定、且右端具有m 组常数向量的n 阶线性代数方程组AX=D 。其中A 为n 阶对称正定矩阵。
        /// </summary>
        /// <param name="a">a[n,n]: 存放对称正定的系数矩阵，返回时其上三角部分存放分解后的U 矩阵</param>
        /// <param name="n">方程组的阶数</param>
        /// <param name="m">方程组右端常数向量的组数</param>
        /// <param name="d">d[n,m]: 存放方程组右端m 组常数向量，返回方程组的m 组解向量</param>
        /// <returns>函数返回标志值。若返回的标志值小千0, 则表示程序工作失败； 若返回的标志值大于0, 则表示正常返回</returns>
        public static int CHLK(ref double[,] a, int n, int m, ref double[,] d)
        {
            //int n = a.GetLength(0);
            //if (a.GetLength(1) != n)
            //{
            //    throw new Exception("a 非方阵");
            //}
            //int m = d.GetLength(1);
            //if (d.GetLength(0) != n)
            //{
            //    throw new Exception("d 不匹配");
            //}
            int i, j, k;
            if ((a[0, 0] + 1.0 == 1.0) || (a[0, 0] < 0.0))
            {
                //throw new Exception("程序工作失败！无解.");
                return -2;
            }
            a[0, 0] = Math.Sqrt(a[0, 0]);
            for (j = 1; j <= n - 1; j++) a[0, j] = a[0, j] / a[0, 0];
            for (i = 1; i <= n - 1; i++)
            {
                for (j = 1; j <= i; j++)
                {
                    a[i, i] = a[i, i] - a[j - 1, i] * a[j - 1, i];
                }
                if ((a[i, i] + 1.0 == 1.0) || (a[i, i] < 0.0))
                {
                    //throw new Exception("程序工作失败！无解.");
                    return -2;
                }
                a[i, i] = Math.Sqrt(a[i, i]);
                if (i != (n - 1))
                {
                    for (j = i + 1; j <= n - 1; j++)
                    {
                        for (k = 1; k <= i; k++)
                            a[i, j] = a[i, j] - a[k - 1, i] * a[k - 1, j];
                        a[i, j] = a[i, j] / a[i, i];
                    }
                }
            }
            for (j = 0; j <= m - 1; j++)
            {
                d[0, j] = d[0, j] / a[0, 0];
                for (i = 1; i <= n - 1; i++)
                {
                    for (k = 1; k <= i; k++)
                        d[i, j] = d[i, j] - a[k - 1, i] * d[k - 1, j];
                    d[i, j] = d[i, j] / a[i, i];
                }
            }
            for (j = 0; j <= m - 1; j++)
            {
                d[n - 1, j] = d[n - 1, j] / a[n - 1, n - 1];
                for (k = n - 1; k >= 1; k--)
                {
                    for (i = k; i <= n - 1; i++)
                    {
                        d[k - 1, j] = d[k - 1, j] - a[k - 1, i] * d[i, j];
                    }
                    d[k - 1, j] = d[k - 1, j] / a[k - 1, k - 1];
                }
            }
            return 2;
        }

        /// <summary>
        /// 用列文逊(Levinson) 递推算法求解n 阶对称托伯利兹(Toeplitz) 型方程组。返回方程组的解向量。
        /// </summary>
        /// <param name="t">存放n 阶T 型矩阵中的元素t0, t1 … tn-1</param>
        /// <param name="b">存放方程组右端的常数向量</param>
        /// <returns>返回方程组的解向量</returns>
        public static double[] TLVS(double[] t, double[] b)
        {
            int n = t.Length;
            if (b.Length != n)
            {
                throw new Exception("b 不匹配");
            }
            int i, j, k;
            double a, beta, q, c, h;
            double[] s = new double[n];
            double[] y = new double[n];
            double[] x = new double[n];
            a = t[0];
            if (Math.Abs(a) + 1.0 == 1.0)
            {
                throw new Exception("程序工作失败！");
            }
            y[0] = 1.0; x[0] = b[0] / a;
            for (k = 1; k <= n - 1; k++)
            {
                beta = 0.0; q = 0.0;
                for (j = 0; j <= k - 1; j++)
                {
                    beta = beta + y[j] * t[j + 1];
                    q = q + x[j] * t[k - j];
                }
                if (Math.Abs(a) + 1.0 == 1.0)
                {
                    throw new Exception("程序工作失败！");
                }
                c = -beta / a; s[0] = c * y[k - 1]; y[k] = y[k - 1];
                if (k != 1)
                    for (i = 1; i <= k - 1; i++)
                        s[i] = y[i - 1] + c * y[k - i - 1];
                a = a + c * beta;
                if (Math.Abs(a) + 1.0 == 1.0)
                {
                    throw new Exception("程序工作失败！");
                }
                h = (b[k] - q) / a;
                for (i = 0; i <= k - 1; i++)
                { x[i] = x[i] + h * s[i]; y[i] = s[i]; }
                x[k] = h * y[k];
            }
            return x;
        }

        /// <summary>
        /// 用高斯－赛德尔(Gauss-Seidel) 迭代法求解系数矩阵具有主对角线占绝对优势的线性代数方程组AX=B 。返回方程组的解向量。
        /// </summary>
        /// <param name="a">存放方程组的系数矩阵</param>
        /// <param name="b">存放方程组右端的常数向量</param>
        /// <param name="eps">给定的精度要求</param>
        /// <returns>返回方程组的解向量</returns>
        public static double[] GSDL(double[,] a, double[] b, double eps)
        {
            int n = a.GetLength(0);
            if (a.GetLength(1) != n)
            {
                throw new Exception("a 非方阵");
            }
            if (b.Length != n)
            {
                throw new Exception("b 不匹配");
            }
            double[] x = new double[n];
            int i, j;
            double p, t, s, q;
            for (i = 0; i <= n - 1; i++)
            {
                p = 0.0; x[i] = 0.0;
                for (j = 0; j <= n - 1; j++)
                    if (i != j) p = p + Math.Abs(a[i, j]);
                if (p >= Math.Abs(a[i, i]))
                {
                    throw new Exception("程序工作失败！系数矩阵不具有主对角线占绝对优势");
                }
            }
            p = eps + 1.0;
            while (p >= eps)
            {
                p = 0.0;
                for (i = 0; i <= n - 1; i++)
                {
                    t = x[i]; s = 0.0;
                    for (j = 0; j <= n - 1; j++)
                        if (j != i) s = s + a[i, j] * x[j];
                    x[i] = (b[i] - s) / a[i, i];
                    q = Math.Abs(x[i] - t) / (1.0 + Math.Abs(x[i]));
                    if (q > p) p = q;
                }
            }
            return x;
        }

        /// <summary>
        /// 用共轭梯度法求解n 阶对称正定方程组AX = B 。返回方程组的解向量。
        /// </summary>
        /// <param name="a">存放对称正定矩阵A</param>
        /// <param name="b">存放方程组右端的常数向量</param>
        /// <param name="eps">控制精度要求</param>
        /// <returns>返回方程组的解向量</returns>
        public static double[] GRAD(double[,] a, double[] b, double eps)
        {
            int n = a.GetLength(0);
            if (a.GetLength(1) != n)
            {
                throw new Exception("a 非方阵");
            }
            if (b.Length != n)
            {
                throw new Exception("b 不匹配");
            }
            int i, k, j;
            double alpha, beta, d, e;
            double[] p = new double[n];
            double[] r = new double[n];
            double[] s = new double[n];
            double[] q = new double[n];
            double[] x = new double[n];
            for (i = 0; i <= n - 1; i++)
            { x[i] = 0.0; p[i] = b[i]; r[i] = b[i]; }
            i = 0;
            while (i <= n - 1)
            {
                for (k = 0; k < n; k++)
                {
                    s[k] = 0.0;
                    for (j = 0; j < n; j++)
                        s[k] = s[k] + a[k, j] * p[j];
                }
                d = 0.0; e = 0.0;
                for (k = 0; k <= n - 1; k++)
                {
                    d = d + p[k] * b[k];
                    e = e + p[k] * s[k];
                }
                alpha = d / e;
                for (k = 0; k <= n - 1; k++)
                    x[k] = x[k] + alpha * p[k];
                for (k = 0; k < n; k++)
                {
                    q[k] = 0.0;
                    for (j = 0; j < n; j++)
                        q[k] = q[k] + a[k, j] * x[j];
                }
                d = 0.0;
                for (k = 0; k <= n - 1; k++)
                {
                    r[k] = b[k] - q[k];
                    d = d + r[k] * s[k];
                }
                beta = d / e; d = 0.0;
                for (k = 0; k <= n - 1; k++) d = d + r[k] * r[k];
                d = Math.Sqrt(d);
                if (d < eps)
                {
                    return x;
                }
                for (k = 0; k <= n - 1; k++)
                    p[k] = r[k] - beta * p[k];
                i = i + 1;
            }
            return x;
        }

        /// <summary>
        /// 用豪斯荷尔德(Householder) 变换求解线性最小二乘问题。返回时存放QR 分解式中的正交矩阵Q
        /// </summary>
        /// <param name="a">存放超定方程组的系数矩阵A, 返回时存放QR 分解式中的R 矩阵</param>
        /// <param name="b">存放方程组右端的常数向量，返回时前n 个分量存放方程组的最小二乘解</param>
        /// <returns>返回时存放QR 分解式中的正交矩阵Q</returns>
        public static double[,] GMQR(double[,] a, double[] b)
        {
            int m = a.GetLength(0);
            int n = a.GetLength(1);
            if (m < n)
            {
                throw new Exception("m应大于或等于n");
            }
            int i, j;
            double d;
            var q = Matrix.MAQR(a);                 //QR分解
            var c = new double[n];
            for (i = 0; i <= n - 1; i++)
            {
                d = 0.0;
                for (j = 0; j <= m - 1; j++)
                    d = d + q[j, i] * b[j];
                c[i] = d;
            }
            b[n - 1] = c[n - 1] / a[n - 1, n - 1];
            for (i = n - 2; i >= 0; i--)
            {
                d = 0.0;
                for (j = i + 1; j <= n - 1; j++)
                    d = d + a[i, j] * b[j];
                b[i] = (c[i] - d) / a[i, i];
            }
            return q;
        }

        /// <summary>
        /// 利用广义逆法求超定方程组AX = B 的最小二乘解。其中A 为m×n(m≥n) 的矩阵，且列线性无关。当m =n 时，即为求线性代数方程组的解。
        /// <para>返回超定方程组的最小二乘解</para>
        /// </summary>
        /// <param name="a">存放m×n的实矩阵A.返回时其对角线给出奇异值（以非递增次序排列），其余元素均为0.</param>
        /// <param name="b">存放超定方程组右端的常数向量</param>
        /// <param name="eps">给定的精度要求</param>
        /// <param name="aa">返回A的广义逆A+</param>
        /// <param name="u">返回左奇异向量U</param>
        /// <param name="v">返回右奇异向量V</param>
        /// <returns>返回超定方程组的最小二乘解</returns>
        public static double[] GMIV(double[,] a, double[] b, double eps, out double[,] aa, out double[,] u, out double[,] v)
        {
            int m = a.GetLength(0);
            int n = a.GetLength(1);
            int ka = Math.Max(m, n) + 1;
            int i, j;
            double[] x = new double[n];
            Matrix.GINV(a, eps, out aa, out u, out v);
            for (i = 0; i <= n - 1; i++)
            {
                x[i] = 0.0;
                for (j = 0; j <= m - 1; j++)
                    x[i] = x[i] + aa[i, j] * b[j];
            }
            return x;
        }

        /// <summary>
        /// 求解病态线性代数方程组AX=B 。返回方程组的解向量
        /// </summary>
        /// <param name="a">存放方程组的系数矩阵</param>
        /// <param name="b">存放方程组右端常数向量</param>
        /// <param name="eps">控制精度要求</param>
        /// <param name="i">尝试迭代次数</param>
        /// <returns>返回方程组的解向量</returns>
        public static double[] BINT(double[,] a, double[] b, double eps, int i = 60)
        {
            int n = a.GetLength(0);
            if (a.GetLength(1)!=n)
            {
                throw new Exception("a 非方阵");
            }
            if (b.Length != n)
            {
                throw new Exception("b 不匹配");
            }
            int j, k;
            double q, qq;
            double[,] p = new double[n, n];
            double[] r = new double[n];
            double[] e = new double[n];
            double[] x = new double[n];
            for (k = 0; k <= n - 1; k++)
                for (j = 0; j <= n - 1; j++) p[k,j] = a[k,j];
            for (k = 0; k <= n - 1; k++) r[k] = b[k];
            GAUS(p, x);
            for (k = 0; k < n; k++) x[k] = r[k];
            q = 1.0 + eps;
            while (q >= eps)
            {
                if (i == 0) throw new Exception("程序工作失败");
                i = i - 1;
                //
                for (k = 0; k < n; k++)
                {
                    e[k] = 0.0;
                    for (j = 0; j < n; j++)
                        e[k] = e[k] + a[k,j] * x[j];
                }
                //
                for (k = 0; k <= n - 1; k++) r[k] = b[k] - e[k];
                for (k = 0; k <= n - 1; k++)
                    for (j = 0; j <= n - 1; j++) p[k,j] = a[k,j];
                GAUS(p,r);
                q = 0.0;
                for (k = 0; k <= n - 1; k++)
                {
                    qq = Math.Abs(r[k]) / (1.0 + Math.Abs(x[k] + r[k]));
                    if (qq > q) q = qq;
                }
                for (k = 0; k <= n - 1; k++) x[k] = x[k] + r[k];
            }
            return x;
        }
    }
}
