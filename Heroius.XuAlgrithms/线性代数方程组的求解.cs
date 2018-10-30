using System;
using System.Collections.Generic;
using System.Text;

namespace Heroius.XuAlgrithms
{
    public static partial class Algrithms
    {
        /// <summary>
        /// 用全选主元高斯（Gauss）消去法求解n阶线性代数方程组AX=B.
        /// </summary>
        /// <param name="a">存放方程组的系数矩阵，返回时将被破坏</param>
        /// <param name="b">存放方程组右端的常数向量，返回方程组的解向量</param>
        public static void GAUS(double[,] a, double[] b)
        {
            int n = a.GetLength(0);
            if (a.GetLength(1)!= n)
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
                        t = Math.Abs(a[i,j]);
                        if (t > d) { d = t; js[k] = j; ss= i; }
                    }
                if (d + 1.0 == 1.0) l = 0;
                else
                {
                    if (js[k] != k)
                        for (i = 0; i <= n - 1; i++)
                        {
                            t = a[i,k];
                            a[i,k] = a[i,js[k]];
                            a[i,js[k]] = t;
                        }
                    if (ss!= k)
                    {
                        for (j = k; j <= n - 1; j++)
                        {
                            t = a[k,j];
                            a[k,j] = a[ss,j];
                            a[ss,j] = t;
                        }
                        for (j = 0; j <= m - 1; j++)
                        {
                            t = b[k,j];
                            b[k,j] = b[ss,j];
                            b[ss,j] = t;
                        }
                    }
                }
                if (l == 0)
                {
                    throw new Exception("系数矩阵奇异！无解. ");
                }
                d = a[k,k];
                for (j = k + 1; j <= n - 1; j++)
                    a[k,j] = a[k,j] / d;
                for (j = 0; j <= m - 1; j++)
                    b[k,j] = b[k,j] / d;
                for (j = k + 1; j <= n - 1; j++)
                    for (i = 0; i <= n - 1; i++)
                    {
                        if (i != k)
                            a[i,j] = a[i,j] - a[i,k] * a[k,j];
                    }
                for (j = 0; j <= m - 1; j++)
                    for (i = 0; i <= n - 1; i++)
                    {
                        if (i != k)
                            b[i,j] = b[i,j] - a[i,k] * b[k,j];
                    }
            }
            for (k = n - 1; k >= 0; k--)
                if (js[k] != k)
                    for (j = 0; j <= m - 1; j++)
                    {
                        t = b[k,j]; b[k,j] = b[js[k],j]; b[js[k],j] = t;
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
            if (ar.GetLength(1)!=n)
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
            int k, i, j,ss=0;
            double p, q, s, d;
            var js = new int[n];
            for (k = 0; k <= n - 2; k++)
            {
                d = 0.0;
                for (i = k; i <= n - 1; i++)
                    for (j = k; j <= n - 1; j++)
                    {
                        p = ar[i,j] * ar[i,j] + ai[i,j] * ai[i,j];
                        if (p > d) { d = p; js[k] = j;ss= i; }
                    }
                if (d + 1.0 == 1.0)
                {
                    throw new Exception("系数矩阵奇异！无解. ");
                }
                if (ss!= k)
                {
                    for (j = k; j <= n - 1; j++)
                    {
                        p = ar[k,j]; ar[k,j] = ar[ss,j]; ar[ss,j] = p;
                        p = ai[k,j]; ai[k,j] = ai[ss,j]; ai[ss,j] = p;
                    }
                    p = br[k]; br[k] = br[ss]; br[ss] = p;
                    p = bi[k]; bi[k] = bi[ss]; bi[ss] = p;
                }
                if (js[k] != k)
                    for (i = 0; i <= n - 1; i++)
                    {
                        p = ar[i,k]; ar[i,k] = ar[i,js[k]]; ar[i,js[k]] = p;
                        p = ai[i,k]; ai[i,k] = ai[i,js[k]]; ai[i,js[k]] = p;
                    }
                for (j = k + 1; j <= n - 1; j++)
                {
                    p = ar[k,j] * ar[k,k]; q = -ai[k,j] * ai[k,k];
                    s = (ar[k,k] - ai[k,k]) * (ar[k,j] + ai[k,j]);
                    ar[k,j] = (p - q) / d; ai[k,j] = (s - p - q) / d;
                }
                p = br[k] * ar[k,k]; q = -bi[k] * ai[k,k];
                s = (ar[k,k] - ai[k,k]) * (br[k] + bi[k]);
                br[k] = (p - q) / d; bi[k] = (s - p - q) / d;
                for (i = k + 1; i <= n - 1; i++)
                {
                    for (j = k + 1; j <= n - 1; j++)
                    {
                        p = ar[i,k] * ar[k,j]; q = ai[i,k] * ai[k,j];
                        s = (ar[i,k] + ai[i,k]) * (ar[k,j] + ai[k,j]);
                        ar[i,j] = ar[i,j] - p + q;
                        ai[i,j] = ai[i,j] - s + p + q;
                    }
                    p = ar[i,k] * br[k]; q = ai[i,k] * bi[k];
                    s = (ar[i,k] + ai[i,k]) * (br[k] + bi[k]);
                    br[i] = br[i] - p + q; bi[i] = bi[i] - s + p + q;
                }
            }
            d = ar[n - 1,n - 1] * ar[n - 1,n - 1] + ai[n - 1,n - 1] * ai[n - 1,n - 1];
            if (d + 1.0 == 1.0)
            {
                throw new Exception("系数矩阵奇异！无解. ");
            }
            p = ar[n - 1,n - 1] * br[n - 1]; q = -ai[n - 1,n - 1] * bi[n - 1];
            s = (ar[n - 1,n - 1] - ai[n - 1,n - 1]) * (br[n - 1] + bi[n - 1]);
            br[n - 1] = (p - q) / d; bi[n - 1] = (s - p - q) / d;
            for (i = n - 2; i >= 0; i--)
                for (j = i + 1; j <= n - 1; j++)
                {
                    p = ar[i,j] * br[j]; q = ai[i,j] * bi[j];
                    s = (ar[i,j] + ai[i,j]) * (br[j] + bi[j]);
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

            int k, i, j,ss=0;
            double p, q, s, d;
            var js = new int[n];
            for (k = 0; k <= n - 1; k++)
            {
                d = 0.0;
                for (i = k; i <= n - 1; i++)
                    for (j = k; j <= n - 1; j++)
                    {
                        p = ar[i,j] * ar[i,j] + ai[i,j] * ai[i,j];
                        if (p > d) { d = p; js[k] = j;ss= i; }
                    }
                if (d + 1.0 == 1.0)
                {
                    throw new Exception("系数矩阵奇异！无解. ");
                }
                if (ss!= k)
                {
                    for (j = k; j <= n - 1; j++)
                    {
                        p = ar[k,j]; ar[k,j] = ar[ss,j]; ar[ss,j] = p;
                        p = ai[k,j]; ai[k,j] = ai[ss,j]; ai[ss,j] = p;
                    }
                    for (j = 0; j <= m - 1; j++)
                    {
                        p = br[k,j]; br[k,j] = br[ss,j]; br[ss,j] = p;
                        p = bi[k,j]; bi[k,j] = bi[ss,j]; bi[ss,j] = p;
                    }
                }
                if (js[k] != k)
                    for (i = 0; i <= n - 1; i++)
                    {
                        p = ar[i,k]; ar[i,k] = ar[i,js[k]]; ar[i,js[k]] = p;
                        p = ai[i,k]; ai[i,k] = ai[i,js[k]]; ai[i,js[k]] = p;
                    }
                for (j = k + 1; j <= n - 1; j++)
                {
                    p = ar[k,j] * ar[k,k];
                    q = -ai[k,j] * ai[k,k];
                    s = (ar[k,k] - ai[k,k]) * (ar[k,j] + ai[k,j]);
                    ar[k,j] = (p - q) / d;
                    ai[k,j] = (s - p - q) / d;
                }
                for (j = 0; j <= m - 1; j++)
                {
                    p = br[k,j] * ar[k,k];
                    q = -bi[k,j] * ai[k,k];
                    s = (ar[k,k] - ai[k,k]) * (br[k,j] + bi[k,j]);
                    br[k,j] = (p - q) / d;
                    bi[k,j] = (s - p - q) / d;
                }
                for (i = 0; i <= n - 1; i++)
                    if (i != k)
                    {
                        for (j = k + 1; j <= n - 1; j++)
                        {
                            p = ar[i,k] * ar[k,j];
                            q = ai[i,k] * ai[k,j];
                            s = (ar[i,k] + ai[i,k]) * (ar[k,j] + ai[k,j]);
                            ar[i,j] = ar[i,j] - p + q;
                            ai[i,j] = ai[i,j] - s + p + q;
                        }
                        for (j = 0; j <= m - 1; j++)
                        {
                            p = ar[i,k] * br[k,j];
                            q = ai[i,k] * bi[k,j];
                            s = (ar[i,k] + ai[i,k]) * (br[k,j] + bi[k,j]);
                            br[i,j] = br[i,j] - p + q;
                            bi[i,j] = bi[i,j] - s + p + q;
                        }
                    }
            }
            for (k = n - 1; k >= 0; k--)
                if (js[k] != k)
                    for (j = 0; j <= m - 1; j++)
                    {
                        p = br[k,j]; br[k,j] = br[js[k],j]; br[js[k],j] = p;
                        p = bi[k,j]; bi[k,j] = bi[js[k],j]; bi[js[k],j] = p;
                    }
        }

        /// <summary>
        /// 用追赶法求解n 阶三对角线方程组AX=D 。
        /// </summary>
        /// <param name="b">以行为主存放三对角矩阵中3 条对角线上的元素</param>
        /// <param name="d">存放方程组右端的常数向社，返回方程组的解向拭</param>
        public static void TRDE(double[] b, double[] d)
        {
            int m = b.Length;
            int n = d.Length;
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


        public static void BAND()
        {

        }
    }
}
