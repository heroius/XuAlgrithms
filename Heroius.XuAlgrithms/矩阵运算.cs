using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Text;

namespace Heroius.XuAlgrithms
{
    /// <summary>
    /// 算法函数
    /// </summary>
    public static partial class Algrithms
    {
        /// <summary>
        /// 求m×n阶矩阵A与n×k阶矩阵B的乘积矩阵C=AB.
        /// </summary>
        /// <param name="A">存放矩阵A的元素</param>
        /// <param name="B">存放矩阵B的元素</param>
        /// <returns>m×k阶乘积矩阵C</returns>
        public static double[,] TRMUL(double[,] A, double[,] B)
        {
            int m = A.GetLength(0), n = A.GetLength(1);
            int k = B.GetLength(1);
            if (A.GetLength(1) != B.GetLength(0))
                throw new Exception("矩阵无法相乘");
            var C = new double[m, k];
            int i, j, t;
            for (i = 0; i < m; i++)
                for (j = 0; j < k; j++)
                {
                    C[i, j] = 0.0;
                    for (t = 0; t < n; t++)
                        C[i, j] = C[i, j] + A[i, t] * B[t, j];
                }
            return C;
        }

        /// <summary>
        /// 求m×n阶矩阵A与n×k阶矩阵B的乘积矩阵C=AB.
        /// </summary>
        /// <param name="ar">矩阵A的实部</param>
        /// <param name="ai">矩阵A的虚部</param>
        /// <param name="br">矩阵B的实部</param>
        /// <param name="bi">矩阵B的虚部</param>
        /// <param name="cr">输出：m*k阶乘积矩阵的实部</param>
        /// <param name="ci">输出：m*k阶乘积矩阵的虚部</param>
        public static void TCMUL(double[,] ar, double[,] ai, double[,] br, double[,] bi, out double[,] cr, out double[,] ci)
        {
            int m = ar.GetLength(0);
            int n = ar.GetLength(1);
            if (ai.GetLength(0) != m || ai.GetLength(1) != n)
            {
                throw new Exception("A实部与虚部矩阵不兼容");
            }
            if (br.GetLength(0) != n)
            {
                throw new Exception("A B 不可相乘");
            }
            int k = br.GetLength(1);
            if (bi.GetLength(0) != n || bi.GetLength(1) != k)
            {
                throw new Exception("B实部与虚部矩阵不兼容");
            }

            cr = new double[m, k];
            ci = new double[m, k];

            int i, j, l;
            double p, q, s;
            for (i = 0; i < m; i++)
                for (j = 0; j < k; j++)
                {
                    cr[i, j] = 0.0; ci[i, j] = 0.0;
                    for (l = 0; l < n; l++)
                    {
                        p = ar[i, l] * br[l, j];
                        q = ai[i, l] * bi[l, j];
                        s = (ar[i, l] + ai[i, l]) * (br[l, j] + bi[l, j]);
                        cr[i, j] = cr[i, j] + p - q;
                        ci[i, j] = ci[i, j] + s - p - q;
                    }
                }
        }

        /// <summary>
        /// 用全选主元高斯-约当(Gauss-Jordan)消去法求n阶实矩阵A的逆矩阵
        /// </summary>
        /// <param name="a">存放矩阵A，返回时存放其逆矩阵</param>
        public static void RINV(double[,] a)
        {
            int n = a.GetLength(0);
            if (a.GetLength(1) != n)
            {
                throw new Exception("输入非n阶矩阵");
            }
            int[] s, js;
            int i, j, k;
            double d, p;
            s = new int[n];
            js = new int[n];
            for (k = 0; k <= n - 1; k++)
            {
                d = 0.0;
                for (i = k; i <= n - 1; i++)
                    for (j = k; j <= n - 1; j++)
                    {
                        p = Math.Abs(a[i, j]);
                        if (p > d) { d = p; s[k] = i; js[k] = j; }
                    }
                if (d + 1.0 == 1.0)
                {
                    throw new Exception("A为奇异矩阵！没有逆矩阵. ");
                }
                if (s[k] != k)
                    for (j = 0; j <= n - 1; j++)
                    {
                        p = a[k, j]; a[k, j] = a[s[k], j]; a[s[k], j] = p;
                    }
                if (js[k] != k)
                    for (i = 0; i <= n - 1; i++)
                    {
                        p = a[i, k]; a[i, k] = a[i, js[k]]; a[i, js[k]] = p;
                    }
                a[k, k] = 1.0 / a[k, k];
                for (j = 0; j <= n - 1; j++)
                    if (j != k) a[k, j] = a[k, j] * a[k, k];
                for (i = 0; i <= n - 1; i++)
                    if (i != k)
                        for (j = 0; j <= n - 1; j++)
                            if (j != k) a[i, j] = a[i, j] - a[i, k] * a[k, j];
                for (i = 0; i <= n - 1; i++)
                    if (i != k) a[i, k] = -a[i, k] * a[k, k];
            }
            for (k = n - 1; k >= 0; k--)
            {
                if (js[k] != k)
                    for (j = 0; j <= n - 1; j++)
                    {
                        p = a[k, j]; a[k, j] = a[js[k], j]; a[js[k], j] = p;
                    }
                if (s[k] != k)
                    for (i = 0; i <= n - 1; i++)
                    {
                        p = a[i, k]; a[i, k] = a[i, s[k]]; a[i, s[k]] = p;
                    }
            }
        }

        /// <summary>
        /// 用全选主元高斯-约当(Gauss-Jordan)消去法求n阶复矩阵A=AR+jAI的逆矩阵
        /// </summary>
        /// <param name="ar">存放矩阵A实部，返回时存放逆矩阵的实部</param>
        /// <param name="ai">存放矩阵A虚部，返回时存放逆矩阵的虚部</param>
        public static void CINV(double[,] ar, double[,] ai)
        {
            int n = ar.GetLength(0);
            if (ar.GetLength(1) != n)
            {
                throw new Exception("输入非n阶矩阵");
            }
            if (ai.GetLength(0) != n || ai.GetLength(1) != n)
            {
                throw new Exception("实部与虚部矩阵不兼容");
            }

            int[] ss, js; int i, j, k;
            double p, q, s, t, d, b;

            ss = new int[n];
            js = new int[n];
            for (k = 0; k <= n - 1; k++)
            {
                d = 0.0;
                for (i = k; i <= n - 1; i++)
                    for (j = k; j <= n - 1; j++)
                    {

                        p = ar[i, j] * ar[i, j] + ai[i, j] * ai[i, j];
                        if (p > d) { d = p; ss[k] = i; js[k] = j; }
                    }
                if (d + 1.0 == 1.0)
                {
                    throw new Exception("矩阵奇异！没有逆矩阵.");
                }
                if (ss[k] != k)
                    for (j = 0; j <= n - 1; j++)
                    {

                        t = ar[k, j]; ar[k, j] = ar[ss[k], j]; ar[ss[k], j] = t;
                        t = ai[k, j]; ai[k, j] = ai[ss[k], j]; ai[ss[k], j] = t;
                    }
                if (js[k] != k)
                    for (i = 0; i <= n - 1; i++)
                    {
                        t = ar[i, k]; ar[i, k] = ar[i, js[k]]; ar[i, js[k]] = t;
                        t = ai[i, k]; ai[i, k] = ai[i, js[k]]; ai[i, js[k]] = t;
                    }
                ar[k, k] = ar[k, k] / d; ai[k, k] = -ai[k, k] / d;
                for (j = 0; j <= n - 1; j++)
                    if (j != k)
                    {
                        p = ar[k, j] * ar[k, k]; q = ai[k, j] * ai[k, k];
                        s = (ar[k, j] + ai[k, j]) * (ar[k, k] + ai[k, k]);
                        ar[k, j] = p - q; ai[k, j] = s - p - q;
                    }
                for (i = 0; i <= n - 1; i++)
                    if (i != k)
                    {
                        for (j = 0; j <= n - 1; j++)
                            if (j != k)
                            {
                                p = ar[k, j] * ar[i, k]; q = ai[k, j] * ai[i, k];
                                s = (ar[k, j] + ai[k, j]) * (ar[i, k] + ai[i, k]);
                                t = p - q; b = s - p - q;
                                ar[i, j] = ar[i, j] - t;
                                ai[i, j] = ai[i, j] - b;
                            }
                    }
                for (i = 0; i <= n - 1; i++)
                    if (i != k)
                    {
                        p = ar[i, k] * ar[k, k]; q = ai[i, k] * ai[k, k];
                        s = (ar[i, k] + ai[i, k]) * (ar[k, k] + ai[k, k]);
                        ar[i, k] = q - p; ai[i, k] = p + q - s;
                    }
            }
            for (k = n - 1; k >= 0; k--)
            {
                if (js[k] != k)
                    for (j = 0; j <= n - 1; j++)
                    {
                        t = ar[k, j]; ar[k, j] = ar[js[k], j]; ar[js[k], j] = t;
                        t = ai[k, j]; ai[k, j] = ai[js[k], j]; ai[js[k], j] = t;
                    }
                if (ss[k] != k)
                    for (i = 0; i <= n - 1; i++)
                    {
                        t = ar[i, k]; ar[i, k] = ar[i, ss[k]]; ar[i, ss[k]] = t;
                        t = ai[i, k]; ai[i, k] = ai[i, ss[k]]; ai[i, ss[k]] = t;
                    }
            }
        }

        /// <summary>
        /// 求n阶对称正定矩阵A的逆矩阵
        /// </summary>
        /// <param name="a">对称正定矩阵A，返回时存放逆矩阵。</param>
        public static void SSGJ(double[,] a)
        {
            int n = a.GetLength(0);
            if (a.GetLength(1) != n)
            {
                throw new Exception("非n阶矩阵");
            }

            int i, j, k, m;
            double w, g;
            double[] b = new double[n];
            for (k = 0; k <= n - 1; k++)
            {
                w = a[0, 0];
                if (Math.Abs(w) + 1.0 == 1.0)
                {
                    throw new Exception("矩阵奇异！");
                }
                m = n - k - 1;
                for (i = 1; i <= n - 1; i++)
                {
                    g = a[i, 0]; b[i] = g / w;
                    if (i <= m) b[i] = -b[i];
                    for (j = 1; j <= i; j++)
                        a[i - 1, j - 1] = a[i, j] + g * b[j];
                }
                a[n - 1, n - 1] = 1.0 / w;
                for (i = 1; i <= n - 1; i++)
                    a[n - 1, i - 1] = b[i];
            }
            for (i = 0; i <= n - 2; i++)
                for (j = i + 1; j <= n - 1; j++)
                    a[i, j] = a[j, i];
        }

        /// <summary>
        /// 用特兰持(Trench)方法求托伯利兹(Toeplitz)矩阵的逆矩阵。
        /// </summary>
        /// <param name="t">存放T型矩阵钟的元素t0,t1...tn-1</param>
        /// <param name="tt">后n-1个元素存放T型矩阵中的元素τ0,τ1...τn-1</param>
        /// <returns>T型矩阵的逆矩阵</returns>
        public static double[,] TRCH(double[] t, double[] tt)
        {
            int n = t.Length;
            if (tt.Length != n)
            {
                throw new Exception("");
            }
            double[,] b = new double[n, n];
            int i, j, k;
            double a, s; double[] c, r, p;
            c = new double[n];
            r = new double[n];
            p = new double[n];
            if (Math.Abs(t[0]) + 1.0 == 1.0)
            {
                throw new Exception("求逆失败！");
            }
            a = t[0]; c[0] = tt[1] / t[0]; r[0] = t[1] / t[0];
            for (k = 0; k <= n - 3; k++)
            {
                s = 0.0;
                for (j = 1; j <= k + 1; j++) s = s + c[k + 1 - j] * tt[j];
                s = (s - tt[k + 2]) / a;
                for (i = 0; i <= k; i++) p[i] = c[i] + s * r[k - i];
                c[k + 1] = -s;
                s = 0.0;
                for (j = 1; j <= k + 1; j++) s = s + r[k + 1 - j] * t[j];
                s = (s - t[k + 2]) / a;
                for (i = 0; i <= k; i++)
                {
                    r[i] = r[i] + s * c[k - i];
                    c[k - i] = p[k - i];
                }
                r[k + 1] = -s;
                a = 0.0;
                for (j = 1; j <= k + 2; j++) a = a + t[j] * c[j - 1];
                a = t[0] - a;
                if (Math.Abs(a) + 1.0 == 1.0)
                {
                    throw new Exception("求逆失败！");
                }
            }
            b[0, 0] = 1.0 / a;
            for (i = 0; i <= n - 2; i++)
            {
                b[0, i + 1] = -r[i] / a; b[i + 1, 0] = -c[i] / a;
            }
            for (i = 0; i <= n - 2; i++)
                for (j = 0; j <= n - 2; j++)
                {
                    b[i + 1, j + 1] = b[i, j] - c[i] * b[0, j + 1];
                    b[i + 1, j + 1] = b[i + 1, j + 1] + c[n - j - 2] * b[0, n - i - 1];
                }
            return b;
        }

        /// <summary>
        /// 用全选主元高斯(Gauss)消去法计算n阶方阵A所对应的行列式值。
        /// </summary>
        /// <param name="a">存放方阵A的元素，返回时被破坏</param>
        /// <returns>行列式的值</returns>
        public static double SDET(double[,] a)
        {
            int n = a.GetLength(0);
            if (a.GetLength(1) != n)
            {
                throw new Exception("非方阵");
            }

            int i, j, k, ss = 0, js = 0;
            double f, q, d;
            f = 1.0;
            double det = 1.0;
            for (k = 0; k <= n - 2; k++)
            {
                q = 0.0;
                for (i = k; i <= n - 1; i++)
                    for (j = k; j <= n - 1; j++)
                    {
                        d = Math.Abs(a[i, j]);
                        if (d > q) { q = d; ss = i; js = j; }
                    }
                if (q + 1.0 == 1.0) det = 0.0;
                else
                {
                    if (ss != k)
                    {
                        f = -f;
                        for (j = k; j <= n - 1; j++)
                        {
                            d = a[k, j]; a[k, j] = a[ss, j]; a[ss, j] = d;
                        }
                    }
                    if (js != k)
                    {
                        f = -f;
                        for (i = k; i <= n - 1; i++)
                        {
                            d = a[i, js]; a[i, js] = a[i, k]; a[i, k] = d;
                        }
                    }
                    det = det * a[k, k];
                    for (i = k + 1; i <= n - 1; i++)
                    {
                        d = a[i, k] / a[k, k];
                        for (j = k + 1; j <= n - 1; j++)
                        {
                            a[i, j] = a[i, j] - d * a[k, j];
                        }
                    }
                }
            }
            det = f * det * a[n - 1, n - 1];
            return det;
        }

        /// <summary>
        /// 用全选主元高斯(Gauss)消去法计算矩阵A的秩
        /// </summary>
        /// <param name="a">矩阵A的元素，返回时被破坏</param>
        /// <returns>矩阵A的秩</returns>
        public static int RANK(double[,] a)
        {
            int m = a.GetLength(0);
            int n = a.GetLength(1);
            int i, j, nn, ss = 0, js = 0, l;
            double q, d;
            nn = m;
            if (m >= n) nn = n;
            int k = 0;
            for (l = 0; l <= nn - 1; l++)
            {
                q = 0.0;
                for (i = l; i <= m - 1; i++)
                    for (j = l; j <= n - 1; j++)
                    {
                        d = Math.Abs(a[i, j]);
                        if (d > q) { q = d; ss = i; js = j; }
                    }
                if (q + 1.0 == 1.0) break;
                k = k + 1;
                if (ss != l)
                {
                    for (j = l; j <= n - 1; j++)
                    {
                        d = a[l, j]; a[l, j] = a[ss, j]; a[ss, j] = d;
                    }
                }
                if (js != l)
                {
                    for (i = l; i <= m - 1; i++)
                    {
                        d = a[i, js]; a[i, js] = a[i, l]; a[i, l] = d;
                    }
                }
                for (i = l + 1; i <= n - 1; i++)
                {
                    d = a[i, l] / a[l, l];
                    for (j = l + 1; j <= n - 1; j++)
                    {
                        a[i, j] = a[i, j] - d * a[l, j];
                    }
                }
            }
            return k;
        }

        /// <summary>
        /// 用乔里斯基(Cholesky)分解法求对称正定矩阵的三角分解，并求行列式值。
        /// <para>对n阶对称正定矩阵A，存在一个实的非奇异下三角矩阵L，使A=LLT</para>
        /// </summary>
        /// <param name="a">存放对称正定矩阵A.返回时其下三角部分存放分解得到的下三角阵L，其余元素均为0.</param>
        /// <returns>行列式的值</returns>
        public static double CHOL(double[,] a)
        {
            int n = a.GetLength(0);
            if (a.GetLength(1) != n)
            {
                throw new Exception("非对称正定矩阵");
            }
            int i, j, k;
            double d;
            double det;
            if ((a[0,0] + 1.0 == 1.0) || (a[0,0] < 0.0))
            {
                throw new Exception("分解失败！");
            }
            a[0,0] = Math.Sqrt(a[0,0]);
            d = a[0,0];
            for (i = 1; i <= n - 1; i++) a[i,0] = a[i,0] / a[0,0];
            for (j = 1; j <= n - 1; j++)
            {
                for (k = 0; k <= j - 1; k++)
                    a[j,j] = a[j,j] - a[j,k] * a[j,k];
                if ((a[j,j] + 1.0 == 1.0) || (a[j,j] < 0.0))
                {
                    throw new Exception("分解失败！");
                }
                a[j,j] = Math.Sqrt(a[j,j]);
                d = d * a[j,j];
                for (i = j + 1; i <= n - 1; i++)
                {
                    for (k = 0; k <= j - 1; k++)
                        a[i,j] = a[i,j] - a[i,k] * a[j,k];
                    a[i,j] = a[i,j] / a[j,j];
                }
            }
            det = d * d;
            for (i = 0; i <= n - 2; i++)
                for (j = i + 1; j <= n - 1; j++) a[i,j] = 0.0;
            return det;
        }

        /// <summary>
        /// 对n阶实矩阵A进行LU分解。即A=LU，L下三角，U上三角。
        /// </summary>
        /// <param name="a">n阶矩阵A，返回时存放Q矩阵。Q=L+U-In</param>
        /// <param name="l">返回下三角矩阵L</param>
        /// <param name="u">返回上三角矩阵U</param>
        public static void LLUU(double[,] a, out double[,] l, out double[,] u)
        {
            int n = a.GetLength(0);
            if (a.GetLength(1)!=n)
            {
                throw new Exception("非矩阵");
            }
            l = new double[n, n];
            u = new double[n, n];
            int i, j, k;
            for (k = 0; k <= n - 2; k++)
            {
                if (Math.Abs(a[k,k]) + 1.0 == 1.0)
                {
                    throw new Exception("分解失败！");
                }
                for (i = k + 1; i <= n - 1; i++) a[i,k] = a[i,k] / a[k,k];
                for (i = k + 1; i <= n - 1; i++)
                {
                    for (j = k + 1; j <= n - 1; j++)
                    {
                        a[i,j] = a[i,j] - a[i,k] * a[k,j];
                    }
                }
            }
            for (i = 0; i <= n - 1; i++)
            {
                for (j = 0; j < i; j++)
                { l[i,j] = a[i,j]; u[i,j] = 0.0; }
                l[i,i] = 1.0; u[i,i] = a[i,i];
                for (j = i + 1; j <= n - 1; j++)
                { l[i,j] = 0.0; u[i,j] = a[i,j]; }
            }
        }

        /// <summary>
        /// 用豪斯荷尔德(Householder)变换对一般m×n阶的实矩阵进行QR分解。
        /// <para>Q为m×m的正交矩阵，R为m×n的上三角矩阵。</para>
        /// </summary>
        /// <param name="a">存放m×n的实矩阵A，返回时其右上三角部分存放QR分解中的上三角矩阵R.</param>
        /// <returns>QR分解中的正交矩阵Q</returns>
        public static double[,] MAQR(double[,] a)
        {
            int m = a.GetLength(0);
            int n = a.GetLength(1);
            double[,] q = new double[m, m];
            int i, j, k, nn, jj;
            double u, alpha, w, t;
            if (m < n)
            {
                throw new Exception("QR分解失败！");
            }
            for (i = 0; i <= m - 1; i++)
                for (j = 0; j <= m - 1; j++)
                {
                    q[i,j] = 0.0;
                    if (i == j) q[i,j] = 1.0;
                }
            nn = n;
            if (m == n) nn = m - 1;
            for (k = 0; k <= nn - 1; k++)
            {

                u = 0.0;
                for (i = k; i <= m - 1; i++)
                {
                    w = Math.Abs(a[i,k]);
                    if (w > u) u = w;
                }
                alpha = 0.0;
                for (i = k; i <= m - 1; i++)
                { t = a[i,k] / u; alpha = alpha + t * t; }
                if (a[k,k] > 0.0) u = -u;
                alpha = u * Math.Sqrt(alpha);
                if (Math.Abs(alpha) + 1.0 == 1.0)
                {
                    throw new Exception("QR分解失败！");
                }
                u = Math.Sqrt(2.0 * alpha * (alpha - a[k,k]));
                if ((u + 1.0) != 1.0)
                {
                    a[k,k] = (a[k,k] - alpha) / u;
                    for (i = k + 1; i <= m - 1; i++) a[i,k] = a[i,k] / u;
                    for (j = 0; j <= m - 1; j++)
                    {
                        t = 0.0;
                        for (jj = k; jj <= m - 1; jj++)
                            t = t + a[jj,k] * q[jj,j];
                        for (i = k; i <= m - 1; i++)
                            q[i,j] = q[i,j] - 2.0 * t * a[i,k];
                    }
                    for (j = k + 1; j <= n - 1; j++)
                    {
                        t = 0.0;
                        for (jj = k; jj <= m - 1; jj++)
                            t = t + a[jj,k] * a[jj,j];
                        for (i = k; i <= m - 1; i++)
                            a[i,j] = a[i,j] - 2.0 * t * a[i,k];
                    }
                    a[k,k] = alpha;
                    for (i = k + 1; i <= m - 1; i++) a[i,k] = 0.0;
                }
            }
            for (i = 0; i <= m - 2; i++)
                for (j = i + 1; j <= m - 1; j++)
                {
                    t = q[i,j]; q[i,j] = q[j,i]; q[j,i] = t;
                }
            return q;
        }

        /// <summary>
        /// 用豪斯荷尔德(Householder)变换以及变形QR算法对一般实矩阵A进行奇异值分解。
        /// <para>   |E 0| T</para>
        /// <para>A=U|   |V</para>
        /// <para>   |0 0|</para>
        /// </summary>
        /// <param name="a">m×n阶实矩阵，返回时其对角线给出奇异值（以非递增次序排列），其余元素为0.</param>
        /// <param name="eps">给定的精度要求</param>
        /// <param name="u">返回左奇异向量U</param>
        /// <param name="v">返回右奇异向量V</param>
        public static void MUAV(double[,] a, double eps, out double[,] u, out double[,] v)
        {
            int m = a.GetLength(0);
            int n = a.GetLength(1);
            int ka = Math.Max(m, n) + 1;
            double[] s = new double[ka];
            double[] e = new double[ka];
            double[] w = new double[ka];
            double[] fg = new double[2];
            double[] cs = new double[2];
            u = new double[m, m];
            v = new double[n, n];
            int i, j, k, l, it, ll, kk, mm, nn, m1, ks;
            double d, dd, t, sm, sm1, em1, sk, ek, b, c, shh;
            it = 60; k = n;
            if (m - 1 < n) k = m - 1;
            l = m;
            if (n - 2 < m) l = n - 2;
            if (l < 0) l = 0;
            ll = k;
            if (l > k) ll = l;
            if (ll >= 1)
            {
                for (kk = 1; kk <= ll; kk++)
                {
                    if (kk <= k)
                    {
                        d = 0.0;
                        for (i = kk; i <= m; i++)
                            d = d + a[i - 1,kk - 1] * a[i - 1,kk - 1];
                        s[kk - 1] = Math.Sqrt(d);
                        if (s[kk - 1] != 0.0)
                        {
                            if (a[kk - 1,kk - 1] != 0.0)
                            {
                                s[kk - 1] = Math.Abs(s[kk - 1]);
                                if (a[kk - 1,kk - 1] < 0.0) s[kk - 1] = -s[kk - 1];
                            }
                            for (i = kk; i <= m; i++)
                            {
                                a[i - 1,kk - 1] = a[i - 1,kk - 1] / s[kk - 1];
                            }
                            a[kk - 1,kk - 1] = 1.0 + a[kk - 1,kk - 1];
                        }
                        s[kk - 1] = -s[kk - 1];
                    }
                    if (n >= kk + 1)
                    {
                        for (j = kk + 1; j <= n; j++)
                        {
                            if ((kk <= k) && (s[kk - 1] != 0.0))
                            {
                                d = 0.0;
                                for (i = kk; i <= m; i++)
                                {
                                    d = d + a[i - 1,kk - 1] * a[i - 1,j - 1];
                                }
                                d = -d / a[kk - 1,kk - 1];
                                for (i = kk; i <= m; i++)
                                {
                                    a[i - 1,j - 1] = a[i - 1,j - 1] + d * a[i - 1,kk - 1];
                                }
                            }
                            e[j - 1] = a[kk - 1,j - 1];
                        }
                    }
                    if (kk <= k)
                    {
                        for (i = kk; i <= m; i++)
                        {
                            u[i - 1,kk - 1] = a[i - 1,kk - 1];
                        }
                    }
                    if (kk <= l)
                    {
                        d = 0.0;
                        for (i = kk + 1; i <= n; i++)
                            d = d + e[i - 1] * e[i - 1];
                        e[kk - 1] = Math.Sqrt(d);
                        if (e[kk - 1] != 0.0)
                        {
                            if (e[kk] != 0.0)
                            {
                                e[kk - 1] = Math.Abs(e[kk - 1]);
                                if (e[kk] < 0.0) e[kk - 1] = -e[kk - 1];
                            }
                            for (i = kk + 1; i <= n; i++)
                                e[i - 1] = e[i - 1] / e[kk - 1];
                            e[kk] = 1.0 + e[kk];
                        }
                        e[kk - 1] = -e[kk - 1];
                        if ((kk + 1 <= m) && (e[kk - 1] != 0.0))
                        {
                            for (i = kk + 1; i <= m; i++) w[i - 1] = 0.0;
                            for (j = kk + 1; j <= n; j++)
                                for (i = kk + 1; i <= m; i++)
                                    w[i - 1] = w[i - 1] + e[j - 1] * a[i - 1,j - 1];
                            for (j = kk + 1; j <= n; j++)
                                for (i = kk + 1; i <= m; i++)
                                {
                                    a[i - 1,j - 1] = a[i - 1,j - 1] - w[i - 1] * e[j - 1] / e[kk];
                                }
                        }
                        for (i = kk + 1; i <= n; i++)
                            v[i - 1,kk - 1] = e[i - 1];
                    }
                }
            }
            mm = n;
            if (m + 1 < n) mm = m + 1;
            if (k < n) s[k] = a[k,k];
            if (m < mm) s[mm - 1] = 0.0;
            if (l + 1 < mm) e[l] = a[l,mm - 1];
            e[mm - 1] = 0.0;
            nn = m;
            if (m > n) nn = n;
            if (nn >= k + 1)
            {
                for (j = k + 1; j <= nn; j++)
                {
                    for (i = 1; i <= m; i++) u[i - 1,j - 1] = 0.0;
                    u[j - 1,j - 1] = 1.0;
                }
            }
            if (k >= 1)
            {
                for (ll = 1; ll <= k; ll++)
                {
                    kk = k - ll + 1;
                    if (s[kk - 1] != 0.0)
                    {
                        if (nn >= kk + 1)
                            for (j = kk + 1; j <= nn; j++)
                            {
                                d = 0.0;
                                for (i = kk; i <= m; i++)
                                {
                                    d = d + u[i - 1,kk - 1] * u[i - 1,j - 1] / u[kk - 1,kk - 1];
                                }
                                d = -d;
                                for (i = kk; i <= m; i++)
                                {
                                    u[i - 1,j - 1] = u[i - 1,j - 1] + d * u[i - 1,kk - 1];
                                }
                            }
                        for (i = kk; i <= m; i++)
                            u[i - 1,kk - 1] = -u[i - 1,kk - 1];
                        u[kk - 1,kk - 1] = 1.0 + u[kk - 1,kk - 1];
                        if (kk - 1 >= 1)
                            for (i = 1; i <= kk - 1; i++) u[i - 1,kk - 1] = 0.0;
                    }
                    else
                    {
                        for (i = 1; i <= m; i++) u[i - 1,kk - 1] = 0.0;
                        u[kk - 1,kk - 1] = 1.0;
                    }
                }
            }
            for (ll = 1; ll <= n; ll++)
            {
                kk = n - ll + 1;
                if ((kk <= l) && (e[kk - 1] != 0.0))
                {
                    for (j = kk + 1; j <= n; j++)
                    {
                        d = 0.0;
                        for (i = kk + 1; i <= n; i++)
                        {
                            d = d + v[i - 1,kk - 1] * v[i - 1,j - 1] / v[kk,kk - 1];
                        }
                        d = -d;
                        for (i = kk + 1; i <= n; i++)
                        {
                            v[i - 1,j - 1] = v[i - 1,j - 1] + d * v[i - 1,kk - 1];
                        }
                    }
                }
                for (i = 1; i <= n; i++) v[i - 1,kk - 1] = 0.0;
                v[kk - 1,kk - 1] = 1.0;
            }
            for (i = 1; i <= m; i++)
                for (j = 1; j <= n; j++) a[i - 1,j - 1] = 0.0;
            m1 = mm; it = 60;
            while (1 == 1)
            {
                if (mm == 0)
                {
                    PPP(m, n, a, s, e, v);
                    return;
                }
                if (it == 0)
                {
                    PPP(m, n, a, s, e, v);
                    throw new Exception("程序工作失败！");
                }
                kk = mm - 1;
                while ((kk != 0) && (Math.Abs(e[kk - 1]) != 0.0))
                {
                    d = Math.Abs(s[kk - 1]) + Math.Abs(s[kk]);
                    dd = Math.Abs(e[kk - 1]);
                    if (dd > eps * d) kk = kk - 1;
                    else e[kk - 1] = 0.0;
                }
                if (kk == mm - 1)
                {
                    kk = kk + 1;
                    if (s[kk - 1] < 0.0)
                    {
                        s[kk - 1] = -s[kk - 1];
                        for (i = 1; i <= n; i++)
                            v[i - 1,kk - 1] = -v[i - 1,kk - 1];
                    }
                    while ((kk != m1) && (s[kk - 1] < s[kk]))
                    {
                        d = s[kk - 1]; s[kk - 1] = s[kk]; s[kk] = d;
                        if (kk < n)
                            for (i = 1; i <= n; i++)
                            {
                                d = v[i - 1,kk - 1];
                                v[i - 1,kk - 1] = v[i - 1,kk];
                                v[i - 1,kk] = d;
                            }
                        if (kk < m)
                            for (i = 1; i <= m; i++)
                            {
                                d = u[i - 1,kk - 1];
                                u[i - 1,kk - 1] = u[i - 1,kk];
                                u[i - 1,kk] = d;
                            }
                        kk = kk + 1;
                    }
                    it = 60;
                    mm = mm - 1;
                }
                else
                {
                    ks = mm;
                    while ((ks > kk) && (Math.Abs(s[ks - 1]) != 0.0))
                    {
                        d = 0.0;
                        if (ks != mm) d = d + Math.Abs(e[ks - 1]);
                        if (ks != kk + 1) d = d + Math.Abs(e[ks - 2]);
                        dd = Math.Abs(s[ks - 1]);
                        if (dd > eps * d) ks = ks - 1;
                        else s[ks - 1] = 0.0;
                    }
                    if (ks == kk)
                    {
                        kk = kk + 1;
                        d = Math.Abs(s[mm - 1]);
                        t = Math.Abs(s[mm - 2]);
                        if (t > d) d = t;
                        t = Math.Abs(e[mm - 2]);
                        if (t > d) d = t;
                        t = Math.Abs(s[kk - 1]);
                        if (t > d) d = t;
                        t = Math.Abs(e[kk - 1]);
                        if (t > d) d = t;
                        sm = s[mm - 1] / d; sm1 = s[mm - 2] / d;
                        em1 = e[mm - 2] / d;
                        sk = s[kk - 1] / d; ek = e[kk - 1] / d;
                        b = ((sm1 + sm) * (sm1 - sm) + em1 * em1) / 2.0;
                        c = sm * em1; c = c * c; shh = 0.0;
                        if ((b != 0.0) || (c != 0.0))
                        {
                            shh = Math.Sqrt(b * b + c);
                            if (b < 0.0) shh = -shh;
                            shh = c / (b + shh);
                        }
                        fg[0] = (sk + sm) * (sk - sm) - shh;
                        fg[1] = sk * ek;
                        for (i = kk; i <= mm - 1; i++)
                        {
                            SSS(fg, cs);
                            if (i != kk) e[i - 2] = fg[0];
                            fg[0] = cs[0] * s[i - 1] + cs[1] * e[i - 1];
                            e[i - 1] = cs[0] * e[i - 1] - cs[1] * s[i - 1];
                            fg[1] = cs[1] * s[i];
                            s[i] = cs[0] * s[i];
                            if ((cs[0] != 1.0) || (cs[1] != 0.0))
                                for (j = 1; j <= n; j++)
                                {
                                    d = cs[0] * v[j - 1,i - 1] + cs[1] * v[j - 1,i];
                                    v[j - 1,i] = -cs[1] * v[j - 1,i - 1] + cs[0] * v[j - 1,i];
                                    v[j - 1,i - 1] = d;
                                }
                            SSS(fg, cs);
                            s[i - 1] = fg[0];
                            fg[0] = cs[0] * e[i - 1] + cs[1] * s[i];
                            s[i] = -cs[1] * e[i - 1] + cs[0] * s[i];
                            fg[1] = cs[1] * e[i];
                            e[i] = cs[0] * e[i];
                            if (i < m)
                                if ((cs[0] != 1.0) || (cs[1] != 0.0))
                                    for (j = 1; j <= m; j++)
                                    {
                                        d = cs[0] * u[j - 1,i - 1] + cs[1] * u[j - 1,i];
                                        u[j - 1,i] = -cs[1] * u[j - 1,i - 1] + cs[0] * u[j - 1,i];
                                        u[j - 1,i - 1] = d;
                                    }
                        }
                        e[mm - 2] = fg[0];
                        it = it - 1;
                    }
                    else
                    {
                        if (ks == mm)
                        {
                            kk = kk + 1;
                            fg[1] = e[mm - 2]; e[mm - 2] = 0.0;
                            for (ll = kk; ll <= mm - 1; ll++)
                            {
                                i = mm + kk - ll - 1;
                                fg[0] = s[i - 1];
                                SSS(fg, cs);
                                s[i - 1] = fg[0];
                                if (i != kk)
                                {
                                    fg[1] = -cs[1] * e[i - 2];
                                    e[i - 2] = cs[0] * e[i - 2];
                                }
                                if ((cs[0] != 1.0) || (cs[1] != 0.0))
                                    for (j = 1; j <= n; j++)
                                    {
                                        d = cs[0] * v[j - 1,i - 1] + cs[1] * v[j - 1,mm - 1];
                                        v[j - 1,mm - 1] = -cs[1] * v[j - 1,i - 1] + cs[0] * v[j - 1,mm - 1];
                                        v[j - 1,i - 1] = d;
                                    }
                            }
                        }
                        else
                        {
                            kk = ks + 1;
                            fg[1] = e[kk - 2];
                            e[kk - 2] = 0.0;
                            for (i = kk; i <= mm; i++)
                            {
                                fg[0] = s[i - 1];
                                SSS(fg, cs);
                                s[i - 1] = fg[0];
                                fg[1] = -cs[1] * e[i - 1];
                                e[i - 1] = cs[0] * e[i - 1];
                                if ((cs[0] != 1.0) || (cs[1] != 0.0))
                                    for (j = 1; j <= m; j++)
                                    {
                                        d = cs[0] * u[j - 1,i - 1] + cs[1] * u[j - 1,kk - 2];
                                        u[j - 1,kk - 2] = -cs[1] * u[j - 1,i - 1] + cs[0] * u[j - 1,kk - 2];
                                        u[j - 1,i - 1] = d;
                                    }
                            }
                        }
                    }
                }
            }
        }

        /// <summary>
        /// 奇异值分解的子程序
        /// </summary>
        static void PPP(int m, int n, double[,] a, double[] s, double[] e, double[,] v)
        {
            int i, j;
            double d;
            if (m >= n) i = n;
            else i = m;
            for (j = 1; j <= i - 1; j++)
            {
                a[j - 1,j - 1] = s[j - 1];
                a[j - 1,j] = e[j - 1];
            }
            a[i - 1,i - 1] = s[i - 1];
            if (m < n) a[i - 1,i] = e[i - 1];
            for (i = 1; i <= n - 1; i++)
                for (j = i + 1; j <= n; j++)
                {
                    d = v[i - 1,j - 1];
                    v[i - 1,j - 1] = v[j - 1,i - 1];
                    v[j - 1,i - 1] = d;
                }
        }
        /// <summary>
        /// 奇异值分解的子程序
        /// </summary>
        static void SSS(double[] fg, double[] cs) {
            double r, d;
            if ((Math.Abs(fg[0]) + Math.Abs(fg[1])) == 0.0)
            {
                cs[0] = 1.0; cs[1] = 0.0; d = 0.0;
            }
            else
            {
                d = Math.Sqrt(fg[0] * fg[0] + fg[1] * fg[1]);
                if (Math.Abs(fg[0]) > Math.Abs(fg[1]))
                {
                    d = Math.Abs(d);
                    if (fg[0] < 0.0) d = -d;
                }
                if (Math.Abs(fg[1]) >= Math.Abs(fg[0]))
                {
                    d = Math.Abs(d);
                    if (fg[1] < 0.0) d = -d;
                }
                cs[0] = fg[0] / d; cs[1] = fg[1] / d;
            }
            r = 1.0;
            if (Math.Abs(fg[0]) > Math.Abs(fg[1])) r = cs[1];
            else
                if (cs[0] != 0.0) r = 1.0 / cs[0];
            fg[0] = d; fg[1] = r;
        }

        /// <summary>
        /// 利用奇异值分解求一般m×n阶实矩阵A的广义逆A+
        /// </summary>
        /// <param name="a">存放m×n的实矩阵A.返回时其对角线给出奇异值（以非递增次序排列），其余元素均为0.</param>
        /// <param name="eps">给定的精度要求</param>
        /// <param name="aa">返回A的广义逆A+</param>
        /// <param name="u">返回左奇异向量U</param>
        /// <param name="v">返回右奇异向量V</param>
        public static void GINV(double[,] a, double eps, out double[,] aa, out double[,] u, out double[,] v)
        {
            int m = a.GetLength(0);
            int n = a.GetLength(1);
            int ka = Math.Max(m, n) + 1;
            aa = new double[n, m];
            u = new double[m, m];
            v = new double[n, n];
            double[] fg = new double[2], cs = new double[2], s = new double[ka], e = new double[ka], w = new double[ka];
            MUAV(a, eps, out u, out v);
            int i, j, k, l;
            j = n;
            if (m < n) j = m;
            j = j - 1;
            k = 0;
            while ((k <= j) && (a[k,k] != 0.0)) k = k + 1;
            k = k - 1;
            for (i = 0; i <= n - 1; i++)
                for (j = 0; j <= m - 1; j++)
                {
                    aa[i,j] = 0.0;
                    for (l = 0; l <= k; l++)
                    {
                        aa[i,j] = aa[i,j] + v[l,i] * u[j,l] / a[l,l];
                    }
                }
        }
    }
}
