using System;
using System.Collections.Generic;
using System.Text;

namespace Heroius.XuAlgrithms
{
    public static partial class Algrithms
    {
        /// <summary>
        /// 用变形QR方法计算实对称三角矩阵的全部特征值与相应的向量
        /// </summary>
        /// <param name="b">存放n阶实对称三角阵的主对角线上的元素</param>
        /// <param name="c">前n-1个元素存放实对称三角阵的次对角线上的元素，返回时存放特征值</param>
        /// <param name="eps">控制精度要求</param>
        /// <param name="Max">允许的最大迭代次数</param>
        /// <returns>存放实对称三角矩阵的特征向量组</returns>
        public static double[,] SSTQ(double[] b, double[] c, double eps, int Max)
        {
            int n = b.Length;
            if (c.Length != n)
            {
                throw new Exception("次对角线数组长度不足");
            }
            double[,] q = new double[n, n];
            for (int ii = 0; ii < n; ii++)
            {
                q[ii, ii] = 1;
            }
            int i, j, k, m, it;
            double d, f, h, g, p, r, e, s;
            c[n - 1] = 0.0; d = 0.0; f = 0.0;
            for (j = 0; j <= n - 1; j++)
            {
                it = 0;
                h = eps * (Math.Abs(b[j]) + Math.Abs(c[j]));
                if (h > d) d = h;
                m = j;
                while ((m <= n - 1) && (Math.Abs(c[m]) > d)) m = m + 1;
                if (m != j)
                {
                    do
                    {
                        if (it == Max)
                        {
                            throw new Exception("程序工作失败！");
                        }
                        it = it + 1;
                        g = b[j];
                        p = (b[j + 1] - g) / (2.0 * c[j]);
                        r = Math.Sqrt(p * p + 1.0);
                        if (p >= 0.0) b[j] = c[j] / (p + r);
                        else b[j] = c[j] / (p - r);
                        h = g - b[j];
                        for (i = j + 1; i <= n - 1; i++) b[i] = b[i] - h;
                        f = f + h; p = b[m]; e = 1.0; s = 0.0;
                        for (i = m - 1; i >= j; i--)
                        {
                            g = e * c[i]; h = e * p;
                            if (Math.Abs(p) >= Math.Abs(c[i]))
                            {
                                e = c[i] / p; r = Math.Sqrt(e * e + 1.0);
                                c[i + 1] = s * p * r; s = e / r; e = 1.0 / r;
                            }
                            else
                            {
                                e = p / c[i]; r = Math.Sqrt(e * e + 1.0);
                                c[i + 1] = s * c[i] * r;
                                s = 1.0 / r; e = e / r;
                            }
                            p = e * b[i] - s * g;
                            b[i + 1] = h + s * (e * g + s * b[i]);
                            for (k = 0; k <= n - 1; k++)
                            {
                                h = q[k, i + 1]; q[k, i + 1] = s * q[k, i] + e * h;
                                q[k, i] = e * q[k, i] - s * h;
                            }
                        }
                        c[j] = s * p; b[j] = e * p;
                    }
                    while (Math.Abs(c[j]) > d);
                }
                b[j] = b[j] + f;
            }
            for (i = 0; i <= n - 1; i++)
            {
                k = i; p = b[i];
                if (i + 1 <= n - 1)
                {
                    j = i + 1;
                    while ((j <= n - 1) && (b[j] <= p))
                    { k = j; p = b[j]; j = j + 1; }
                }
                if (k != i)
                {
                    b[k] = b[i]; b[i] = p;
                    for (j = 0; j <= n - 1; j++)
                    {
                        p = q[j, i]; q[j, i] = q[j, k]; q[j, k] = p;
                    }
                }
            }
            return q;
        }


        public static double[,] STRQ(double[,] a)
        {
            int n = a.GetLength(0);
            if (a.GetLength(1) != n)
            {
                throw new Exception("非矩阵");
            }
            double[,] q = new double[n, n];
            double[] b = new double[n];
            double[] c = new double[n];

            int i, j, k;
            double h, f, g, h2;
            for (i = 0; i <= n - 1; i++)
                for (j = 0; j <= n - 1; j++) q[i,j] = a[i,j];
            for (i = n - 1; i >= 1; i--)
            {
                h = 0.0;
                if (i > 1)
                    for (k = 0; k <= i - 1; k++) h = h + q[i,k] * q[i,k];
                if (h + 1.0 == 1.0)
                {
                    c[i] = 0.0;
                    if (i == 1) c[i] = q[i,i - 1];
                    b[i] = 0.0;
                }
                else
                {
                    c[i] = Math.Sqrt(h);
                    if (q[i,i - 1] > 0.0) c[i] = -c[i];
                    h = h - q[i,i - 1] * c[i];
                    q[i,i - 1] = q[i,i - 1] - c[i];
                    f = 0.0;
                    for (j = 0; j <= i - 1; j++)
                    {
                        q[j,i] = q[i,j] / h;
                        g = 0.0;
                        for (k = 0; k <= j; k++) g = g + q[j,k] * q[i,k];
                        if (j + 1 <= i - 1)
                            for (k = j + 1; k <= i - 1; k++)
                                g = g + q[k,j] * q[i,k];
                        c[j] = g / h;
                        f = f + g * q[j,i];
                    }
                    h2 = f / (h + h);
                    for (j = 0; j <= i - 1; j++)
                    {
                        f = q[i,j];
                        g = c[j] - h2 * f;
                        c[j] = g;
                        for (k = 0; k <= j; k++)
                        {
                            q[j,k] = q[j,k] - f * c[k] - g * q[i,k];
                        }
                    }
                    b[i] = h;
                }
            }
            for (i = 0; i <= n - 2; i++) c[i] = c[i + 1];
            c[n - 1] = 0.0;
            b[0] = 0.0;
            for (i = 0; i <= n - 1; i++)
            {
                if ((b[i] != 0.0) && (i - 1 >= 0))
                    for (j = 0; j <= i - 1; j++)
                    {
                        g = 0.0;
                        for (k = 0; k <= i - 1; k++) g = g + q[i,k] * q[k,j];
                        for (k = 0; k <= i - 1; k++)
                        {
                            q[k,j] = q[k,j] - g * q[k,i];
                        }
                    }
                b[i] = q[i,i]; q[i,i] = 1.0;
                if (i - 1 >= 0)
                    for (j = 0; j <= i - 1; j++)
                    {
                        q[i,j] = 0.0; q[j,i] = 0.0;
                    }
            }
            return q;
        }
    }
}
