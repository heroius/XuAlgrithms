using System;
using System.Collections.Generic;
using System.Text;

namespace Heroius.XuAlgrithms
{
    /// <summary>
    /// 矩阵特征值与特征向量的计算
    /// </summary>
    public static partial class MatrixEigen
    {
        /// <summary>
        /// 用变形QR方法计算实对称三角矩阵的全部特征值与相应的向量
        /// </summary>
        /// <param name="b">存放n阶实对称三角阵的主对角线上的元素，返回时存放全部特征值</param>
        /// <param name="c">前n-1个元素存放实对称三角阵的次对角线上的元素</param>
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

        /// <summary>
        /// 用豪斯荷尔德(Householder)变换将n阶实对称矩阵约化为对称三对角阵。
        /// </summary>
        /// <param name="a">存放n阶实对称矩阵</param>
        /// <param name="q">返回豪斯荷尔德变换的乘积矩阵Q。</param>
        /// <param name="b">返回对称三角阵中的主对角线元素。</param>
        /// <param name="c">前n-1个元素返回对称三角阵中的次对角线元素。</param>
        public static void STRQ(double[,] a, out double[,] q, out double[] b, out double[] c)
        {
            int n = a.GetLength(0);
            if (a.GetLength(1) != n)
            {
                throw new Exception("非矩阵");
            }
            q = new double[n, n];
            b = new double[n];
            c = new double[n];

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
        }

        /// <summary>
        /// 用初等相似变换将一般实矩阵约化为上H矩阵，即赫申伯格(Hessenberg)矩阵。
        /// </summary>
        /// <param name="a">存放一般实矩阵A，返回上H矩阵。</param>
        public static void HHBG(double[,] a)
        {
            int n = a.GetLength(0);
            if (a.GetLength(1) != n)
            {
                throw new Exception("非矩阵");
            }
            int i = 0, j, k;
            double d, t;
            for (k = 1; k <= n - 2; k++)
            {
                d = 0.0;
                for (j = k; j <= n - 1; j++)
                {
                    t = a[j,k - 1];
                    if (Math.Abs(t) > Math.Abs(d)) { d = t; i = j; }
                }
                if (Math.Abs(d) + 1.0 != 1.0)
                {
                    if (i != k)
                    {
                        for (j = k - 1; j <= n - 1; j++)
                        {
                            t = a[i,j]; a[i,j] = a[k,j]; a[k,j] = t;
                        }
                        for (j = 0; j <= n - 1; j++)
                        {
                            t = a[j,i]; a[j,i] = a[j,k]; a[j,k] = t;
                        }
                    }
                    for (i = k + 1; i <= n - 1; i++)
                    {
                        t = a[i,k - 1] / d; a[i,k - 1] = 0.0;
                        for (j = k; j <= n - 1; j++)
                            a[i,j] = a[i,j] - t * a[k,j];
                        for (j = 0; j <= n - 1; j++)
                            a[j,k] = a[j,k] + t * a[j,i];
                    }
                }
            }
        }

        /// <summary>
        /// 用带原点位移的双重步QR方法计算实上H矩阵的全部特征值。
        /// </summary>
        /// <param name="a">存放上H矩阵A</param>
        /// <param name="eps">控制精度要求</param>
        /// <param name="Max">控制最大迭代次数</param>
        /// <param name="u">返回n个特征值的实部</param>
        /// <param name="v">返回n个特征值的虚部</param>
        public static void HHQR(double[,] a, double eps, int Max, out double[] u, out double[] v)
        {
            int n = a.GetLength(0);
            if (a.GetLength(1) != n)
            {
                throw new Exception("非矩阵");
            }
            u = new double[n];
            v = new double[n];

            int m, it, i, j, k, l;
            double b, c, w, g, xy, p, q, r, x, s, e, f, z, y;
            it = 0; m = n;
            while (m != 0)
            {
                l = m - 1;
                while ((l > 0) && (Math.Abs(a[l,l - 1]) > eps *
                      (Math.Abs(a[l - 1,l - 1]) + Math.Abs(a[l,l])))) l = l - 1;
                if (l == m - 1)
                {
                    u[m - 1] = a[m - 1,m - 1]; v[m - 1] = 0.0;
                    m = m - 1; it = 0;
                }
                else if (l == m - 2)
                {
                    b = -(a[m - 1,m - 1] + a[m - 2,m - 2]);
                    c = a[m - 1,m - 1] * a[m - 2,m - 2] - a[m - 1,m - 2] * a[m - 2,m - 1];
                    w = b * b - 4.0 * c;
                    y = Math.Sqrt(Math.Abs(w));
                    if (w > 0.0)
                    {
                        xy = 1.0;
                        if (b < 0.0) xy = -1.0;
                        u[m - 1] = (-b - xy * y) / 2.0;
                        u[m - 2] = c / u[m - 1];
                        v[m - 1] = 0.0; v[m - 2] = 0.0;
                    }
                    else
                    {
                        u[m - 1] = -b / 2.0; u[m - 2] = u[m - 1];
                        v[m - 1] = y / 2.0; v[m - 2] = -v[m - 1];
                    }
                    m = m - 2; it = 0;
                }
                else
                {
                    if (it >= Max)
                    {
                        throw new Exception("程序工作失败！");
                    }
                    it = it + 1;
                    for (j = l + 2; j <= m - 1; j++) a[j,j - 2] = 0.0;
                    for (j = l + 3; j <= m - 1; j++) a[j,j - 3] = 0.0;
                    for (k = l; k <= m - 2; k++)
                    {
                        if (k != l)
                        {
                            p = a[k,k - 1]; q = a[k + 1,k - 1];
                            r = 0.0;
                            if (k != m - 2) r = a[k + 2,k - 1];
                        }
                        else
                        {
                            x = a[m - 1,m - 1] + a[m - 2,m - 2];
                            y = a[m - 2,m - 2] * a[m - 1,m - 1] - a[m - 2,m - 1] * a[m - 1,m - 2];
                            p = a[l,l] * (a[l,l] - x) + a[l,l + 1] * a[l + 1,l] + y;
                            q = a[l + 1,l] * (a[l,l] + a[l + 1,l + 1] - x);
                            r = a[l + 1,l] * a[l + 2,l + 1];
                        }
                        if ((Math.Abs(p) + Math.Abs(q) + Math.Abs(r)) != 0.0)
                        {
                            xy = 1.0;
                            if (p < 0.0) xy = -1.0;
                            s = xy * Math.Sqrt(p * p + q * q + r * r);
                            if (k != l) a[k,k - 1] = -s;
                            e = -q / s; f = -r / s; x = -p / s;
                            y = -x - f * r / (p + s);
                            g = e * r / (p + s);
                            z = -x - e * q / (p + s);
                            for (j = k; j <= m - 1; j++)
                            {
                                p = x * a[k,j] + e * a[k + 1,j];
                                q = e * a[k,j] + y * a[k + 1,j];
                                r = f * a[k,j] + g * a[k + 1,j];
                                if (k != m - 2)
                                {
                                    p = p + f * a[k + 2,j];
                                    q = q + g * a[k + 2,j];
                                    r = r + z * a[k + 2,j]; a[k + 2,j] = r;
                                }
                                a[k + 1,j] = q; a[k,j] = p;
                            }
                            j = k + 3;
                            if (j >= m - 1) j = m - 1;
                            for (i = l; i <= j; i++)
                            {
                                p = x * a[i,k] + e * a[i,k + 1];
                                q = e * a[i,k] + y * a[i,k + 1];
                                r = f * a[i,k] + g * a[i,k + 1];
                                if (k != m - 2)
                                {
                                    p = p + f * a[i,k + 2];
                                    q = q + g * a[i,k + 2];
                                    r = r + z * a[i,k + 2]; a[i,k + 2] = r;
                                }
                                a[i,k + 1] = q; a[i,k] = p;
                            }
                        }
                    }
                }
            }
        }

        /// <summary>
        /// 用雅可比(Jacobi)方法求实对称矩阵的全部特征值与相应的特征向量。
        /// </summary>
        /// <param name="a">存放n阶实对称矩阵，返回时对角线上存放n个特征值</param>
        /// <param name="eps">控制精度要求</param>
        /// <param name="Max">控制最大迭代次数</param>
        /// <returns></returns>
        public static double[,] JCBI(double[,] a, double eps, int Max)
        {
            int n = a.GetLength(0);
            if (a.GetLength(1) != n)
            {
                throw new Exception("非矩阵");
            }

            double[,] v = new double[n, n];
            int i, j, p = 0, q = 0, l;
            double fm, cn, sn, omega, x, y, d;
            l = 1;
            for (i = 0; i <= n - 1; i++)
            {
                v[i,i] = 1.0;
                for (j = 0; j <= n - 1; j++)
                    if (i != j) v[i,j] = 0.0;
            }
            while (1 == 1)
            {
                fm = 0.0;
                for (i = 1; i <= n - 1; i++)
                    for (j = 0; j <= i - 1; j++)
                    {
                        d = Math.Abs(a[i,j]);
                        if ((i != j) && (d > fm))
                        {
                            fm = d; p = i; q = j;
                        }
                    }
                if (fm < eps) return v;
                if (l > Max)
                {
                    throw new Exception("程序工作失败！");
                }
                l = l + 1;
                x = -a[p,q]; y = (a[q,q] - a[p,p]) / 2.0;
                omega = x / Math.Sqrt(x * x + y * y);
                if (y < 0.0) omega = -omega;
                sn = 1.0 + Math.Sqrt(1.0 - omega * omega);
                sn = omega / Math.Sqrt(2.0 * sn);
                cn = Math.Sqrt(1.0 - sn * sn);
                fm = a[p,p];
                a[p,p] = fm * cn * cn + a[q,q] * sn * sn + a[p,q] * omega;
                a[q,q] = fm * sn * sn + a[q,q] * cn * cn - a[p,q] * omega;
                a[p,q] = 0.0; a[q,p] = 0.0;
                for (j = 0; j <= n - 1; j++)
                    if ((j != p) && (j != q))
                    {
                        fm = a[p,j];
                        a[p,j] = fm * cn + a[q,j] * sn;
                        a[q,j] = -fm * sn + a[q,j] * cn;
                    }
                for (i = 0; i <= n - 1; i++)
                    if ((i != p) && (i != q))
                    {
                        fm = a[i,p];
                        a[i,p] = fm * cn + a[i,q] * sn;
                        a[i,q] = -fm * sn + a[i,q] * cn;
                    }
                for (i = 0; i <= n - 1; i++)
                {
                    fm = v[i,p];
                    v[i,p] = fm * cn + v[i,q] * sn;
                    v[i,q] = -fm * sn + v[i,q] * cn;
                }
            }
        }

        /// <summary>
        /// 用雅可比(Jacobi)过关法求实对称矩阵的全部特征值与相应的特征向量
        /// </summary>
        /// <param name="a">存放n阶实对称矩阵，返回时对角线上存放n个特征值</param>
        /// <param name="eps">控制精度要求</param>
        /// <returns>返回特征向量，其中第j列与第j个特征值对应的特征向量</returns>
        public static double[,] JCBJ(double[,] a, double eps)
        {
            int n = a.GetLength(0);
            if (a.GetLength(1) != n)
            {
                throw new Exception("非矩阵");
            }
            double[,] v = new double[n, n];

            int i, j, p, q;
            double ff, fm, cn, sn, omega, x, y, d;
            for (i = 0; i <= n - 1; i++)
            {
                v[i,i] = 1.0;
                for (j = 0; j <= n - 1; j++)
                    if (i != j) v[i,j] = 0.0;
            }
            ff = 0.0;
            for (i = 1; i <= n - 1; i++)
                for (j = 0; j <= i - 1; j++)
                { d = a[i,j]; ff = ff + d * d; }
            ff = Math.Sqrt(2.0 * ff);
            loop0:
            ff = ff / (1.0 * n);
            loop1:
            for (i = 1; i <= n - 1; i++)
                for (j = 0; j <= i - 1; j++)
                {
                    d = Math.Abs(a[i,j]);
                    if (d > ff)
                    {
                        p = i; q = j;
                        goto loop;
                    }
                }
            if (ff < eps) return v;
            goto loop0;
            loop:
            x = -a[p,q]; y = (a[q,q] - a[p,p]) / 2.0;
            omega = x / Math.Sqrt(x * x + y * y);
            if (y < 0.0) omega = -omega;
            sn = 1.0 + Math.Sqrt(1.0 - omega * omega);
            sn = omega / Math.Sqrt(2.0 * sn);
            cn = Math.Sqrt(1.0 - sn * sn);
            fm = a[p,p];
            a[p,p] = fm * cn * cn + a[q,q] * sn * sn + a[p,q] * omega;
            a[q,q] = fm * sn * sn + a[q,q] * cn * cn - a[p,q] * omega;
            a[p,q] = 0.0; a[q,p] = 0.0;
            for (j = 0; j <= n - 1; j++)
                if ((j != p) && (j != q))
                {
                    fm = a[p,j];
                    a[p,j] = fm * cn + a[q,j] * sn;
                    a[q,j] = -fm * sn + a[q,j] * cn;
                }
            for (i = 0; i <= n - 1; i++)
                if ((i != p) && (i != q))
                {
                    fm = a[i,p];
                    a[i,p] = fm * cn + a[i,q] * sn;
                    a[i,q] = -fm * sn + a[i,q] * cn;
                }
            for (i = 0; i <= n - 1; i++)
            {
                fm = v[i,p];
                v[i,p] = fm * cn + v[i,q] * sn;
                v[i,q] = -fm * sn + v[i,q] * cn;
            }
            goto loop1;
        }
    }
}
