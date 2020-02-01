using System;
using System.Collections.Generic;
using System.Text;

namespace Heroius.XuAlgrithms
{
    /// <summary>
    /// 数学变换与滤波
    /// </summary>
    public static class Transformation
    {
        /// <summary>
        /// 根据函数f(x)在区间[0,2π]上的2n+1个等距点
        /// <para>xi=2π(i+0.5)/(2n+1), i=0,1,...,2n</para>
        /// <para>处的函数值 fi=f(xi), 求傅里叶（Fourier）级数</para>
        /// <para>f(x)=0.5a0+Σ(ak*Math.Cos(kx)+bk*Math.Sin(kx))</para>
        /// <para>的前2n+1个系数ak(k=0,1,...,n)和bk(k=0,1,...,n)的近似值</para>
        /// </summary>
        /// <param name="f">f[2n+1]: 存放区间[0,2π]内的2n+1个等距点处的函数值</param>
        /// <param name="n">等距点数为2n+1</param>
        /// <param name="a">a[n+1]: 返回傅里叶级数中的系数ak(k=0,1,...,n)</param>
        /// <param name="b">b[n+1}: 返回傅里叶级数中的系数bk(k=0,1,...,n)</param>
        public static void FOUR(double[] f, int n, out double[] a, out double[] b)
        {
            int i, j;
            double t, c, s, c1, s1, u1, u2, u0;
            a = new double[n + 1];
            b = new double[n + 1];
            t = 6.283185306 / (2.0 * n + 1.0);
            c = Math.Cos(t); s = Math.Sin(t);
            t = 2.0 / (2.0 * n + 1.0); c1 = 1.0; s1 = 0.0;
            for (i = 0; i <= n; i++)
            {
                u1 = 0.0; u2 = 0.0;
                for (j = 2 * n; j >= 1; j--)
                {
                    u0 = f[j] + 2.0 * c1 * u1 - u2;
                    u2 = u1; u1 = u0;
                }
                a[i] = t * (f[0] + u1 * c1 - u2);
                b[i] = t * u1 * s1;
                u0 = c * c1 - s * s1; s1 = c * s1 + s * c1; c1 = u0;
            }
            return;
        }

        /// <summary>
        /// 用 FFT 计算离散傅里叶 (Fourier) 变换 。 
        /// </summary>
        /// <param name="pr">pr[n]: 当l=0时，存放n个采样输入的实部，返回离散傅里叶变换的模；
        /// 当l=1时，存放傅里叶变换的n个实部，返回逆傅里叶变换的模</param>
        /// <param name="pi">pi[n]: 当l=0时，存放n个采样输入的虚部，返回离散傅里叶变换的幅角；
        /// 当l=1时，存放傅里叶变换的n个虚部，返回逆傅里叶变换的幅角。
        /// 其中幅角的单位为度</param>
        /// <param name="n">采样点数</param>
        /// <param name="k">满足 n=2^k</param>
        /// <param name="fr">fr[n]: 当l=0时，返回傅里叶变换的n个实部；当l=1时，返回逆傅里叶变换的n个实部</param>
        /// <param name="fi">fi[n]: 当l=0时，返回傅里叶变换的n个虚部；当l=1时，返回逆傅里叶变换的n个虚部</param>
        /// <param name="l">当l=0时，表示要求本函数计算傅里叶变换；
        /// 当l=1时，表示要求本函数计算逆傅里叶变换</param>
        /// <param name="il">当il=0时，表示不要求本函数计算傅里叶变换或逆傅里叶变换的模与幅角；
        /// 当il=1时，表示要求本函数计算傅里叶变换或逆傅里叶变换的模与幅角</param>
        public static void KFFT(ref double[] pr, ref double[] pi, int n, int k, out double[] fr, out double[] fi, int l, int il)
        {
            int it, m, _is, i, j, nv, l0;
            double p, q, s, vr, vi, poddr, poddi;
            fr = new double[n];
            fi = new double[n];
            for (it = 0; it <= n - 1; it++)
            {
                m = it; _is = 0;
                for (i = 0; i <= k - 1; i++)
                { j = m / 2; _is = 2 * _is + (m - 2 * j); m = j; }
                fr[it] = pr[_is]; fi[it] = pi[_is];
            }
            pr[0] = 1.0; pi[0] = 0.0;
            p = 6.283185306 / (1.0 * n);
            pr[1] = Math.Cos(p); pi[1] = -Math.Sin(p);
            if (l != 0) pi[1] = -pi[1];
            for (i = 2; i <= n - 1; i++)
            {
                p = pr[i - 1] * pr[1]; q = pi[i - 1] * pi[1];
                s = (pr[i - 1] + pi[i - 1]) * (pr[1] + pi[1]);
                pr[i] = p - q; pi[i] = s - p - q;
            }
            for (it = 0; it <= n - 2; it = it + 2)
            {
                vr = fr[it]; vi = fi[it];
                fr[it] = vr + fr[it + 1]; fi[it] = vi + fi[it + 1];
                fr[it + 1] = vr - fr[it + 1]; fi[it + 1] = vi - fi[it + 1];
            }
            m = n / 2; nv = 2;
            for (l0 = k - 2; l0 >= 0; l0--)
            {
                m = m / 2; nv = 2 * nv;
                for (it = 0; it <= (m - 1) * nv; it = it + nv)
                    for (j = 0; j <= (nv / 2) - 1; j++)
                    {
                        p = pr[m * j] * fr[it + j + nv / 2];
                        q = pi[m * j] * fi[it + j + nv / 2];
                        s = pr[m * j] + pi[m * j];
                        s = s * (fr[it + j + nv / 2] + fi[it + j + nv / 2]);
                        poddr = p - q; poddi = s - p - q;
                        fr[it + j + nv / 2] = fr[it + j] - poddr;
                        fi[it + j + nv / 2] = fi[it + j] - poddi;
                        fr[it + j] = fr[it + j] + poddr;
                        fi[it + j] = fi[it + j] + poddi;
                    }
            }
            if (l != 0)
                for (i = 0; i <= n - 1; i++)
                {
                    fr[i] = fr[i] / (1.0 * n);
                    fi[i] = fi[i] / (1.0 * n);
                }
            if (il != 0)
                for (i = 0; i <= n - 1; i++)
                {
                    pr[i] = Math.Sqrt(fr[i] * fr[i] + fi[i] * fi[i]);
                    if (Math.Abs(fr[i]) < 0.000001 * Math.Abs(fi[i]))
                    {
                        if ((fi[i] * fr[i]) > 0) pi[i] = 90.0;
                        else pi[i] = -90.0;
                    }
                    else
                        pi[i] = Math.Atan(fi[i] / fr[i]) * 360.0 / 6.283185306;
                }
            return;
        }

        /// <summary>
        /// 计算给定序列 的沃什 (Walsh) 变换序列
        /// </summary>
        /// <param name="p">p[n]: 存放长度n=2^k的给定输入序列</param>
        /// <param name="n">输入序列的长度</param>
        /// <param name="k">满足 n=2^k</param>
        /// <param name="x">x[n]: 返回输入序列pi(i=0,1,...,n-1)的沃什变换序列</param>
        public static void KFWT(double[] p, int n, int k, double[] x)
        {
            int m, l, it, ii, i, j, _is;
            double q;
            m = 1; l = n; it = 2;
            x[0] = 1; ii = n / 2; x[ii] = 2;
            for (i = 1; i <= k - 1; i++)
            {
                m = m + m; l = l / 2; it = it + it;
                for (j = 0; j <= m - 1; j++)
                    x[j * l + l / 2] = it + 1 - x[j * l];
            }
            for (i = 0; i <= n - 1; i++)
            { ii = (int)(x[i] - 1); x[i] = p[ii]; }
            l = 1;
            for (i = 1; i <= k; i++)
            {
                m = n / (2 * l) - 1;
                for (j = 0; j <= m; j++)
                {
                    it = 2 * l * j;
                    for (_is = 0; _is <= l - 1; _is++)
                    {
                        q = x[it + _is] + x[it + _is + l];
                        x[it + _is + l] = x[it + _is] - x[it + _is + l];
                        x[it + _is] = q;
                    }
                }
                l = 2 * l;
            }
            return;
        }

        /// <summary>
        /// 用五点三次平滑公式对等距点上的观测数据进行平滑
        /// </summary>
        /// <param name="n">给定等距观测点数。要求 n≥5 </param>
        /// <param name="y">y[n]: 存放n个等距观测点上的观测数据</param>
        /// <param name="yy">yy[n]: 返回n个等距观测点上的平滑结果</param>
        public static void KSPT(int n, double[] y, out double[] yy)
        {
            int i;
            yy = new double[n];
            if (n < 5)
            { for (i = 0; i <= n - 1; i++) yy[i] = y[i]; }
            else
            {
                yy[0] = 69.0 * y[0] + 4.0 * y[1] - 6.0 * y[2] + 4.0 * y[3] - y[4];
                yy[0] = yy[0] / 70.0;
                yy[1] = 2.0 * y[0] + 27.0 * y[1] + 12.0 * y[2] - 8.0 * y[3];
                yy[1] = (yy[1] + 2.0 * y[4]) / 35.0;
                for (i = 2; i <= n - 3; i++)
                {
                    yy[i] = -3.0 * y[i - 2] + 12.0 * y[i - 1] + 17.0 * y[i];
                    yy[i] = (yy[i] + 12.0 * y[i + 1] - 3.0 * y[i + 2]) / 35.0;
                }
                yy[n - 2] = 2.0 * y[n - 5] - 8.0 * y[n - 4] + 12.0 * y[n - 3];
                yy[n - 2] = (yy[n - 2] + 27.0 * y[n - 2] + 2.0 * y[n - 1]) / 35.0;
                yy[n - 1] = -y[n - 5] + 4.0 * y[n - 4] - 6.0 * y[n - 3];
                yy[n - 1] = (yy[n - 1] + 4.0 * y[n - 2] + 69.0 * y[n - 1]) / 70.0;
            }
            return;
        }

        /// <summary>
        ///  对实践离散点上的采样数据进行卡尔曼 (Kalman) 滤波
        /// </summary>
        /// <param name="n">动态系统的维数</param>
        /// <param name="m">观测系统的维数</param>
        /// <param name="k">观测序列的长度</param>
        /// <param name="f">f[n,n]: 系统状态转移矩阵</param>
        /// <param name="q">q[n,n]: 模型噪声Wk的协方差阵</param>
        /// <param name="r">r[m,m]: 观测噪声Vk的协方差阵</param>
        /// <param name="h">h[m,n]: 观测矩阵</param>
        /// <param name="y">y[k,m]: 观测向量序列，其中y(i,j)(i=0,1,...k-1;j=0,1,...,m-1)表示第i时刻的观测向量的第j个分量</param>
        /// <param name="x">x[k,n]: 其中x(0,j)(j=0,1,...,n-1)存放给定的初值，
        /// 其余各行返回状态向量估值序列。
        /// X(i,j)(i=0,1,...,k-1;j=0,1,...,n-1)表示第i时刻的状态向量估值的第j个分量</param>
        /// <param name="p">p[n,n]: 存放初值P0，返回最后时刻的估计误差协方差阵</param>
        /// <param name="g">g[n,m]: 返回最后时刻的稳定增益矩阵</param>
        /// <returns></returns>
        public static int LMAN(int n, int m, int k, double[,] f, double[,] q, double[,] r, double[,] h, double[,] y, ref double[,] x, ref double[,] p, out double[,] g)
        {
            int i, j, kk, ii, l, jj, js = 0;
            double[] pp = Utility.C.Convert(p);
            double[] ff = Utility.C.Convert(f);
            double[] hh = Utility.C.Convert(h);
            double[] rr = Utility.C.Convert(r);
            double[] qq = Utility.C.Convert(q);
            double[] xx = Utility.C.Convert(x);
            double[] yy = Utility.C.Convert(y);
            double[] gg = new double[n * m];
            double[] e = new double[m * m], a, b;
            l = m;
            if (l < n) l = n;
            a = new double[l * l];
            b = new double[l * l];
            for (i = 0; i <= n - 1; i++)
                for (j = 0; j <= n - 1; j++)
                {
                    ii = i * l + j; a[ii] = 0.0;
                    for (kk = 0; kk <= n - 1; kk++)
                        a[ii] = a[ii] + pp[i * n + kk] * ff[j * n + kk];
                }
            for (i = 0; i <= n - 1; i++)
                for (j = 0; j <= n - 1; j++)
                {
                    ii = i * n + j; pp[ii] = qq[ii];
                    for (kk = 0; kk <= n - 1; kk++)
                        pp[ii] = pp[ii] + ff[i * n + kk] * a[kk * l + j];
                }
            for (ii = 2; ii <= k; ii++)
            {
                for (i = 0; i <= n - 1; i++)
                    for (j = 0; j <= m - 1; j++)
                    {
                        jj = i * l + j; a[jj] = 0.0;
                        for (kk = 0; kk <= n - 1; kk++)
                            a[jj] = a[jj] + pp[i * n + kk] * hh[j * n + kk];
                    }
                for (i = 0; i <= m - 1; i++)
                    for (j = 0; j <= m - 1; j++)
                    {
                        jj = i * m + j; e[jj] = rr[jj];
                        for (kk = 0; kk <= n - 1; kk++)
                            e[jj] = e[jj] + hh[i * n + kk] * a[kk * l + j];
                    }

                var ee = Utility.C.Convert(e, m, m);
                js = Matrix.RINV(ref ee, m);
                e = Utility.C.Convert(ee);

                if (js == 0)
                {
                    p = Utility.C.Convert(pp, n, n);
                    g = Utility.C.Convert(gg, n, m);
                    x = Utility.C.Convert(xx, k, n);
                    return (js);
                }
                for (i = 0; i <= n - 1; i++)
                    for (j = 0; j <= m - 1; j++)
                    {
                        jj = i * m + j; gg[jj] = 0.0;
                        for (kk = 0; kk <= m - 1; kk++)
                            gg[jj] = gg[jj] + a[i * l + kk] * e[j * m + kk];
                    }
                for (i = 0; i <= n - 1; i++)
                {
                    jj = (ii - 1) * n + i; xx[jj] = 0.0;
                    for (j = 0; j <= n - 1; j++)
                        xx[jj] = xx[jj] + ff[i * n + j] * xx[(ii - 2) * n + j];
                }
                for (i = 0; i <= m - 1; i++)
                {
                    jj = i * l; b[jj] = yy[(ii - 1) * m + i];
                    for (j = 0; j <= n - 1; j++)
                        b[jj] = b[jj] - hh[i * n + j] * xx[(ii - 1) * n + j];
                }
                for (i = 0; i <= n - 1; i++)
                {
                    jj = (ii - 1) * n + i;
                    for (j = 0; j <= m - 1; j++)
                        xx[jj] = xx[jj] + gg[i * m + j] * b[j * l];
                }
                if (ii < k)
                {
                    for (i = 0; i <= n - 1; i++)
                        for (j = 0; j <= n - 1; j++)
                        {
                            jj = i * l + j; a[jj] = 0.0;
                            for (kk = 0; kk <= m - 1; kk++)
                                a[jj] = a[jj] - gg[i * m + kk] * hh[kk * n + j];
                            if (i == j) a[jj] = 1.0 + a[jj];
                        }
                    for (i = 0; i <= n - 1; i++)
                        for (j = 0; j <= n - 1; j++)
                        {
                            jj = i * l + j; b[jj] = 0.0;
                            for (kk = 0; kk <= n - 1; kk++)
                                b[jj] = b[jj] + a[i * l + kk] * pp[kk * n + j];
                        }
                    for (i = 0; i <= n - 1; i++)
                        for (j = 0; j <= n - 1; j++)
                        {
                            jj = i * l + j; a[jj] = 0.0;
                            for (kk = 0; kk <= n - 1; kk++)
                                a[jj] = a[jj] + b[i * l + kk] * ff[j * n + kk];
                        }
                    for (i = 0; i <= n - 1; i++)
                        for (j = 0; j <= n - 1; j++)
                        {
                            jj = i * n + j; pp[jj] = qq[jj];
                            for (kk = 0; kk <= n - 1; kk++)
                                pp[jj] = pp[jj] + ff[i * n + kk] * a[j * l + kk];
                        }
                }
            }
            p = Utility.C.Convert(pp, n, n);
            g = Utility.C.Convert(gg, n, m);
            x = Utility.C.Convert(xx, k, n);
            return (js);
        }

        /// <summary>
        /// 对等间隔的量测数据进行滤波估值
        /// </summary>
        /// <param name="n">量测数据的点数</param>
        /// <param name="x">x[n]: n个等间隔点上的量测值</param>
        /// <param name="t">采样周期</param>
        /// <param name="a">滤波器结构参数α</param>
        /// <param name="b">滤波器结构参数β</param>
        /// <param name="c">滤波器结构参数γ</param>
        /// <param name="y">y[n]: 返回n个等间隔点上的滤波估值</param>
        public static void KABG(int n, double[] x, double t, double a, double b, double c, out double[] y)
        {
            int i;
            double s1, ss, v1, vv, a1, aa;
            aa = 0.0; vv = 0.0; ss = 0.0;
            y = new double[n];
            for (i = 0; i <= n - 1; i++)
            {
                s1 = ss + t * vv + t * t * aa / 2.0;
                v1 = vv + t * aa; a1 = aa;
                ss = s1 + a * (x[i] - s1); y[i] = ss;
                vv = v1 + b * (x[i] - s1);
                aa = a1 + 2.0 * c * (x[i] - s1) / (t * t);
            }
            return;
        }
    }
}
