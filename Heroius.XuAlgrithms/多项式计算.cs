using System;
using System.Collections.Generic;
using System.Text;

namespace Heroius.XuAlgrithms
{
    /// <summary>
    /// 多项式计算
    /// </summary>
    public static class Polynomial
    {
        /// <summary>
        /// 计算多项式p(x)=Σax在指定点x处的函数值
        /// </summary>
        /// <param name="a">存放n-l次多项式的n个系数</param>
        /// <param name="n">多项式的项数</param>
        /// <param name="x">指定的自变扯值</param>
        /// <returns>返回多项式值p(x)</returns>
        public static double PLYV(double[] a, int n, double x)
        {
            double u = a[n - 1];
            for (int i = n - 2; i >= 0; i--)
            {
                u = u * x + a[i];
            }
            return u;
        }

        /// <summary>
        /// 利用系数预处理法对多项式p(x)=Σax同时对多个x值进行求值，其中n=2^k (k≥1) 。
        /// </summary>
        /// <param name="a">存放n-l次多项式的n个系数，返回时其值将被改变</param>
        /// <param name="n">多项式的项数</param>
        /// <param name="x">存放给定的m 个自变量值</param>
        /// <param name="m">给定自变量的个数</param>
        /// <param name="p">返回时存放与给定m个自变量值对应的多项式值</param>
        public static void PLYS(ref double[] a, int n, double[] x, int m, out double[] p)
        {
            int i, j, mm, nn, ll, t, s, kk, k;
            double[] b = new double[2 * n];
            double y = a[n - 1], z;
            p = new double[m];
            for (i = 0; i < n; i++)
            {
                b[i] = a[i] / y;
            }
            k = Convert.ToInt32(Math.Round(Math.Log(n - 0.5) / Math.Log(2) + 1));
            nn = 1;
            for (i = 0; i < k; i++)
            {
                nn *= 2;
            }
            for (i = n; i < nn - 1; i++)
            {
                b[i] = 0;
            }
            b[nn - 1] = 1;
            t = nn;
            s = 1;
            for (i = 1; i < k; i++)
            {
                t /= 2;
                mm = -t;
                for (j = 0; j <= s; j++)
                {
                    mm = mm + t + t;
                    b[mm - 1] = b[mm - 1] - 1;
                    for (kk = 2; kk <= t; kk++)
                    {
                        b[mm - kk] = b[mm - kk] - b[mm - 1] * b[mm + t - kk];
                    }
                }
                s += s;
            }
            for (kk = 1; kk <= m; kk++)
            {
                for (i = 0; i <= (nn - 2) / 2; i++)
                {
                    a[i] = x[kk - 1] + b[2 * i];
                }
                mm = 1;
                z = x[kk - 1];
                for (i = 1; i <= k - 1; i++)
                {
                    mm += mm;
                    ll = mm + mm;
                    z *= z;
                    for (j = 0; j < nn; j += ll)
                    {
                        a[j / 2] = a[j / 2] + a[(j + mm) / 2] * (z + b[j + mm - 1]);
                    }
                }
                z = z * z / x[kk - 1];
                if (nn != n) a[0] = a[0] - z;
                p[kk - 1] = a[0] * y;
            }
        }

        /// <summary>
        /// 计算二维多项式p(x,y)=ΣΣaxy在给定点(x,y)处的函数值
        /// </summary>
        /// <param name="a">存放二维多项式的系数</param>
        /// <param name="m">自变量x的最高次数为m-l</param>
        /// <param name="n">自变量y的最高次数为n-l</param>
        /// <param name="x">给定的一对自变址值</param>
        /// <param name="y">给定的一对自变址值</param>
        /// <returns>返回函数值p(x,y)</returns>
        public static double BPLY(double[,] a, int m, int n, double x, double y)
        {
            double u = 0, s, xx = 1;
            double[] aa = Utility.C.Convert(a);
            for (int i = 0; i < m; i++)
            {
                s = aa[i * n + n - 1] * xx;
                for (int j = n-2; j >= 0; j--)
                {
                    s = s * y + aa[i * n + j] * xx;
                }
                u = u + s;
                xx = xx * x;
            }
            return u;
        }

        /// <summary>
        /// 计算复系数多项式 p(z)=Σaz 在给定复数z 时的函数值。
        /// </summary>
        /// <param name="ar">存放多项式系数的实部</param>
        /// <param name="ai">存放多项式系数的虚部</param>
        /// <param name="n">多项式的项数，其最高次数为n-l</param>
        /// <param name="x">给定复数z的实部</param>
        /// <param name="y">给定复数z的虚部</param>
        /// <param name="u">返回多项式值p(z)的实部</param>
        /// <param name="v">返回多项式值p(z)的虚部</param>
        public static void CPLY(double[] ar, double[] ai, int n, double x, double y, out double u, out double v)
        {
            int i;
            double p, q, s = ar[n-1], t=ai[n-1];
            for (i = n - 2; i >= 0; i--)
            {
                Complex.CMUL(s, t, x, y, out p, out q);
                s = p + ar[i];
                t = q + ai[i];
            }
            u = s;
            v = t;
        }

        /// <summary>
        /// 求两个多项式 P(x)=Σpx Q(x)=Σqx 的乘积多项式 S(x)=P(x)Q(x)=Σsx
        /// </summary>
        /// <param name="p">存放多项式P(x)的系数</param>
        /// <param name="m">多项式P(x)的项数，其最高次数为m-1</param>
        /// <param name="q">存放多项式Q(x)的系数</param>
        /// <param name="n">多项式Q(x)的项数，其最高次数为n-1</param>
        /// <param name="s">返回乘积多项式的系数</param>
        /// <param name="k">返回多项式乘积S(x)的项数，其最高次数为k-1，k=m+n-1</param>
        public static void PMUL(double[] p, int m, double[] q, int n, out double[] s, out int k)
        {
            k = m + n - 1;
            s = new double[k];
            for (int i = 0; i <=k-1; i++)
            {
                s[i] = 0;
            }
        }

        /// <summary>
        /// 求两个复系数多项式 P(z)=Σpz Q(z)=Σqz 的乘积多项式 S(z)=P(z)Q(z)=Σsz
        /// </summary>
        /// <param name="pr">存放多项式P(z)系数的实部</param>
        /// <param name="pi">存放多项式P(z)系数的虚部</param>
        /// <param name="m">多项式P(z)的项数，其最高次数为m-1</param>
        /// <param name="qr">存放多项式Q(z)系数的实部</param>
        /// <param name="qi">存放多项式Q(z)系数的虚部</param>
        /// <param name="n">多项式Q(z)的项数，其最高次数为n-1</param>
        /// <param name="sr">返回乘积多项式系数的实部</param>
        /// <param name="si">返回乘积多项式系数的虚部</param>
        /// <param name="k">返回多项式乘积S(z)的项数，其最高次数为k-1，k=m+n-1</param>
        public static void CPML(double[] pr, double[] pi, int m, double[] qr, double[] qi, int n, out double[] sr, out double[] si, out int k)
        {
            k = m + n - 1;
            sr = new double[k];
            si = new double[k];
            double a, b, c, d, u, v;
            for (int i = 0; i < k; i++)
            {
                sr[i] = 0; si[i] = 0;
            }
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    a = pr[i];
                    b = pi[i];
                    c = qr[j];
                    d = qi[j];
                    Complex.CMUL(a, b, c, d, out u, out v);
                    sr[i + j] += u;
                    si[i + j] += v;
                }
            }
        }

        /// <summary>
        /// 求多项式P(x)=Σpx被多项式Q(x)=Σqx除得的商多项式S(x)和余多项式R(x)。
        /// </summary>
        /// <param name="p">存放多项式P(x)的系数，返回时其中的值将被破坏</param>
        /// <param name="m">多项式P(x)的项数，其最高次数为m-1</param>
        /// <param name="q">存放多项式Q(x)的系数</param>
        /// <param name="n">多项式Q(x)的项数，其最高次数为n-1</param>
        /// <param name="s">返回商多项式S(x)的系数</param>
        /// <param name="k">商多项式S(x)的项数，其最高次数为k-1。其中k=m-n+l</param>
        /// <param name="r">返回余多项式R(x)的系数</param>
        /// <param name="l">余多项式R(x)的项数，其最高次数为l-1。其中l=n-l</param>
        public static void PDIV(double[] p, int m, double[] q, int n, out double[] s, out int k, out double[] r, out int l)
        {
            int i, j, mm, ll;
            k = m - n + 1;
            l = n - 1;
            s = new double[k];
            r = new double[l];
            for (i = 0; i <= k - 1; i++) s[i] = 0;
            if (q[n - 1] == 0) return;
            ll = m - 1;
            for (i = k; i >= 1; i--)
            {
                s[i - 1] = p[ll] / q[n - 1];
                mm = ll;
                for (j = 1; j <= n - 1; j++)
                {
                    p[mm - 1] = p[mm - 1] - s[i - 1] * q[n - j - 1];
                    mm--;
                }
                ll--;
            }
            for (i = 0; i <= l - 1; i++) r[i] = p[i];
        }

        /// <summary>
        /// 求复系数多项式P(z)=Σpz被复系数多项式Q(z)=Σqz除得的商多项式S(z)和余多项式R(z)。
        /// </summary>
        /// <param name="pr">存放多项式P(z)的系数的实部，返回时其中的值将被破坏</param>
        /// <param name="pi">存放多项式P(z)的系数的虚部，返回时其中的值将被破坏</param>
        /// <param name="m">多项式P(z)的项数，其最高次数为m-1</param>
        /// <param name="qr">存放多项式Q(z)的系数的实部</param>
        /// <param name="qi">存放多项式Q(z)的系数的虚部</param>
        /// <param name="n">多项式Q(z)的项数，其最高次数为n-1</param>
        /// <param name="sr">返回商多项式S(z)的系数的实部</param>
        /// <param name="si">返回商多项式S(z)的系数的虚部</param>
        /// <param name="k">商多项式S(z)的项数，其最高次数为k-1。其中k=m-n+l</param>
        /// <param name="rr">返回余多项式R(z)的系数的实部</param>
        /// <param name="ri">返回余多项式R(z)的系数的虚部</param>
        /// <param name="l">余多项式R(z)的项数，其最高次数为l-1。其中l=n-l</param>
        public static void CPDV(double[] pr, double[] pi, int m, double[] qr, double[] qi, int n, out double[] sr, out double[] si, int k, out double[] rr, out double[] ri, out int l)
        {
            int i, j, mm, ll;
            k = m - n + 1;
            l = n - 1;
            sr = new double[k];
            si = new double[k];
            rr = new double[l];
            ri = new double[l];
            double a, b, c, d, u, v;
            for (i = 0; i <= k - 1; i++)
            {
                sr[i] = 0;
                si[i] = 0;
            }
            d = qr[n - 1] * qr[n - 1] + qi[n - 1] * qi[n - 1];
            if (d == 0) return;
            ll = m - 1;
            for (i = k; i >= 1; i--)
            {
                a = pr[ll]; b = pi[ll];c = qr[n - 1];d = qi[n - 1];
                Complex.CDIV(a, b, c, d, out u, out v);
                sr[i - 1] = u;si[i - 1] = v;
                mm = ll;
                for (j = 1; j <= n - 1; j++)
                {
                    a = sr[i - 1];
                    b = si[i - 1];
                    c = qr[n - j - 1];
                    d = qi[n - j - 1];
                    Complex.CMUL(a, b, c, d, out u, out v);
                    pr[mm - 1] = pr[mm - 1] - u;
                    pi[mm - 1] = pi[mm - 1] - v;
                    mm--;
                }
                ll--;
            }
            for (i = 0; i < l; i++)
            {
                rr[i] = pr[i];
                ri[i] = pi[i];
            }
        }
    }
}
