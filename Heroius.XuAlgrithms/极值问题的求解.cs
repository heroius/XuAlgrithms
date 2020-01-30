using System;
using System.Collections.Generic;
using System.Text;

namespace Heroius.XuAlgrithms
{
    /// <summary>
    /// 极值问题的求解
    /// </summary>
    public static class Extremum
    {
        /// <summary>
        /// 用连分式法求目标函数f(x)的极值点
        /// </summary>
        /// <param name="x">x[2]: x(0)存放极值点初值，返回极值点x；x(1)返回极值点x处的函数值</param>
        /// <param name="eps">控制精度要求，一般取10^-10~10^-35之间的数</param>
        /// <param name="k">允许迭代次数的最大值</param>
        /// <param name="js">js[2]: js(0)返回实际迭代次数，若js(0)=k（允许迭代次数的最大值），
        /// 则有可能没有满足精度要求，返回的极值点只作为参考。
        /// js(1)返回标志：若js(1)<0，则说明返回的极值点为极小值点；
        /// 若js(1)>0，则说明返回的极值点为极大值点。</param>
        /// <param name="func">指向计算目标函数f(x)值与倒数f'(x)值的函数
        /// <para>第一参数：输入x的值</para>
        /// <para>返回double[2]，其中：(0)为f(x)；(1)为f'(x)</para>
        /// </param>
        public static void MAX1(ref double[] x, double eps, int k, out int[] js, Func<double, double[]> func)
        {
            double xx, h1, h2=0, dx;
            double[]
                y = new double[10],
                b = new double[10],
                z = new double[2];
            js = new int[2];
            js[0] = 0;
            int i, j, jt = 1, m;
            while (jt==1)
            {
                j = 0;
                while (j<=7)
                {
                    if (j <= 2) xx = x[0] + 0.01 * j;
                    else xx = h2;
                    z = func(xx);
                    if (Math.Abs(z[1])<eps)
                    {
                        jt = 0;
                        j = 10;
                    }
                    else
                    {
                        h1 = z[1];
                        h2 = xx;
                        if(j==0)
                        {
                            y[0] = h1;
                            b[0] = h2;
                        }
                        else
                        {
                            y[j] = h1;
                            m = 0;
                            i = 0;
                            while ((m==0)&&(j<=j-1))
                            {
                                if (h2 == b[i]) m = 1;
                                else h2 = (h1 - y[i]) / (h2 - b[i]);
                                i++;
                            }
                            b[j] = h2;
                            if (m != 0) b[j] = 1e35;
                            h2 = 0;
                            for (i = j - 1; i >= 0; i--)
                                h2 = -y[i] / (b[i + 1] + h2);
                            h2 += b[0];
                        }
                        j++;
                    }
                }
                x[0] = h2;
                js[0]++;
                if (js[0] == k) jt = 0;
            }
            xx = x[0];
            z = func(xx);
            x[1] = z[0];
            if (Math.Abs(x[0]) <= 1) dx = 1e-5;
            else dx = Math.Abs(x[0] * 1e-5);
            xx = x[0] - dx;
            z = func(xx);
            h1 = z[0];
            xx = x[0] + dx;
            z = func(xx);
            h2 = z[0];
            js[1] = -1;
            if ((h1 + h2 - 2 * x[1]) <= 0) js[1] = 1;
        }

        /// <summary>
        /// 用连分式法求多元函数的极值点与极值点处的函数值。
        /// </summary>
        /// <param name="x">x[n+1]: 前n个分量存放极值点初值，返回极值点的n个坐标；最后一个分量返回极值点处的函数值</param>
        /// <param name="n">自变量个数</param>
        /// <param name="eps">控制精度要求</param>
        /// <param name="k">允许最大迭代次数</param>
        /// <param name="js">js[2]: </param>
        /// <param name="func">指向计箕目标函数值以及各偏导数值的函数名
        /// <para>第一参数：x[] 存放自变量</para>
        /// <para>第二参数：n 自变量个数</para>
        /// <para>第三参数：j 指定计算标记。当j=0 计算f(x)；其他值n计算第n个变量的偏导数值</para>
        /// <para>返回值：根据第三参数j的计算结果</para>
        /// </param>
        public static void MAXN(ref double[] x, int n, double eps, int k, out double[] js, Func<double[], int, int, double> func)
        {
            int i, j, m, l, jt, il;
            double p, zz=0, t, h1, h2, f, dx;
            double[] y = new double[10], b = new double[10];
            js = new double[2];
            js[0] = 0;
            jt = 1; h2 = 0.0;
            while (jt == 1)
            {
                t = 0.0; js[0] = js[0] + 1;
                for (i = 1; i <= n; i++)
                {
                    f = func(x,n,i);
                    t = t + Math.Abs(f);
                }
                if (t < eps) jt = 0;
                else
                {
                    for (i = 0; i <= n - 1; i++)
                    {
                        il = 5;
                        while (il != 0)
                        {
                            j = 0; t = x[i]; il = il - 1;
                            while (j <= 7)
                            {
                                if (j <= 2) zz = t + j * 0.01;
                                else zz = h2;
                                x[i] = zz;
                                f = func(x, n, i + 1);
                                if (Math.Abs(f) + 1.0 == 1.0) { j = 10; il = 0; }
                                else
                                {
                                    h1 = f; h2 = zz;
                                    if (j == 0) { y[0] = h1; b[0] = h2; }
                                    else
                                    {
                                        y[j] = h1; m = 0; l = 0;
                                        while ((m == 0) && (l <= j - 1))
                                        {
                                            p = h2 - b[l];
                                            if (Math.Abs(p) + 1.0 == 1.0) m = 1;
                                            else h2 = (h1 - y[l]) / p;
                                            l = l + 1;
                                        }
                                        b[j] = h2;
                                        if (m != 0) b[j] = 1.0e+35;
                                        h2 = 0.0;
                                        for (l = j - 1; l >= 0; l--)
                                            h2 = -y[l] / (b[l + 1] + h2);
                                        h2 = h2 + b[0];
                                    }
                                    j = j + 1;
                                }
                            }
                            x[i] = h2;
                        }
                        x[i] = zz;
                    }
                    if (js[0] == 10) jt = 0;
                }
            }
            js[1] = 1;
            f = func(x, n, 0);
            x[n] = f;
            dx = 0.00001; t = x[0];
            x[0] = t + dx; h1 = func(x, n, 0);
            x[0] = t - dx; h2 = func(x, n, 0);
            x[0] = t;
            t = h1 + h2 - 2.0 * f;
            if (t > 0.0) js[1] = -1;
            j = 1; jt = 1;
            while (jt == 1)
            {
                j = j + 1; dx = 0.00001; jt = 0;
                t = x[j - 1];
                x[j - 1] = t + dx; h2 = func(x, n, 0);
                x[j - 1] = t - dx; h1 = func(x, n, 0);
                x[j - 1] = t; t = h1 + h2 - 2.0 * f;
                if ((t * js[1] < 0.0) && (j < n)) jt = 1;
            }
            if (t * js[1] > 0.0) js[1] = 0;
        }

        /// <summary>
        /// 求解不等式约束条件下的线性规划问题。
        /// </summary>
        /// <param name="a" xml:space="preserve">a[m,m+n]: 存放矩阵A。
        /// <para>a00&#x20;&#x20;&#x20;&#x20;&#x20;&#x20;a01&#x20;&#x20;&#x20;&#x20;&#x20;&#x20;...&#x20;&#x20;&#x20;&#x20;a0,n-1&#x20;&#x20;&#x20;&#x20;&#x20;&#x20;1&#x20;&#x20;&#x20;&#x20;&#x20;&#x20;0&#x20;&#x20;&#x20;&#x20;&#x20;&#x20;...&#x20;&#x20;&#x20;&#x20;0</para>
        /// <para>a10&#x20;&#x20;&#x20;&#x20;&#x20;&#x20;a11&#x20;&#x20;&#x20;&#x20;&#x20;&#x20;...&#x20;&#x20;&#x20;&#x20;a1,n-1&#x20;&#x20;&#x20;&#x20;&#x20;&#x20;0&#x20;&#x20;&#x20;&#x20;&#x20;&#x20;1&#x20;&#x20;&#x20;&#x20;&#x20;&#x20;...&#x20;&#x20;&#x20;&#x20;0</para>
        /// <para>...&#x20;&#x20;&#x20;&#x20;&#x20;&#x20;...&#x20;&#x20;&#x20;&#x20;&#x20;&#x20;&#x20;&#x20;&#x20;&#x20;&#x20;&#x20;&#x20;...&#x20;&#x20;&#x20;&#x20;&#x20;&#x20;&#x20;&#x20;&#x20;...&#x20;&#x20;&#x20;&#x20;...&#x20;&#x20;&#x20;&#x20;...&#x20;&#x20;...</para>
        /// <para>am-1,0&#x20;&#x20;&#x20;am-1,1&#x20;&#x20;&#x20;...&#x20;&#x20;&#x20;&#x20;am-1,n-1&#x20;&#x20;&#x20;&#x20;1&#x20;&#x20;&#x20;&#x20;&#x20;&#x20;0&#x20;&#x20;&#x20;&#x20;&#x20;&#x20;...&#x20;&#x20;&#x20;&#x20;1</para>
        /// </param>
        /// <param name="b">b[m]: 存放不等式约束条件的右端项值b0,b1,...,bm-1</param>
        /// <param name="c">c[m+n]: 存放目标函数中的系数，其中后m个分量为0</param>
        /// <param name="m">不等式约束条件的个数</param>
        /// <param name="n">变量个数</param>
        /// <param name="x">x[m+n]: 前n个分量返回目标函数f的极小值点的n个坐标，
        /// 第n+1个分量返回目标函数f的极小值，其余为本函数的工作单元</param>
        /// <returns>函数返回标志值。
        /// 若返回的标志值为0, 则表示目标函数值无界；
        /// 若返回的标志值不为0, 则表示正常返回</returns>
        public static int LPLQ(double[,] a, double[] b, double[] c, int m, int n, out double[] x)
        {
            int i, mn, k, j;
            double s, z, dd, y;
            int[] jjs = new int[m];
            double[,]
                p = new double[m, m],
                d = new double[m, m + n];
            x = new double[m + n];

            for (i = 0; i <= m - 1; i++) jjs[i] = n + i;
            mn = m + n; s = 0.0;
            while (true)
            {
                for (i = 0; i <= m - 1; i++)
                    for (j = 0; j <= m - 1; j++) 
                        p[i,j] = a[i,jjs[j]];
                if (Matrix.RINV(ref p, m) == 0)
                {
                    x[n] = s; return 1;
                }
                Matrix.TRMUL(p, a, m, m, mn, out d);
                for (i = 0; i <= mn - 1; i++) x[i] = 0.0;
                for (i = 0; i <= m - 1; i++)
                {
                    s = 0.0;
                    for (j = 0; j <= m - 1; j++) s = s + p[i,j] * b[j];
                    x[jjs[i]] = s;
                }
                k = -1; dd = 1.0e-35;
                for (j = 0; j <= mn - 1; j++)
                {
                    z = 0.0;
                    for (i = 0; i <= m - 1; i++) z = z + c[jjs[i]] * d[i,j];
                    z = z - c[j];
                    if (z > dd) { dd = z; k = j; }
                }
                if (k == -1)
                {
                    s = 0.0;
                    for (j = 0; j <= n - 1; j++) s = s + c[j] * x[j];
                    x[n] = s;
                    return 1;
                }
                j = -1;
                dd = 1.0e+20;
                for (i = 0; i <= m - 1; i++)
                    if (d[i,k] >= 1.0e-20)
                    {
                        y = x[jjs[i]] / d[i,k];
                        if (y < dd) { dd = y; j = i; }
                    }
                if (j == -1)
                {
                    x[n] = s;
                    return 0;
                }
                jjs[j] = k;
            }
        }

        /// <summary>
        /// 用单形调优法求解无约束条件下的n 维极值问题。
        /// </summary>
        /// <param name="n">变量个数</param>
        /// <param name="d">初始单形中任意两顶点间的距离</param>
        /// <param name="u">扩张系数μ。一般取 1.2＜μ＜2.0</param>
        /// <param name="v">收缩系数λ。一般取 0.0＜λ＜1.0</param>
        /// <param name="x">x[n+1]: 前n个分量返回极小值点的n个坐标，最后一个分量返回极小值。</param>
        /// <param name="eps">控制精度要求</param>
        /// <param name="k">允许的最大迭代次数</param>
        /// <param name="xx">xx[n,n+1]: 返回最后单形的n+1个顶点坐标(x0i,x1i,...,xn-1,i)(i=0,1,...,n)</param>
        /// <param name="f">f[n+1]: 返回最后单形的n+1个顶点的目标函数值</param>
        /// <param name="s">计算目标函数值的函数
        /// <para>第一参数：x[n] 自变量的值</para>
        /// <para>第二参数：n 自变量个数</para>
        /// <para>返回值：y 目标函数值</para>
        /// </param>
        /// <returns>函数返回实际迭代次数。若实际迭代次数等于允许的最大迭代次数k，则有可能未达到精度要求，返回的极值点只作为参考</returns>
        public static int JSIM(int n, double d, double u, double v, out double[] x, double eps, int k, out double[,] xx, out double[] f, Func<double[], int, double> s)
        {
            //
        }
    }
}
