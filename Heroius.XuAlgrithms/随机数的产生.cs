using System;
using System.Collections.Generic;
using System.Text;

namespace Heroius.XuAlgrithms
{
    /// <summary>
    /// 随机数的产生
    /// </summary>
    public static class Random
    {
        /// <summary>
        /// 产生0到1之间均匀分布的一个随机数
        /// </summary>
        /// <param name="r">随机数种子</param>
        /// <returns>返回0到1之间均匀分布的一个随机数</returns>
        public static double RND1(double r)
        {
            double s = 65536, u = 2053, v = 13849;
            double m = r / s;
            r = r - m * s;
            r = u * r + v;
            m = r / s;
            r = r - m * s;
            return r / s;
        }

        /// <summary>
        /// 产生0~1之间均匀分布的随机数序列。
        /// </summary>
        /// <param name="r">随机数种子</param>
        /// <param name="n">随机数序列长度</param>
        /// <returns>返回随机数序列</returns>
        public static double[] RNDS(double r, int n)
        {
            double[] p = new double[n];
            double s = 65536, u = 2053, v = 13849;
            double m;
            for (int i = 0; i < n; i++)
            {
                r = u * r + v;
                m = r / s;
                r = r - m * s;
                p[i] = r / s;
            }
            return p;
        }

        /// <summary>
        /// 产生给定区间[a,b] 内均匀分布的一个随机整数
        /// </summary>
        /// <param name="a">随机整数所在的区间左</param>
        /// <param name="b">随机整数所在的区间右</param>
        /// <param name="r">向随机数种子，应为大于零的奇数</param>
        /// <returns>区间[a,b] 内均匀分布的一个随机整数</returns>
        public static int RAB1(int a, int b, int r)
        {
            int k = b - a + 1, l = 2;
            int m, i, p;
            while (l<k)
            {
                l++;
            }
            m = 4 * l;
            k = r;
            i = 1;
            while (i<=1)
            {
                k *= 5;
                k = k % m;
                l = k / 4 + a;
                if (l <= b)
                {
                    p = l;
                    i++;
                }
            }
            r = k;
            return p;
        }
    }
}
