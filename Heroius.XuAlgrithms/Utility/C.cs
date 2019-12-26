using System;
using System.Collections.Generic;
using System.Text;

namespace Heroius.XuAlgrithms.Utility
{
    /// <summary>
    /// 提供C语言特性的相关方法，以保持在C#中的兼容性
    /// </summary>
    public static class C
    {
        /// <summary>
        /// 将二维数组转为一维数组
        /// </summary>
        /// <param name="a">二维数组</param>
        /// <returns>一维数组</returns>
        public static double[] Convert(double[,] a)
        {
            int i = a.GetLength(0);
            int j = a.GetLength(1);
            double[] r = new double[i * j];
            int itr = 0;
            for (int ii = 0; ii < i; ii++)
            {
                for (int jj = 0; jj < j; jj++)
                {
                    r[itr] = a[ii, jj];
                    itr++;
                }
            }
            return r;
        }

        /// <summary>
        /// 将一维数组转为二维数组
        /// </summary>
        /// <param name="a">一维数组</param>
        /// <param name="i">二维数组的第一维度长</param>
        /// <param name="j">二维数组的第二维度长</param>
        /// <returns>二维数组</returns>
        public static double[,] Convert(double[] a, int i, int j)
        {
            double[,] r = new double[i, j];
            int itr = 0;
            for (int ii = 0; ii < i; ii++)
            {
                for (int jj = 0; jj < j; jj++)
                {
                    r[ii, jj] = a[itr];
                    itr++;
                }
            }
            return r;
        }
    }
}
