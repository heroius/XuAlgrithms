using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Text;

namespace Heroius.XuAlgrithms
{
    /// <summary>
    /// 求m×n阶矩阵A与n×k阶矩阵B的乘积矩阵C=AB.
    /// </summary>
    [Description("实矩阵相乘")]
    public class TRMUL
    {
        int m, n, k;
        double[,] a, b, c;

        /// <summary>
        /// 求m×n阶矩阵A与n×k阶矩阵B的乘积矩阵C=AB.
        /// </summary>
        /// <param name="mm">m</param>
        /// <param name="nn">n</param>
        /// <param name="kk">k</param>
        public TRMUL(int mm, int nn, int kk)
        {
            m = mm; n = nn; k = kk;
            a = new double[m, n];
            b = new double[n, k];
            c = new double[m, k];
        }

        /// <summary>
        /// 输入参数
        /// </summary>
        /// <param name="aa">存放矩阵A的元素</param>
        /// <param name="bb">存放矩阵B的元素</param>
        public void Input(double[,] aa, double[,] bb)
        {
            int i, j;
            for (i = 0; i < m; i++)                       //读入矩阵A
                for (j = 0; j < n; j++) a[i, j] = aa[i, j];
            for (i = 0; i < n; i++)                       //读入矩阵B
                for (j = 0; j < k; j++) b[i, j] = bb[i, j];
        }

        /// <summary>
        /// 执行C=AB
        /// </summary>
        public void MUL()
        {
            int i, j, t;
            for (i = 0; i < m; i++)
                for (j = 0; j < k; j++)
                {
                    c[i, j] = 0.0;
                    for (t = 0; t < n; t++)
                        c[i, j] = c[i, j] + a[i, t] * b[t, j];
                }
        }

        /// <summary>
        /// 获取乘积矩阵C
        /// </summary>
        /// <returns></returns>
        public double[,] Output()
        {
            return c;
        }
    }
}
