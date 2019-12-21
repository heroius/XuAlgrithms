using System;
using System.Collections.Generic;
using System.Text;

namespace Heroius.XuAlgrithms
{
    public class NonlinearEquations
    {
        /// <summary>
        /// 用对分法搜索方程f(x)在区间[a,b]内的实根
        /// </summary>
        /// <param name="a">求根区间的左端点</param>
        /// <param name="b">求根区间的右端点</param>
        /// <param name="h">搜索求根时采用的步长</param>
        /// <param name="eps">控制精度要求</param>
        /// <param name="m">在区间[a,b]内实根个数的预估值</param>
        /// <param name="func">计算方程左端函数值f(x)的函数</param>
        /// <returns>在区间[a,b]内实际搜索到的的实根</returns>
        public static double[] DHRT(double a, double b, double h, double eps, int m, Func<double, double> func)
        {
            int js;
            double z, y, z1, y1, z0, y0;
            int n = 0; z = a; y = func(z);
            double[] x = new double[m];
            while ((z <= b + h / 2.0) && (n != m))
            {
                if (Math.Abs(y) < eps)
                {
                    n = n + 1; x[n - 1] = z;
                    z = z + h / 2.0; y = func(z);
                }
                else
                {
                    z1 = z + h; y1 = func(z1);
                    if (Math.Abs(y1) < eps)
                    {
                        n = n + 1; x[n - 1] = z1;
                        z = z1 + h / 2.0; y = func(z);
                    }
                    else if (y * y1 > 0.0)
                    { y = y1; z = z1; }
                    else
                    {
                        js = 0;
                        while (js == 0)
                        {
                            if (Math.Abs(z1 - z) < eps)
                            {
                                n = n + 1; x[n - 1] = (z1 + z) / 2.0;
                                z = z1 + h / 2.0; y = func(z);
                                js = 1;
                            }
                            else
                            {
                                z0 = (z1 + z) / 2.0; y0 = func(z0);
                                if (Math.Abs(y0) < eps)
                                {
                                    x[n] = z0; n = n + 1; js = 1;
                                    z = z0 + h / 2.0; y = func(z);
                                }
                                else if ((y * y0) < 0.0)
                                { z1 = z0; y1 = y0; }
                                else { z = z0; y = y0; }
                            }
                        }
                    }
                }
            }
            return x;
        }
    }
}
