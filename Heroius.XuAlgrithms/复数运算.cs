using System;
using System.Collections.Generic;
using System.Text;

namespace Heroius.XuAlgrithms
{
    /// <summary>
    /// 复数运算
    /// </summary>
    public static class Complex
    {
        /// <summary>
        /// 计算两个复数的乘积u+jv=(a+jb)(c+jd)，其中j=√(-1)。
        /// </summary>
        /// <param name="a">表示复数a+jb</param>
        /// <param name="b">表示复数a+jb</param>
        /// <param name="c">表示复数c+jd</param>
        /// <param name="d">表示复数c+jd</param>
        /// <param name="u">指向返回的复数乘积u+jv</param>
        /// <param name="v">指向返回的复数乘积u+jv</param>
        public static void CMUL(double a, double b, double c, double d, out double u, out double v)
        {
            double p = a * c;
            double q = b * d;
            double s = (a + b) * (c + d);
            u = p - q;
            v = s - p - q;
        }

        /// <summary>
        /// 计算两个复数的商u+jv=(a+jb)/(c+jd)，其中j=√(-1)。
        /// </summary>
        /// <param name="a">表示复数a+jb</param>
        /// <param name="b">表示复数a+jb</param>
        /// <param name="c">表示复数c+jd</param>
        /// <param name="d">表示复数c+jd</param>
        /// <param name="u">指向返回的复数商u+jv</param>
        /// <param name="v">指向返回的复数商u+jv</param>
        public static void CDIV(double a, double b, double c, double d, out double u, out double v)
        {
            double
                p = a * c,
                q = -b * d,
                s = (a + b) * (c - d),
                w = c * c + d * d;
            if (w==0)
            {
                u = 1e35 * a / Math.Abs(a);
                v = 1e35 * b / Math.Abs(b);
            }
            else
            {
                u = (p - q) / w;
                v = (s - p - q) / w;
            }
        }

        /// <summary>
        /// 计算复变量的整数次幕u+jv=(x+jy)^n，其中n为整数，j=√(-1)。
        /// </summary>
        /// <param name="x">表示复数x+jy</param>
        /// <param name="y">表示复数x+jy</param>
        /// <param name="n">幕的次数</param>
        /// <param name="u">返回的复数u+jv</param>
        /// <param name="v">返回的复数u+jv</param>
        public static void POWR(double x, double y, int n, out double u, out double v)
        {
            double q = Math.Atan2(y, x);
            double r = Math.Sqrt(x * x + y * y);
            if (r!=0)
            {
                r = n * Math.Log(r);
                r = Math.Exp(r);
            }
            u = r * Math.Cos(n * q);
            v = r * Math.Sin(n * q);
        }

        /// <summary>
        /// 计算复变量n次方根 u+jv=(x+jy)^(1/n)，其中n为正整数，j=√(-1)。
        /// </summary>
        /// <param name="x">表示复数x+jy</param>
        /// <param name="y">表示复数x+jy</param>
        /// <param name="n">幕的次数</param>
        /// <param name="u">返回的复数z=x+jy的n次方根的n个值的实部</param>
        /// <param name="v">返回的复数z=x+jy的n次方根的n个值的虚部</param>
        public static void NTRT(double x, double y, int n, out double[] u, out double[] v)
        {
            u = new double[n];
            v = new double[n];
            if (n < 1) throw new Exception("n应为正整数");
            double
                q = Math.Atan2(y, x),
                r = Math.Sqrt(x * x + y * y);
            if (r != 0)
            {
                r = (1.0 / n) * Math.Log(r);
                r = Math.Exp(r);
            }
            for(int k = 0; k < n; k++)
            {
                var t = (2.0 * k * Math.PI + q) / n;
                u[k] = r * Math.Cos(t);
                v[k] = r * Math.Sin(t);
            }
        }

        /// <summary>
        /// 计算复变量的指数 u+jv=e^(x+jy)，其中 j=√(-1)。
        /// </summary>
        /// <param name="x">表示复数x+jy</param>
        /// <param name="y">表示复数x+jy</param>
        /// <param name="u">返回e^(x+jy)的实部</param>
        /// <param name="v">返回e^(x+jy)的虚部</param>
        public static void CEXP(double x, double y, out double u, out double v)
        {
            double p = Math.Exp(x);
            u = p * Math.Cos(y);
            v = p * Math.Sin(y);
        }

        /// <summary>
        /// 计算复变量的自然对数 u+jv=ln(x+jy)，其中 j=√(-1)。
        /// </summary>
        /// <param name="x">表示复数x+jy</param>
        /// <param name="y">表示复数x+jy</param>
        /// <param name="u">返回ln(x+jy)的实部</param>
        /// <param name="v">返回ln(x+jy)的虚部</param>
        public static void CLOG(double x, double y, out double u, out double v)
        {
            double p = Math.Log(Math.Sqrt(x * x + y * y));
            u = p;
            v = Math.Atan2(y, x);
        }

        /// <summary>
        /// 计算复变量的正弦值 u+jv=sin(x+jy)，其中 j=√(-1)。
        /// </summary>
        /// <param name="x">表示复数x+jy</param>
        /// <param name="y">表示复数x+jy</param>
        /// <param name="u">返回sin(x+jy)的实部</param>
        /// <param name="v">返回sin(x+jy)的虚部</param>
        public static void CSIN(double x, double y, out double u, out double v)
        {
            double
                p = Math.Exp(y),
                q = Math.Exp(-y);
            u = Math.Sin(x) * (p + q) / 2;
            v = Math.Cos(x) * (p - q) / 2;
        }

        /// <summary>
        /// 计算复变量的余弦值 u+jv=cos(x+jy)，其中 j=√(-1)。
        /// </summary>
        /// <param name="x">表示复数x+jy</param>
        /// <param name="y">表示复数x+jy</param>
        /// <param name="u">返回cos(x+jy)的实部</param>
        /// <param name="v">返回cos(x+jy)的虚部</param>
        public static void CCOS(double x, double y, out double u, out double v)
        {
            double
                p = Math.Exp(y),
                q = Math.Exp(-y);
            u = Math.Cos(x) * (p + q) / 2;
            v = -Math.Sin(x) * (p - q) / 2;
        }
    }
}
