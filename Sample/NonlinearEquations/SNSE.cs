using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.NonlinearEquations
{
    public class SNSE : SampleItem
    {
        public SNSE() : base("NonlinearEquations.SNSE") { }

        public override string Execute()
        {
            double[] x = new double[] { 1.5, 6.5, -5 };
            int js = 500;
            var i = Heroius.XuAlgrithms.NonlinearEquations.SNSE(3, 0.000001, x, js, Demo);
            StringBuilder builder = new StringBuilder();
            if (i>0 && i< js)
            {
                for (i = 0; i <= 2; i++)
                {
                    builder.AppendLine($"x({i}) = {x[i].ToString("E")}");
                }
            }
            return builder.ToString();
        }

        /// <summary>
        /// 给定一个n维点x，
        /// </summary>
        /// <param name="x">给定的点</param>
        /// <param name="y">输出偏导数值</param>
        /// <param name="n">方程数</param>
        /// <returns></returns>
        public double Demo(double[] x, double[] y, int n)
        {
            double f1 = x[0] - 5.0 * x[1] * x[1] + 7.0 * x[2] * x[2] + 12;
            double f2 = 3.0 * x[0] * x[1] + x[0] * x[2] - 11 * x[0];
            double f3 = 2.0 * x[1] * x[2] + 40 * x[0];
            double z = f1 * f1 + f2 * f2 + f3 * f3;//目标函数在x处的取值
            
            //对x[0]求导
            double df1 = 1;
            double df2 = 3 * x[1] + x[2] - 11;
            double df3 = 40;
            y[0] = 2.0 * (f1 * df1 + f2 * df2 + f3 * df3);
            
            df1 = 10 * x[1];
            df2 = 3.0 * x[0];
            df3 = 2 * x[2];
            y[1] = 2.0 * (f1 * df1 + f2 * df2 + f3 * df3);

            df1 = 14 * x[2];
            df2 = x[0];
            df3 = 2 * x[1];
            y[2] = 2.0 * (f1 * df1 + f2 * df2 + f3 * df3);

            return z;
        }
    }
}
