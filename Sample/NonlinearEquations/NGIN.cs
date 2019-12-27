using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.NonlinearEquations
{
    public class NGIN : SampleItem
    {
        public NGIN() : base("NonlinearEquations.NGIN") { }

        public override string Execute()
        {
            int m = 2, n = 2, ka = 3;
            double eps1 = 0.000001, eps2 = 0.000001;
            double[] x = new double[] { 0.5, -1 };
            var i = Heroius.XuAlgrithms.NonlinearEquations.NGIN(m, n, eps1, eps2, x, ka, F, S);
            StringBuilder builder = new StringBuilder();
            builder.AppendLine($"i={i}");
            for (i = 0; i <= 1; i++)
            {
                builder.AppendLine($"x({i}) = {x[i].ToString("E")}");
            }
            return builder.ToString();
        }

        public void F(int m, int n, double[] x, ref double[] d)
        {
            d[0] = x[0] * x[0] + 10 * x[0] * x[1] + 4 * x[1] * x[1] + 0.7401006;
            d[1] = x[0] * x[0] - 3 * x[0] * x[1] + 2 * x[1] * x[1] - 1.0201228;
        }

        public void S(int m, int n, double[] x, ref double[] p)
        {
            double[,] p_t = new double[2, 2];
            p_t[0,0] = 2 * x[0] + 10 * x[1];
            p_t[0,1] = 10 * x[0] + 8 * x[1];
            p_t[1, 0] = 2 * x[0] - 3 * x[1];
            p_t[1, 1] = -3 * x[0] + 4 * x[1];
            p = Heroius.XuAlgrithms.Utility.C.Convert(p_t);
        }
    }
}
