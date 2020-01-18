using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.ODEs
{
    public class MRSN : SampleItem
    {
        public MRSN() : base("ODEs.MRSN") { }

        public override string Execute()
        {
            int n = 2, k = 11;
            double[] y = new double[] { 0.0, 1.0 };
            double t = 0.0, h = 0.1, eps = 0.00001, x;
            double[,] z;
            Heroius.XuAlgrithms.ODEs.MRSN(t, h, n, y, eps, k, out z, Demo);
            StringBuilder builder = new StringBuilder();
            for (int i = 0; i < k; i++)
            {
                x = i * h;
                builder.AppendLine($"t={x}\t\ty0={z[0, i]}\t\ty1={z[1, i]}");
            }
            return builder.ToString();
        }

        public double[] Demo(double t, double[] y, int n)
        {
            double q = 60.0 * (0.06 + t * (t - 0.6));
            return new double[] {
                q*y[1],
                -q*y[0]
            };
        }
    }
}
