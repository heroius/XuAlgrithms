using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.ODEs
{
    public class GIL : SampleItem
    {
        public GIL() : base("ODEs.GIL") { }

        public override string Execute()
        {
            int n = 3, k = 11;
            double[] q = new double[] { 0,0,0};
            double[] y = new double[] { 0.0, 1.0, 1.0 };
            double t = 0.0, h = 0.1, eps = 0.000001, x;
            double[,] z;

            StringBuilder builder = new StringBuilder();
            for (int i = 0; i < n; i++)
            {
                builder.Append($"y({i})={y[i]}");
            }
            for (int j = 1; j < k; j++)
            {
                Heroius.XuAlgrithms.ODEs.GIL(t, h, ref y, n, eps, ref q, Demo);
                t = j * h;
                builder.AppendLine($"t={t}");
                for (int i = 0; i < n; i++)
                {
                    builder.Append($"y({i})={y[i]}\t\t");
                }
                builder.AppendLine();
            }
            return builder.ToString();
        }

        public double[] Demo(double t, double[] y, int n)
        {
            return new double[] {
                y[1], -y[0], -y[2]
            };
        }
    }
}
