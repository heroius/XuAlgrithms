using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.ODEs
{
    public class TNR : SampleItem
    {
        public TNR() : base("ODEs.TNR") { }

        public override string Execute()
        {
            int n = 3;
            double[] y = new double[] { 1.0, 0.0, -1 };
            double t = 0.0, h = 0.001;

            StringBuilder builder = new StringBuilder();
            builder.AppendLine($"t={t}");
            for (int i = 0; i < n; i++)
            {
                builder.Append($"y({i})={y[i]}\t\t");
            }
            builder.AppendLine();
            for (int i = 0; i < 10; i++)
            {
                Heroius.XuAlgrithms.ODEs.TNR(t, h, n, ref y, Demo);
                t += h;
                builder.AppendLine($"t={t}");

                for (int j = 0; j < n; j++)
                {
                    builder.Append($"y({j})={y[j]}\t\t");
                }
                builder.AppendLine();
            }
            return builder.ToString();
        }

        public double[] Demo(double t, double[] y, int n)
        {
            return new double[] {
                -21*y[0]+19*y[1]-20*y[2],
                19*y[0]-21*y[1]+20*y[2],
                40*y[0] -40*y[1]-40*y[2]
            };
        }
    }
}
