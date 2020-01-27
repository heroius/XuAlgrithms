using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.ODEs
{
    public class ADMS : SampleItem
    {
        public ADMS() : base("ODEs.ADMS") { }

        public override string Execute()
        {
            int n = 3, k = 11;
            double[] y = new double[] { 0, 1.0, 1 };
            double t = 0.0, h = 0.05, eps = 0.0001;
            double[,] z;

            StringBuilder builder = new StringBuilder();
            Heroius.XuAlgrithms.ODEs.ADMS(t, h, n, y, eps, k, out z, Demo);
            for (int i = 0; i < k; i++)
            {
                t = i * h;
                builder.AppendLine($"t={t}");
                for (int j = 0; j < 3; j++)
                {
                    builder.Append($"y({j}) ={ z[j, i]}\t\t");
                }
                builder.AppendLine();
            }
            return builder.ToString();
        }

        public double[] Demo(double t, double[] y, int n)
        {
            return new double[] {
                y[1],
                -y[0],
                -y[2]
            };
        }
    }
}
