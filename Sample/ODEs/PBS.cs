using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.ODEs
{
    public class PBS : SampleItem
    {
        public PBS() : base("ODEs.PBS") { }

        public override string Execute()
        {
            int n = 2;
            double[] y = new double[] { 1.0, 0.0 };
            double t = 0.0, h = 0.1, eps = 0.000001;

            StringBuilder builder = new StringBuilder();
            builder.AppendLine($"t={t}\t\ty(0)={y[0]}\t\ty(1)={y[1]}");
            for (int i = 0; i < 10; i++)
            {
                Heroius.XuAlgrithms.ODEs.PBS(t, h, n, ref y, eps, Demo);
                t += h;
                builder.AppendLine($"t={t}\t\ty(0)={y[0]}\t\ty(1)={y[1]}");
            }
            return builder.ToString();
        }

        public double[] Demo(double t, double[] y, int n)
        {
            return new double[] {
                -y[1],
                y[0]
            };
        }
    }
}
