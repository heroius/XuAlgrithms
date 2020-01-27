using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.ODEs
{
    public class GJFQ : SampleItem
    {
        public GJFQ() : base("ODEs.GJFQ") { }

        public override string Execute()
        {
            int n = 2, k = 11;
            double[] y = new double[] { 1.0, 0.0 };
            double t = 0.0, h = 0.1, eps = 0.00001;
            double[,] z;

            StringBuilder builder = new StringBuilder();
            Heroius.XuAlgrithms.ODEs.GJFQ(t, h, n, y, eps, k, out z, Demo);
            for (int i = 0; i < k; i++)
            {
                t = i * h;
                builder.AppendLine($"t={t}\t\ty(0)={z[0,i]}\t\ty(1)={z[1,i]}");
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
