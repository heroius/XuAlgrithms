using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.DataProcessing
{
    public class SQT1 : SampleItem
    {
        public SQT1() : base("DataProcessing.SQT1") { }

        public override string Execute()
        {
            int n = 11;
            double[] x = new double[] { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
            double[] y = new double[] { 2.75, 2.84, 2.965, 3.01, 3.20, 3.25, 3.38, 3.43, 3.55, 3.66, 3.74 };
            Heroius.XuAlgrithms.DataProcessing.SQT1(x, y, n, out double[] a, out double[] dt);

            StringBuilder builder = new StringBuilder();
            builder.AppendLine($"n={n}");
            builder.AppendLine($"a={a[1]}\t\tb={a[0]}");
            builder.AppendLine($"q={dt[0]}\t\ts={dt[1]}\t\tp={dt[2]}");
            builder.AppendLine($"umax={dt[3]}\t\tumin={dt[4]}\t\tu={dt[5]}");
            return builder.ToString();
        }
    }
}
