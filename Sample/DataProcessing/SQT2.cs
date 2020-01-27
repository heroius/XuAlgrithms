using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.DataProcessing
{
    public class SQT2 : SampleItem
    {
        public SQT2() : base("DataProcessing.SQT2") { }

        public override string Execute()
        {
            int m = 3, n = 5;
            double[,] x = new double[,] { { 1.1, 1.0, 1.2, 1.1, 0.9 }, { 2.0, 2.0, 1.8, 1.9, 2.1 }, { 3.2, 3.2, 3.0, 2.9, 2.9 } };
            double[] y = new double[] { 10.1, 10.2, 10.0, 10.1, 10.0 };
            Heroius.XuAlgrithms.DataProcessing.SQT2(x, y, m, n, out double[] a, out double[] dt, out double[] v);

            StringBuilder builder = new StringBuilder();
            for (int i = 0; i < n-1; i++)
            {
                builder.AppendLine($"a({i})={a[i]}");
            }
            builder.AppendLine($"q={dt[0]}\t\ts={dt[1]}\t\tr={dt[2]}");
            for (int i = 0; i < m; i++)
            {
                builder.AppendLine($"v({i})={v[i]}");
            }
            builder.AppendLine($"u={dt[3]}");
            return builder.ToString();
        }
    }
}
