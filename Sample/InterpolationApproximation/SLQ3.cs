using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.InterpolationApproximation
{
    public class SLQ3 : SampleItem
    {
        public SLQ3() : base("InterpolationApproximation.SLQ3") { }

        public override string Execute()
        {
            int i, j;
            double u = 0.9, v = 0.8, w;
            double[] x = new double[6], y = new double[5];
            double[,] z = new double[6, 5];
            for (i = 0; i <= 5; i++) x[i] = 0.2 * i;
            for (j = 0; j <= 4; j++) y[j] = 0.25 * j;
            for (i = 0; i <= 5; i++)
                for (j = 0; j <= 4; j++)
                    z[i, j] = Math.Exp(-(x[i] - y[j]));
            w = Heroius.XuAlgrithms.InterpolationApproximation.SLQ3(x, y, z, 6, 5, u, v);
            StringBuilder builder = new StringBuilder();
            builder.AppendLine($"x={u}, y={v}, z={w}");
            u = 0.3; v = 0.9;
            w = Heroius.XuAlgrithms.InterpolationApproximation.SLQ3(x, y, z, 6, 5, u, v);
            builder.AppendLine($"x={u}, y={v}, z={w}");
            return builder.ToString();
        }
    }
}
