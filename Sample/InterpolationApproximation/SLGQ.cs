using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.InterpolationApproximation
{
    public class SLGQ : SampleItem
    {
        public SLGQ() : base("InterpolationApproximation.SLGQ") { }

        public override string Execute()
        {
            int i, j;
            double u = 0.35, v = 0.65, w;
            double[] x = new double[11], y = new double[11];
            double[,] z = new double[11,11];
            for (i = 0; i <= 10; i++)
            { x[i] = 0.1 * i; y[i] = x[i]; }
            for (i = 0; i <= 10; i++)
                for (j = 0; j <= 10; j++)
                    z[i,j] = Math.Exp(-(x[i] - y[j]));
            w = Heroius.XuAlgrithms.InterpolationApproximation.SLGQ(x, y, z, 11, 11, u, v);
            StringBuilder builder = new StringBuilder();
            builder.AppendLine($"x={u}, y={v}, z={w}");
            u = 0.45; v = 0.55;
            w = Heroius.XuAlgrithms.InterpolationApproximation.SLGQ(x, y, z, 11, 11, u, v);
            builder.AppendLine($"x={u}, y={v}, z={w}");
            return builder.ToString();
        }
    }
}
