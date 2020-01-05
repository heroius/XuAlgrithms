using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.InterpolationApproximation
{
    public class PIR2 : SampleItem
    {
        public PIR2() : base("InterpolationApproximation.PIR2") { }

        public override string Execute()
        {
            int i, j;
            double[] x = new double[11], y = new double[21], dt;
            double[,] z = new double[11, 21], a;
            for (i = 0; i <= 10; i++) x[i] = 0.2 * i;
            for (i = 0; i <= 20; i++) y[i] = 0.1 * i;
            for (i = 0; i <= 10; i++)
                for (j = 0; j <= 20; j++)
                    z[i, j] = Math.Exp(x[i] * x[i] - y[j] * y[j]);
            Heroius.XuAlgrithms.InterpolationApproximation.PIR2(x, y, z, 11, 21, out a, 6, 5, out dt);
            StringBuilder builder = new StringBuilder();
            builder.AppendLine("MATA(i,j) IS:");
            for (i = 0; i <= 5; i++)
            {
                for (j = 0; j <= 4; j++)
                    builder.AppendLine(a[i, j].ToString());
            }
            for (i = 0; i <= 2; i++)
                builder.AppendLine($"dt({i})={dt[i]}");
            return builder.ToString();
        }
    }
}
