using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.InterpolationApproximation
{
    public class SPL3 : SampleItem
    {
        public SPL3() : base("InterpolationApproximation.SPL3") { }

        public override string Execute()
        {
            string fm = "0.000000";
            int n, m, i, j;
            double u, s;
            double[] x = new double[37], y = new double[37], dy, ddy;
            double[] t = new double[36], z, dz, ddz;
            for (i = 0; i <= 36; i++)
            {
                x[i] = i * 6.2831852 / 36.0;
                y[i] = Math.Sin(x[i]);
            }
            for (i = 0; i <= 35; i++)
                t[i] = (0.5 + i) * 6.2831852 / 36.0;
            n = 37; m = 36;
            s = Heroius.XuAlgrithms.InterpolationApproximation.SPL3(x, y, n, out dy, out ddy, t, m, out z, out dz, out ddz);
            StringBuilder builder = new StringBuilder();
            builder.AppendLine("x(i)\t\ty(i)=sin(x)\t\tdy(i)=cos(x)\t\tddy(i)=-sin(x)");
            builder.AppendLine($"{x[0]}\t\t{y[0].ToString(fm)}\t\t{dy[0].ToString(fm)}\t\t{ddy[0].ToString(fm)}");
            for (i = 0; i <= 35; i++)
            {
                u = t[i] * 36.0 / 0.62831852;
                builder.AppendLine($"{u}\t\t{z[i].ToString(fm)}\t\t{dz[i].ToString(fm)}\t\t{ddz[i].ToString(fm)}");
                u = x[i + 1] * 36.0 / 0.62831852;
                j = i + 1;
                builder.AppendLine($"{u}\t\t{y[j].ToString(fm)}\t\t{dy[j].ToString(fm)}\t\t{ddy[j].ToString(fm)}");
            }
            builder.AppendLine($"s={s.ToString("E")}");
            return builder.ToString();
        }
    }
}
