using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.ODEs
{
    public class GEAR : SampleItem
    {
        public GEAR() : base("ODEs.GEAR") { }

        public override string Execute()
        {
            int i, j, k, jjs;
            double a, b, hmax, h;
            double[]
                y = new double[3],
                t = new double[30],
                hmin = new double[4] { 0.0001, 0.0001, 0.00001, 0.00001 },
                eps = new double[4] { 0.001, 0.0001, 0.0001, 0.00001 };
            double[,] z = new double[3, 30];
            a = 0; b = 1.0; h = 0.01; hmax = 0.1;
            StringBuilder builder = new StringBuilder();
            for (k = 0; k <= 3; k++)
            {
                y[0] = 1.0; y[1] = 0.0; y[2] = -1.0;
                jjs = Heroius.XuAlgrithms.ODEs.GEAR(a, b, hmin[k], hmax, h, eps[k], 3, y, 30, out t, out z, DemoS, DemoF);
                builder.AppendLine();
                builder.AppendLine($"h={h}");
                builder.AppendLine($"hmin={hmin[k]}");
                builder.AppendLine($"hmax={hmax}");
                builder.AppendLine($"eps={eps[k]}");
                builder.AppendLine($"jjs={jjs}");
                for (i = 0; i < 30; i++)
                {
                    builder.AppendLine($"t={t[i]}");
                    for (j = 0; j < 3; j++)
                        builder.Append($"y({j})={z[j, i]}\t\t");
                    builder.AppendLine();
                }
            }
            return builder.ToString();
        }

        public double[] DemoF(double t, double[] y, int n)
        {
            return new double[] {
                -21*y[0]+19*y[1]-20*y[2],
                19*y[0]-21*y[1]+20*y[2],
                40*y[0] -40*y[1]-40*y[2]
            };
        }

        public double[,] DemoS(double t, double[] y, int n)
        {
            double[,] p = new double[3, 3];
            p[0, 0] = -21.0; p[0, 1] = 19.0; p[0, 2] = -20.0;
            p[1, 0] = 19.0; p[1, 1] = -21.0; p[1, 2] = 20.0;
            p[2, 0] = 40.0; p[2, 1] = -40.0; p[2, 2] = -40.0;
            return p;
        }
    }
}
