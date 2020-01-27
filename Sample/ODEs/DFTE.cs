using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.ODEs
{
    public class DFTE : SampleItem
    {
        public DFTE() : base("ODEs.DFTE") { }

        public override string Execute()
        {
            double a = 2, b = 3, ya = 0, yb = 0;
            double[]
                y;
            Heroius.XuAlgrithms.ODEs.DFTE(a, b, ya, yb, 11, out y, Demo);
            StringBuilder builder = new StringBuilder();
            for (int i = 0; i < 11; i++)
            {
                builder.AppendLine($"y({i})={y[i]}");
            }
            return builder.ToString();
        }

        public double[] Demo(double x)
        {
            return new double[] {
                -1,0,2/(x*x),1/x
            };
        }
    }
}
