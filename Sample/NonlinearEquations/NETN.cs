using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.NonlinearEquations
{
    public class NETN : SampleItem
    {
        public NETN() : base("NonlinearEquations.NETN") { }

        public override string Execute()
        {
            double[] x = new double[] { 1, 1, 1 };
            double t = 0.1, h = 0.1, e = 0.0000001;
            int k = 100;
            var i = Heroius.XuAlgrithms.NonlinearEquations.NETN(3, e, t, h, x, k, Demo);
            StringBuilder builder = new StringBuilder();
            builder.AppendLine($"i={i}");
            for (i = 0; i <= 2; i++)
            {
                builder.AppendLine($"x({i}) = {x[i].ToString("E")}");
            }
            return builder.ToString();
        }


        public void Demo(double[] x, double[] y, int n)
        {
            y[0] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - 1;
            y[1] = 2.0 * x[0] * x[0] + x[1] * x[1] - 4 * x[2];
            y[2] = 3.0 * x[0] * x[0] - 4 * x[1] + x[2]*x[2];
        }
    }
}
