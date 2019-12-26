using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.NonlinearEquations
{
    public class SRRT : SampleItem
    {
        public SRRT() : base("NonlinearEquations.SRRT") { }

        public override string Execute()
        {
            double[] a = new double[] { -20.0, 7.0, -7.0, 1.0, 3.0, -5.0, 1.0 };
            double[] xr, xi;
            var i = Heroius.XuAlgrithms.NonlinearEquations.SRRT(a, 6, out xr, out xi);
            StringBuilder builder = new StringBuilder();
            if (i>0)
            {
                for (i = 0; i <= 5; i++)
                {
                    builder.AppendLine($"x({i}) = {xr[i].ToString("E")}, j {xi[i].ToString("E")}");
                }
            }
            return builder.ToString();
        }
    }
}
