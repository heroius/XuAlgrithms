using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.NonlinearEquations
{
    public class CSRT : SampleItem
    {
        public CSRT() : base("NonlinearEquations.CSRT") { }

        public override string Execute()
        {
            double[] ar = new double[] { 0.1, 21.33, 4.9,0.0,3.0,1.0 };
            double[] ai = new double[] {-100.0,0.0,-19.0,-0.01,2.0,0.0 };
            double[] xr, xi;
            var i = Heroius.XuAlgrithms.NonlinearEquations.CSRT(ar, ai, 5, out xr, out xi);
            StringBuilder builder = new StringBuilder();
            if (i>0)
            {
                for (i = 0; i <= 4; i++)
                {
                    builder.AppendLine($"x({i}) = {xr[i].ToString("E")}, j {xi[i].ToString("E")}");
                }
            }
            return builder.ToString();
        }
    }
}
