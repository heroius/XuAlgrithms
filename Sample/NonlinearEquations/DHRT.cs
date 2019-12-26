using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.NonlinearEquations
{
    public class DHRT : SampleItem
    {
        public DHRT(): base("NonlinearEquations.DHRT") { }

        public override string Execute()
        {
            double[] x;
            var k = Heroius.XuAlgrithms.NonlinearEquations.DHRT(-2, 5, 0.2, 0.000001, 6, Demo, out x);
            StringBuilder sb = new StringBuilder();
            sb.Append($"M = {k}");
            for (int i = 0; i < k; i++)
            {
                sb.AppendLine();
                sb.Append($"x{i} = {x[i]}");
            }
            return sb.ToString();
        }

        public double Demo(double x)
        {
            return Math.Pow(x, 6) - 5* Math.Pow(x, 5) + 3*Math.Pow(x, 4) + Math.Pow(x, 3) - 7*Math.Pow(x, 2) + 7*x - 20;
        }
    }
}
