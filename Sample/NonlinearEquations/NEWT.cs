using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.NonlinearEquations
{
    public class NEWT : SampleItem
    {
        public NEWT(): base("NonlinearEquations.NEWT") { }

        public override string Execute()
        {
            double x = 1.5;
            var k = Heroius.XuAlgrithms.NonlinearEquations.NEWT(ref x, 0.000001, 60, Demo);
            return $"k = {k}, x = {x}";
        }

        public Tuple<double, double> Demo(double x)
        {
            var y = Math.Pow(x, 3) - Math.Pow(x, 2) - 1;
            var dy = 3 * x * x - 2 * x;
            return new Tuple<double, double>(y, dy);
        }
    }
}
