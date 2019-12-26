using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.NonlinearEquations
{
    public class PQRT:SampleItem
    {
        public PQRT() : base("NonlinearEquations.PQRT") { }

        public override string Execute()
        {
            double x = 1.0;
            var k = Heroius.XuAlgrithms.NonlinearEquations.PQRT(ref x, 0.0000001, Demo);
            return $"k = {k}, x = {x}";
        }

        public double Demo(double x)
        {
            return x*x * x-x*x-1;
        }
    }
}
