using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.NonlinearEquations
{
    public class ATKN:SampleItem
    {
        public ATKN() : base("NonlinearEquations.ATKN") { }

        public override string Execute()
        {
            double x = 0.0;
            var k = Heroius.XuAlgrithms.NonlinearEquations.ATKN(ref x, 0.0000001, 20, Demo);
            return $"k = {k}, x = {x}";
        }

        public double Demo(double x)
        {
            return 6 - x * x;
        }
    }
}
