using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.Samples
{
    public class TRCH:SampleItem
    {
        public TRCH() : base("TRCH") { }

        public override string Execute()
        {
            double[] t = new double[] { 10, 5, 4, 3, 2, 1 };
            double[] tt = new double[] { 0, -1, -2, -3, -4, -5 };
            var b = Heroius.XuAlgrithms.Algrithms.TRCH(t,tt);
            return Utility.MakeMatrixString(b);
        }
    }
}
