using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.Matrix
{
    public class TRCH:SampleItem
    {
        public TRCH() : base("Matrix.TRCH") { }

        public override string Execute()
        {
            double[] t = new double[] { 10, 5, 4, 3, 2, 1 };
            double[] tt = new double[] { 0, -1, -2, -3, -4, -5 };
            var b = Heroius.XuAlgrithms.Matrix.TRCH(t,tt);
            return Utility.MakeMatrixString(b);
        }
    }
}
