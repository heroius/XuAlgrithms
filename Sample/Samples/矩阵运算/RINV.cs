using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.Samples
{
    public class RINV: SampleItem
    {
        public RINV() : base("RINV") { }

        public override string Execute()
        {
            double[,] a = new double[,] {
                { 0.2368, 0.2471, 0.2568, 1.2671 },
                { 1.1161, 0.1254, 0.1397, 0.1490 },
                { 0.1582, 1.1675, 0.1768, 0.1871 },
                { 0.1968, 0.2071, 1.2168, 0.2277 }
            };
            Heroius.XuAlgrithms.Algrithms.RINV(a);
            return Utility.MakeMatrixString(a);
        }
    }
}
