using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.Samples
{
    public class CINV:SampleItem
    {
        public CINV() : base("CINV") { }

        public override string Execute()
        {
            double[,] ar = new double[,] {
                { 0.2368, 0.2471, 0.2568, 1.2671 },
                { 1.1161, 0.1254, 0.1397, 0.1490 },
                { 0.1582, 1.1675, 0.1768, 0.1871 },
                { 0.1968, 0.2071, 1.2168, 0.2277 }
            };
            double[,] ai = new double[,] {
                { 0.1345, 0.1678, 0.1875, 1.1161 },
                { 1.2671, 0.2017, 0.7024, 0.2721 },
                { -0.2836, -1.1967, 0.3558, -0.2078 },
                { 0.3576, -1.2345, 2.1185, 0.4773 }
            };
            Heroius.XuAlgrithms.Algrithms.CINV(ar, ai);
            return $"real part:\r\n{Utility.MakeMatrixString(ar)}imaginary part:\r\n{Utility.MakeMatrixString(ai)}";
        }
    }
}
