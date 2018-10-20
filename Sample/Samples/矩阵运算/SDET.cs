using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.Samples
{
    public class SDET:SampleItem
    {
        public SDET() : base("SDET") { }

        public override string Execute()
        {
            double[,] a = new double[,] {
                {1,2,3,4 },
                {5,6,7,8 },
                {9,10,11,12 },
                {13,14,15,16 }
            };
            double[,] b = new double[,] {
                {3,-3,-2,4 },
                {5,-5,1,8 },
                {11,8,5,-7 },
                {5,-1,-3,-1 }
            };
            var deta = Heroius.XuAlgrithms.Algrithms.SDET(a);
            var detb = Heroius.XuAlgrithms.Algrithms.SDET(b);
            return $"det(a)={deta}\r\ndet(b)={detb}";
        }
    }
}
