using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.Samples
{
    public class GINV:SampleItem
    {
        public GINV() : base("GINV") { }

        public override string Execute()
        {
            double[,] a = new double[,] {
                {1,2,3,4 },
                {6,7,8,9 },
                {1,2,13,0 },
                {16,17,8,9 },
                {2,4,3,4 }
            };
            double eps = 0.000001;
            Heroius.XuAlgrithms.Algrithms.GINV(a, eps, out double[,] aa, out double[,] u, out double[,] v);
            return Utility.MakeMatrixString(aa);
        }
    }
}
