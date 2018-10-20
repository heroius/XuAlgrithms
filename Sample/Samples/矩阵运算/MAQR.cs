using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.Samples
{
    public class MAQR:SampleItem
    {
        public MAQR() : base("MAQR") { }

        public override string Execute()
        {
            double[,] a = new double[,] {
                {1,1,-1 },
                {2,1,0 },
                {1,-1,0 },
                {-1,2,1 }
            };
            var q = Heroius.XuAlgrithms.Algrithms.MAQR(a);
            return $"Q:\r\n{Utility.MakeMatrixString(q)}\r\nR:{Utility.MakeMatrixString(a)}";
        }
    }
}
