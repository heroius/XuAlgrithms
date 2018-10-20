using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.Samples
{
    public class RANK:SampleItem
    {
        public RANK() : base("RANK") { }

        public override string Execute()
        {
            double[,] a = new double[,] {
                {1,2,3,4 },
                {5,6,7,8 },
                {9,10,11,12 },
                {13,14,15,16 },
                {17,18,19,20 }
            };
            return Heroius.XuAlgrithms.Algrithms.RANK(a).ToString();
        }
    }
}
