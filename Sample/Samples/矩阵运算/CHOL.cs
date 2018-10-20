using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.Samples
{
    public class CHOL:SampleItem
    {
        public CHOL() : base("CHOL") { }

        public override string Execute()
        {
            double[,] a = new double[,] {
                {5,7,6,5 },
                { 7,10,8,7 },
                {6,8,10,9 },
                {5,7,9,10 }
            };
            return Heroius.XuAlgrithms.Algrithms.CHOL(a).ToString();
        }
    }
}
