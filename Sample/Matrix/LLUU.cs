using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.Matrix
{
    public class LLUU:SampleItem
    {
        public LLUU() : base("Matrix.LLUU") { }

        public override string Execute()
        {
            double[,] a = new double[,] {
                {2,4,4,2 },
                {3,3,12,6 },
                {2,4,-1,2 },
                {4,2,1,1 }
            };
            Heroius.XuAlgrithms.Matrix.LLUU(a, out double[,] l, out double[,] u);
            return $"lower triangle:\r\n{Utility.MakeMatrixString(l)}\r\nupper triangle:\r\n{Utility.MakeMatrixString(u)}";
        }
    }
}
