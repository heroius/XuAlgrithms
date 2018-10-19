using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.Samples
{
    public class TRMUL : SampleItem
    {
        public TRMUL() : base("TRMUL") { }

        public override string Execute()
        {
            double[,] a = new double[,] {
                { 0.01,3,-2,0,4 },
                {-2,-1,5,-7,2 },
                {0,8,4,1,-5 },
                {3,-3,0.002,-4,1 }
            };
            double[,] b = new double[,] {
                {4,5,-1 },
                {2,-2,6 },
                {7,8,1 },
                {0,3,-5 },
                {9,8,-6 }
            };
            Heroius.XuAlgrithms.TRMUL core = new Heroius.XuAlgrithms.TRMUL(4, 5, 3);
            core.Input(a, b);
            core.MUL();
            var c = core.Output();
            StringBuilder builder = new StringBuilder();
            int row = c.GetLength(0);
            int col = c.GetLength(1);
            for (int i = 0; i < row; i++)
            {
                for (int j = 0; j < col; j++)
                {
                    builder.Append(c[i, j]);
                    builder.Append("\t");
                }
                builder.Append("\r\n");
            }
            return builder.ToString();
        }
    }
}
