using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Heroius.Extension;

namespace Sample.MatrixEigen
{
    public class SSTQ:SampleItem
    {
        public SSTQ() : base("MatrixEigen.SSTQ") { }

        public override string Execute()
        {
            double[] b = new double[]{ 9.2952, 11.6267, 10.9604, 6.11765, 15 };
            double[] c = new double[] {-0.749485, -4.49627, -2.15704, 7.14143,0 };
            double eps = 0.000001;
            int Max = 60;
            var q = Heroius.XuAlgrithms.MatrixEigen.SSTQ(b, c, eps, Max);
            return $"{b.Select(n => n.ToString()).Merge("\t")}\r\n{Utility.MakeMatrixString(q)}";
        }
    }
}
