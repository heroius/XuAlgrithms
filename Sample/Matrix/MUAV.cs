using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.Matrix
{
    public class MUAV:SampleItem
    {
        public MUAV() : base("Matrix.MUAV") { }

        public override string Execute()
        {
            double[,] a = new double[,] {
                { 1, 1, -1 },
                { 2, 1, 0 },
                { 1, -1, 0 },
                {-1, 2, 1 }
            };
            double[,] b = new double[,] {
                {1,1,-1,-1 },
                {2,1,0,2 },
                {1,-1,0,1 }
            };
            double eps = 0.000001;
            Heroius.XuAlgrithms.Matrix.MUAV(a, eps, out double[,] au, out double[,] av);
            Heroius.XuAlgrithms.Matrix.MUAV(b, eps, out double[,] bu, out double[,] bv);
            return $"matrix au:\r\n{Utility.MakeMatrixString(au)}\r\nmatrix av:\r\n{Utility.MakeMatrixString(av)}\r\nmatrix bu:\r\n{Utility.MakeMatrixString(bu)}\r\nmatrix bv:\r\n{Utility.MakeMatrixString(bv)}";
        }
    }
}
