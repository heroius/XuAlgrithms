namespace Sample.Matrix
{
    public class TCMUL:SampleItem
    {
        public TCMUL() : base("Matrix.TCMUL") { }

        public override string Execute()
        {
            double[,] ar = new double[,]
            {
                { 1, 2, 3, -2 },
                { 1, 5, 1, 3 },
                { 0, 4, 2, -1 }
            };
            double[,] ai = new double[,] {
                { 1, -1, 2, 1 },
                { -1, -1, 2, 0 },
                { -3, -1, 2, 2 }
            };
            double[,] br = new double[,]
            {
                { 1, 4, 5, -2 },
                { 3, 0, 2, -1 },
                { 6, 3, 1, 2 },
                { 2, -3, -2, 1 }
            };
            double[,] bi = new double[,] {
                { -1, -1, 1, 1 },
                { 2, 1, 0, 5 },
                { -3, 2, 1, -1 },
                { -1, -2, 1, -2 }
            };
            Heroius.XuAlgrithms.Matrix.TCMUL(ar, ai, br, bi, out double[,] cr, out double[,] ci);
            return $"real part:\r\n{Utility.MakeMatrixString(cr)}imaginary part:\r\n{Utility.MakeMatrixString(ci)}";
        }
    }
}
