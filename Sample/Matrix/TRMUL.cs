namespace Sample.Matrix
{
    public class TRMUL : SampleItem
    {
        public TRMUL() : base("Matrix.TRMUL") { }

        public override string Execute()
        {
            double[,] a = new double[,] {
                { 1,3,-2,0,4 },
                {-2,-1,5,-7,2 },
                {0,8,4,1,-5 },
                {3,-3,2,-4,1 }
            };
            double[,] b = new double[,] {
                {4,5,-1 },
                {2,-2,6 },
                {7,8,1 },
                {0,3,-5 },
                {9,8,-6 }
            };
            Heroius.XuAlgrithms.Matrix.TRMUL(a, b, 4, 5, 3, out double[,] c);
            return Utility.MakeMatrixString(c);
        }
    }
}
