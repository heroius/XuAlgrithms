namespace Sample.Samples
{
    public class TRMUL : SampleItem
    {
        public TRMUL() : base("TRMUL") { }

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
            var c = Heroius.XuAlgrithms.Algrithms.TRMUL(a, b);
            return Utility.MakeMatrixString(c);
        }
    }
}
