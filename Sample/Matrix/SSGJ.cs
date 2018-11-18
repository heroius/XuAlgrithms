namespace Sample.Matrix
{
    public class SSGJ:SampleItem
    {
        public SSGJ() : base("Matrix.SSGJ") { }

        public override string Execute()
        {
            double[,] a = new double[,] {
                { 5,7,6,5 },
                { 7,10,8,7 },
                { 6,8,10,9 },
                { 5,7,9,10 }
            };
            Heroius.XuAlgrithms.Matrix.SSGJ(a);
            return Utility.MakeMatrixString(a);
        }
    }
}
