namespace Sample.Samples
{
    public class SSGJ:SampleItem
    {
        public SSGJ() : base("SSGJ") { }

        public override string Execute()
        {
            double[,] a = new double[,] {
                { 5,7,6,5 },
                { 7,10,8,7 },
                { 6,8,10,9 },
                { 5,7,9,10 }
            };
            Heroius.XuAlgrithms.Algrithms.SSGJ(a);
            return Utility.MakeMatrixString(a);
        }
    }
}
