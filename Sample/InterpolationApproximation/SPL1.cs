using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.InterpolationApproximation
{
    public class SPL1 : SampleItem
    {
        public SPL1() : base("InterpolationApproximation.SPL1") { }

        public override string Execute()
        {
            int n = 12, m = 8;

            double[] x = new double[12] {
                0.52 , 8.0 , 17.95 , 28.65 , 50.65,104.6 ,
156.6, 260.7 , 364.4,468.0 , 507.0,520.0
            };
            double[] y = new double[12]
            {
5.28794,13.84,20.2,24.9,31.1,36.5,
36.6,31.0,20.9,7.8,1.5,0.2
            };
            double[] t = new double[8]
            {
                4.0,14.0,30.0,60.0,130.0,230.0,
450.0,515.0
            };

            double[] dy = new double[12], ddy, z, dz, ddz;
            dy[0] = 1.86548;
            dy[11] = -0.046115;
            double s =
            Heroius.XuAlgrithms.InterpolationApproximation.SPL1(x, y, n, ref dy, out ddy, t, m, out z, out dz, out ddz);
            StringBuilder builder = new StringBuilder();
            builder.AppendLine("x(i)    y(i)    dy(i)   ddy(i)");
            for (int i = 0; i < n; i++)
            {
                builder.AppendLine($"{x[i].ToString("E")} {y[i].ToString("E")}  {dy[i].ToString("E")}  {ddy[i].ToString("E")}");
            }
            builder.AppendLine($"s={s.ToString("E")}"); 
            builder.AppendLine("t(i)    z(i)    dz(i)   ddz(i)");
            for (int i = 0; i < m; i++)
            {
                builder.AppendLine($"{t[i].ToString("E")} {z[i].ToString("E")}  {dz[i].ToString("E")}  {ddz[i].ToString("E")}");
            }
            return builder.ToString();
        }
    }
}
