using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.InterpolationApproximation
{
    public class SPL : SampleItem
    {
        public SPL() : base("InterpolationApproximation.SPL") { }

        public override string Execute()
        {
            double t;
            double[] x = new double[11] { 
                -1, -0.95, -0.75, -0.55, -0.3, 0, 0.2, 0.45, 0.6, 0.8, 1
            };
            double[] y = new double[11]
            {
                0.0384615, 0.0424403, 0.06639, 0.116788,
                0.307692, 1, 0.5, 0.164948, 0.1, 0.0588236, 0.0384615
            };
            double[] s;
            int k = -1, n = 11;
            t = -0.85;
            Heroius.XuAlgrithms.InterpolationApproximation.SPL(x, y, n, k, t, out s);
            StringBuilder builder = new StringBuilder();
            builder.AppendLine($"x={t}, f(x)={s[4]}");
            builder.AppendLine($"s0={s[0]}, s1={s[1]}, s2={s[2]}, s3={s[3]}");
            t = 0.15;
            Heroius.XuAlgrithms.InterpolationApproximation.SPL(x, y, n, k, t, out s);
            builder.AppendLine($"x={t}, f(x)={s[4]}");
            builder.AppendLine($"s0={s[0]}, s1={s[1]}, s2={s[2]}, s3={s[3]}");
            return builder.ToString();
        }
    }
}
