using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.InterpolationApproximation
{
    public class PQS : SampleItem
    {
        public PQS() : base("InterpolationApproximation.PQS") { }

        public override string Execute()
        {
            double t, z;
            double[] x = new double[10] { 
                -1, -0.8, -0.65, -0.4, -0.3,
                0, 0.2, 0.45, 0.8, 1
            };
            double[] y = new double[10]
            {
                0.0384615, 0.0588236, 0.0864865, 0.2, 0.307692,
                1, 0.5, 0.164948, 0.0588236, 0.0384615
            };
            t = -0.85;
            z = Heroius.XuAlgrithms.InterpolationApproximation.PQS(x, y, 10, t);
            StringBuilder builder = new StringBuilder();
            builder.AppendLine($"t={t}, z={z}");
            t = 0.25;
            z = Heroius.XuAlgrithms.InterpolationApproximation.PQS(x, y, 10, t);
            builder.AppendLine($"t={t}, z={z}");
            return builder.ToString();
        }
    }
}
