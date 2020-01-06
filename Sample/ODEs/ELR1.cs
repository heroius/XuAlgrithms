using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.ODEs
{
    public class ELR1 : SampleItem
    {
        public ELR1() : base("ODEs.ELR1") { }

        public override string Execute()
        {
            int n = 3, k = 11;
            double[] y = new double[3] { -1, 0, 1 }; double[,] z = new double[3,11];
            double t=0,h=0.01,x;
            Heroius.XuAlgrithms.ODEs.ELR1(t, y, 3, h, 11, out z, Demo);
            StringBuilder builder = new StringBuilder();
            for (int i = 0; i <= 10; i++)
            {
                x = i * h;
                builder.AppendLine($"t={x}");
                for (int j = 0; j <= 2; j++)
                    builder.Append($"y({j})={z[j,i].ToString("E")}    ");
                builder.AppendLine();
            }
            return builder.ToString();
        }

        public double[] Demo(double t, double[] y, int n)
        {
            return new double[3] { 
                y[1],
                -y[0],
                -y[2]
            };
        }
    }
}
