using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample.DataProcessing
{
    public class RHIS : SampleItem
    {
        public RHIS() : base("DataProcessing.RHIS") { }

        public override string Execute()
        {
            int n = 100, m = 10;
            int[] g, q;
            double[] dt = new double[3];
            double x0=192, h=2;
            double[] x = new double[] { 
                193.199,195.673,195.757,196.051,196.092,196.596,
                196.579,196.763,196.847,197.267,197.392,197.477,
                198.189,193.850,198.944,199.070,199.111,199.153,
                199.237,199.698,199.572,199.614,199.824,199.908,
                200.188,200.160,200.243,200.285,200.453,200.704,
                200.746,200.830,200.872,200.914,200.956,200.998,
                200.998,201.123,201.208,201.333,201.375,201.543,
                201.543,201.584,201.711,201.878,201.919,202.004,
                202.004,202.088,202.172,202.172,202.297,202.339,
                202.381,202.507,202.591,202.716,202.633,202.884,
                203.051,203.052,203.094,203.094,203.177,203.178,
                203.219,203.764,203.765,203.848,203.890,203.974,
                204.184,204.267,204.352,204.352,204.729,205.106,
                205.148,205.231,205.357,205.400,205.483,206.070,
                206.112,206.154,206.155,206.615,206.657,206.993,
                207.243,207.621,208.124,208.375,208.502,208.628,
                208.670,208.711,210.012,211.394
            };
            Heroius.XuAlgrithms.DataProcessing.RHIS(x, n, x0, h, m, out dt, out g, out q);

            StringBuilder builder = new StringBuilder();
            builder.AppendLine($"n={n}");
            builder.AppendLine($"x0={x0}\t\th={h}\t\tm={m}");
            builder.AppendLine($"xa={dt[0]}\t\ts={dt[1]}\t\tt={dt[2]}");

            int k = 1, z = 0, j;
            double s;
            char[] a = new char[50];
            for (int i = 0; i < m; i++)
            {
                if (q[i]>z)
                {
                    z = q[i];
                }
            }
            while (z>50)
            {
                z = z / 2;
                k = 2 * k;
            }
            for (int i = 0; i < m; i++)
            {
                s= x0 + (i + 0.5) * h;
                for (j = 0; j < 50; j++)
                {
                    a[j] = ' ';
                }
                j = q[i] / k;
                for (z = 0; z < j; z++)
                {
                    a[z] = 'X';
                }
                j = g[i] / k;
                if ((j > 0) && (j < 50)) a[j] = '*';
                builder.Append($"{s}\t\t{ q[i]}\t\t");
                for (j = 0; j < 50; j++)
                {
                    builder.Append(a[j]);
                }
                builder.AppendLine();
            }
            return builder.ToString();
        }
    }
}
