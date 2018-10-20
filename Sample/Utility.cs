using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample
{
    public static class Utility
    {
        public static string MakeMatrixString(double[,] m)
        {
            StringBuilder builder = new StringBuilder();
            int row = m.GetLength(0);
            int col = m.GetLength(1);
            for (int i = 0; i < row; i++)
            {
                for (int j = 0; j < col; j++)
                {
                    builder.Append(m[i, j]);
                    builder.Append("\t");
                }
                builder.Append("\r\n");
            }
            return builder.ToString();
        }
    }
}
