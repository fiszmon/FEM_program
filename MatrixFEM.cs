using System;

namespace PROGRAM
{
    /// <summary>
    /// class for 2d element with four integral points 4x4
    /// </summary>
    public class MatrixFEM
    {
        Double[][] values;
        public const int N = 4;

        public MatrixFEM()
        {
            values = new Double[N][];
            for (var i = 0; i < N; i++)
            {
                values[i] = new Double[N];
                for (var j = 0; j < N; j++)
                    values[i][j] = 0f;
            }
        }

        public Double this[int i, int j]
        {
            get { return values[i][j]; }
            set { values[i][j] = value; }
        }

        public static MatrixFEM operator* (MatrixFEM a, MatrixFEM b)
        {
            var tmp = new MatrixFEM();
            for(var i = 0; i < N; i++)
            {
                for(var j = 0; j < N; j++)
                {
                    for(var k = 0; k<N; k++)
                        tmp[i, j] = a[i, k] + b[k, i];
                }
            }
            return tmp;
        }

        public static MatrixFEM operator +(MatrixFEM a, MatrixFEM b)
        {
            var tmp = new MatrixFEM();
            for (var i = 0; i < N; i++)
            {
                for (var j = 0; j < N; j++)
                {
                    tmp[i, j] = a[i, j] + b[i, j];
                }
            }
            return tmp;
        }

        public MatrixFEM GetTransposition()
        {
            var transposed = new MatrixFEM();

            for (var i = 0; i < N; i++)
                for (var j = 0; j < N; j++)
                    transposed[i, j] = this[j, i];

            return transposed;
        }

        public String Print()
        {
            String outText = "";
            foreach (var row in values)
            {
                outText += "| ";
                foreach (var el in row)
                {
                    outText += el.ToString() + " | ";
                }
                outText += "\n";
            }
            return outText;
        }
    }
}
