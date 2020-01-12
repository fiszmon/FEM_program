using System;
using System.Linq;

namespace PROGRAM
{
    public class Element
    {
        public Int64 ID;
        public int[] Nodes;
        public const int N = 4;
        public MatrixFEM H, C, Hbc;
        public Double[] detJ, P;

        public bool[] AreaBC = Enumerable.Repeat(false, 4).ToArray();


        public Element(Int64 id, int[] points)
        {
            ID = id;
            Nodes = points;
            H = new MatrixFEM();
            Hbc = new MatrixFEM();
            C = new MatrixFEM();
            detJ = Enumerable.Repeat(0d, 4).ToArray();
            P = Enumerable.Repeat(0d, 4).ToArray();
        }

        public void ResetMatrixes()
        {
            H = new MatrixFEM();
            Hbc = new MatrixFEM();
            C = new MatrixFEM();
            detJ = Enumerable.Repeat(0d, 4).ToArray();
            P = Enumerable.Repeat(0d, 4).ToArray();
        }

        public void Add_Pvector(Double[] vP)
        {
            for (var i = 0; i < 4; i++)
                P[i] += vP[i];
        }
    }
}
