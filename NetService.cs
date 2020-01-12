using System;

namespace PROGRAM
{
    public static class NetService
    {


        public static void InitializeGlogalMatrixes(Net net)
        {
            Net.globalH = new Double[net.GetGlobalData().nN][];
            Net.globalHbc = new Double[net.GetGlobalData().nN][];
            Net.globalC = new Double[net.GetGlobalData().nN][];
            Net.globalP = new Double[net.GetGlobalData().nN];
            for (var i = 0; i < net.GetGlobalData().nN; i++)
            {
                for (var j = 0; j < net.GetGlobalData().nN; j++)
                {
                    Net.globalH[i] = new Double[net.GetGlobalData().nN];
                    Net.globalHbc[i] = new Double[net.GetGlobalData().nN];
                    Net.globalC[i] = new Double[net.GetGlobalData().nN];
                    for (var k = 0; k < net.GetGlobalData().nN; k++)
                    {
                        Net.globalH[i][k] = 0f;
                        Net.globalHbc[i][k] = 0f;
                        Net.globalC[i][k] = 0f;
                    }
                }
                Net.globalP[i] = 0f;
            }
        }

        public static void AgregateToGlobalWithElement(Element el)
        {
            for (var i = 0; i < 4; i++)
            {
                Net.globalP[el.Nodes[i]] += el.P[i];
                for (var j = 0; j < 4; j++)
                {
                    Net.globalH[el.Nodes[i]][el.Nodes[j]] += el.H[i, j];
                    Net.globalHbc[el.Nodes[i]][el.Nodes[j]] += el.Hbc[i, j];
                    Net.globalC[el.Nodes[i]][el.Nodes[j]] += el.C[i, j];
                }
            }
        }

        private static Double[][] matrixH_C_dT = null;
        public static Double[][] Count_H_C_dT(int nN, Double step)
        {
            if (matrixH_C_dT != null)
                return matrixH_C_dT;

            matrixH_C_dT = new Double[nN][];
            for (var i = 0; i < nN; i++)
            {
                matrixH_C_dT[i] = new Double[nN];
                for (var j = 0; j < nN; j++)
                {
                    matrixH_C_dT[i][j] = Net.globalH[i][j] + Net.globalC[i][j]/step + Net.globalHbc[i][j];
                }
            }
            return matrixH_C_dT;
        }

        public static Double[] Count_P_CdT_T0(int nN, Double step, Grid grid)
        {
            for (int i = 0; i < nN; i++)
            {
                for (int j = 0; j < nN; j++)
                {
                    Net.globalC[i][j] /= step;
                    Net.globalC[j][i] *= grid.Nodes[i].T;
                }
            }
            for (int i = 0; i < nN; i++)
                for (int j = 0; j < nN; j++)
                    Net.globalP[i] += Net.globalC[i][j];

            return Net.globalP;
        }

        public static Double[] Count_Final_vector(int nN)
        {
            Double m, s, e;
            e = (Double)Math.Pow(10, -12);
            Double[] finalVec = new Double[nN];

            Double[][] tmp = new Double[nN][];
            for (int i = 0; i < nN; i++)
            {
                tmp[i] = new Double[nN + 1];
                for (int j = 0; j < nN; j++)
                {
                    tmp[i][j] = matrixH_C_dT[i][j];
                }
                tmp[i][nN] = Net.globalP[i];
            }

            for (int i = 0; i < nN - 1; i++)
            {
                for (int j = i + 1; j < nN; j++)
                {
                    if (Math.Abs(tmp[i][i]) < e)
                    {
                        Console.WriteLine("Can not divide by 0");
                        break;
                    }

                    m = -tmp[j][i] / tmp[i][i];
                    for (int k = 0; k < nN + 1; k++)
                    {
                        tmp[j][k] += m * tmp[i][k];
                    }
                }
            }

            for (int i = nN - 1; i >= 0; i--)
            {
                s = tmp[i][nN];
                for (int j = nN - 1; j >= 0; j--)
                {
                    s -= tmp[i][j] * finalVec[j];
                }
                if (Math.Abs(tmp[i][i]) < e)
                {
                    Console.WriteLine("Can not divide by 0");
                    break;
                }
                finalVec[i] = s / tmp[i][i];
            }
            return finalVec;
        }

        public static void PrepareToNextStep(int nN, Grid grid, Double[] T_next)
        {
            for (int i = 0; i < nN; i++)
            {
                grid.Nodes[i].T = T_next[i];
                for (int j = 0; j < nN; j++)
                {
                    Net.globalH[i][j] = 0;
                    Net.globalC[i][j] = 0;
                    Net.globalHbc[i][j] = 0;
                }
                Net.globalP[i] = 0f;
            }
            matrixH_C_dT = null;
            foreach (var el in grid.Elements)
                el.ResetMatrixes();
        }

        public static Double FindMax(Double[] tab)
        {
            var max = tab[0];
            foreach (var v in tab)
                if (max < v) max = v;
            return max;
        }

        public static Double FindMin(Double[] tab)
        {
            var min = tab[0];
            foreach (var v in tab)
                if (min > v) min = v;
            return min;
        }

        public static String PrintNodesTemperature(Net net)
        {
            var print = "Nodes Temperature:\n";
            for (int j = 0; j < net.GetGlobalData().nN; j++)
            {
                print+=net.GetGrid().Nodes[j].T + " ";
                if ((j+1)%net.GetGlobalData().nH == 0)
                    print += "\n";
            }
            print += "\n";
            return print;
        }

        public static String PrintGlobalH()
        {
            return PrintMatrix(Net.globalH);
        }

        public static String PrintGlobalHbc()
        {
            return PrintMatrix(Net.globalHbc);
        }

        public static String PrintGlobalC()
        {
            return PrintMatrix(Net.globalC);
        }

        public static String PrintMatrix(Double[][] matrix)
        {
            var print = "";
            foreach (var row in matrix)
            {
                foreach (var el in row)
                    print += String.Format("{0:0.000}", el) + " ";
                print += '\n';
            }
            return print;
        }

        public static String PrintMatrix(Double[] matrix)
        {
            var print = "";
            foreach (var el in matrix)
                print += String.Format("{0:0.000}", el) + " ";
            print += '\n';
            return print;
        }
    }
}
