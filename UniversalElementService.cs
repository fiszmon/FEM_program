using System;
using System.Collections.Generic;

namespace PROGRAM
{
    public enum UniversalDeriative { Without, DerXi, DerEta };
    public enum IntegralPointsNames { Pc1, Pc2, Pc3, Pc4, All };
    public enum MatrixHVariables { dN_dx_T, dN_dy_T, dN_dx_T_detJ, dN_dy_T_detJ, dN_dx_dy_detJ_K, H }
    public enum BCForHbc { Pow1, Pow2, Pow3, Pow4}

    public static class UniversalElementService
    {
        public static MatrixFEM GenerateMatrixFEM_WithStDerative(UniversalDeriative universalDeriative = UniversalDeriative.Without)
        {
            MatrixFEM matrix = new MatrixFEM();

            for (var i = 0; i < MatrixFEM.N; ++i)
            {
                matrix[0, i] = func_N1(UniversalElement.IntegralPoints[i].Xi, UniversalElement.IntegralPoints[i].Eta, universalDeriative);
                matrix[1, i] = func_N2(UniversalElement.IntegralPoints[i].Xi, UniversalElement.IntegralPoints[i].Eta, universalDeriative);
                matrix[2, i] = func_N3(UniversalElement.IntegralPoints[i].Xi, UniversalElement.IntegralPoints[i].Eta, universalDeriative);
                matrix[3, i] = func_N4(UniversalElement.IntegralPoints[i].Xi, UniversalElement.IntegralPoints[i].Eta, universalDeriative);
            }

            return matrix;
        }

        public static Double func_N1(Double xi, Double eta, UniversalDeriative universalDeriative = UniversalDeriative.Without)
        {
            Double value = 0;
            switch (universalDeriative)
            {
                case UniversalDeriative.Without:
                    value = 0.25f * (1 - xi) * (1 - eta);
                    break;
                case UniversalDeriative.DerEta:
                    value = -0.25f * (1 - xi);
                    break;
                case UniversalDeriative.DerXi:
                    value = -0.25f * (1 - eta);
                    break;
            }

            return value;
        }
        public static Double func_N2(Double xi, Double eta, UniversalDeriative universalDeriative = UniversalDeriative.Without)
        {
            Double value = 0;
            switch (universalDeriative)
            {
                case UniversalDeriative.Without:
                    value = 0.25f * (1 + xi) * (1 - eta);
                    break;
                case UniversalDeriative.DerEta:
                    value = -0.25f * (1 + xi);
                    break;
                case UniversalDeriative.DerXi:
                    value = 0.25f * (1 - eta);
                    break;
            }

            return value;
        }
        static Double func_N3(Double xi, Double eta, UniversalDeriative universalDeriative = UniversalDeriative.Without)
        {
            Double value = 0;
            switch (universalDeriative)
            {
                case UniversalDeriative.Without:
                    value = 0.25f * (1 + xi) * (1 + eta);
                    break;
                case UniversalDeriative.DerEta:
                    value = 0.25f * (1 + xi);
                    break;
                case UniversalDeriative.DerXi:
                    value = 0.25f * (1 + eta);
                    break;
            }

            return value;
        }
        static Double func_N4(Double xi, Double eta, UniversalDeriative universalDeriative = UniversalDeriative.Without)
        {
            Double value = 0;
            switch (universalDeriative)
            {
                case UniversalDeriative.Without:
                    value = 0.25f * (1 - xi) * (1 + eta);
                    break;
                case UniversalDeriative.DerEta:
                    value = 0.25f * (1 - xi);
                    break;
                case UniversalDeriative.DerXi:
                    value = -0.25f * (1 + eta);
                    break;
            }

            return value;
        }

        public static MatrixFEM Count_J1_1(Grid grid, MatrixFEM NderEta, MatrixFEM NderXi, Element element) // obliczanie Jakobianu przekształcenia
        {
            Node[] nodes = GetElementNodes(grid, element);
            MatrixFEM J = new MatrixFEM();
            for (var i = 0; i < nodes.Length; i++)
            {
                J[0, i] = NderXi[0, 0] * nodes[0].X + NderXi[1, 0] * nodes[1].X + NderXi[2, 0] * nodes[2].X + NderXi[3, 0] * nodes[3].X;
                J[1, i] = NderXi[0, 0] * nodes[0].Y + NderXi[1, 0] * nodes[1].Y + NderXi[2, 0] * nodes[2].Y + NderXi[3, 0] * nodes[3].Y;
                J[2, i] = NderEta[0, 0] * nodes[0].X + NderEta[1, 0] * nodes[1].X + NderEta[2, 0] * nodes[2].X + NderEta[3, 0] * nodes[3].X;
                J[3, i] = NderEta[0, 0] * nodes[0].Y + NderEta[1, 0] * nodes[1].Y + NderEta[2, 0] * nodes[2].Y + NderEta[3, 0] * nodes[3].Y;
            }
            return J;
        }

        private static Double[] detJ = null;
        public static Double[] Count_detJ(MatrixFEM J, bool count = false) // wyznacznik macierzy Jakobianu przekształcenia J
        {
            if (detJ != null && !count)
                return detJ;

            detJ = new Double[4];
            for (var i = 0; i < 4; i++)
                detJ[i] = J[0, i] * J[3, i] - J[i, 1] * J[2, i];
            return detJ;
        }

        private static MatrixFEM J1_1_1 = null;
        public static MatrixFEM Count_J1_1_1(MatrixFEM J, Double[] detJ, bool count = false) //podzielenie macierzy Jakobianu J przez tablice wyznaczników detJ
        {
            if (J1_1_1 != null && !count)
                return J1_1_1;

            J1_1_1 = new MatrixFEM();
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    switch (i)
                    {
                        case 0:
                            J1_1_1[i, j] = J[3, i] / detJ[i];
                            break;
                        case 1:
                            J1_1_1[i, j] = J[1, i] / detJ[i];
                            break;
                        case 2:
                            J1_1_1[i, j] = J[2, i] / detJ[i];
                            break;
                        case 3:
                            J1_1_1[i, j] = J[0, i] / detJ[i];
                            break;
                    }
                }
            }
            return J1_1_1;
        }

        private static MatrixFEM dN_dx = null;
        public static MatrixFEM Count_dN_dx(MatrixFEM NderEta, MatrixFEM NderXi, bool count = false, MatrixFEM J = null) // macierz pochodnych funkcji kształtu po dx
        {
            if (dN_dx != null && !count)
                return dN_dx;

            if (J == null)
                J = J1_1_1;
            dN_dx = new MatrixFEM();
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    dN_dx[i, j] = J[0, i] * NderXi[i, j] + J[1, i] * NderEta[i, j];

            return dN_dx;
        }

        private static MatrixFEM dN_dy = null;
        public static MatrixFEM Count_dN_dy(MatrixFEM NderEta, MatrixFEM NderXi, bool count = false, MatrixFEM J = null) // macierz pochodnych funkcji kształtu po dx
        {
            if (dN_dy != null && !count)
                return dN_dy;

            if (J == null)
                J = J1_1_1;
            dN_dy = new MatrixFEM();
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    dN_dy[i, j] = J[2, i] * NderXi[i, j] + J[3, i] * NderEta[i, j];

            return dN_dy;
        }

        public static Dictionary<MatrixHVariables, Dictionary<IntegralPointsNames, MatrixFEM>> Count_matrixH
            (Double K, MatrixFEM dN_dx = null, MatrixFEM dN_dy = null, Double[] detJ = null)
        {
            if (dN_dx == null || dN_dy == null || detJ == null)
            {
                dN_dx = UniversalElementService.dN_dx;
                dN_dy = UniversalElementService.dN_dy;
                detJ = UniversalElementService.detJ;
            }

            //MatrixFEM dN_dxT = dN_dx.GetTransposition(), dN_dyT = dN_dy.GetTransposition();

            Dictionary<IntegralPointsNames, MatrixFEM> dictdNdx_dNdxT = new Dictionary<IntegralPointsNames, MatrixFEM>();
            Dictionary<IntegralPointsNames, MatrixFEM> dictdNdy_dNdyT = new Dictionary<IntegralPointsNames, MatrixFEM>();
            Dictionary<IntegralPointsNames, MatrixFEM> dictdNdx_dNdxT_detJ = new Dictionary<IntegralPointsNames, MatrixFEM>();
            Dictionary<IntegralPointsNames, MatrixFEM> dictdNdy_dNdyT_detJ = new Dictionary<IntegralPointsNames, MatrixFEM>();
            Dictionary<IntegralPointsNames, MatrixFEM> dictdNdx_dNdxT_dNdy_dNdyT_detJ_K = new Dictionary<IntegralPointsNames, MatrixFEM>();
            var dMatrixH = new Dictionary<IntegralPointsNames, MatrixFEM>();
            dMatrixH[IntegralPointsNames.All] = new MatrixFEM();

            for (var i = 0; i < 4; i++)
            {
                dictdNdx_dNdxT[(IntegralPointsNames)i] = new MatrixFEM();
                dictdNdy_dNdyT[(IntegralPointsNames)i] = new MatrixFEM();
                dictdNdx_dNdxT_detJ[(IntegralPointsNames)i] = new MatrixFEM();
                dictdNdy_dNdyT_detJ[(IntegralPointsNames)i] = new MatrixFEM();
                dictdNdx_dNdxT_dNdy_dNdyT_detJ_K[(IntegralPointsNames)i] = new MatrixFEM();
            }

            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                {
                    for (var k = 0; k < 4; k++)
                    {
                        dictdNdx_dNdxT[(IntegralPointsNames)k][i, j] = dN_dx[i, k] * dN_dx[j, k];
                        dictdNdy_dNdyT[(IntegralPointsNames)k][i, j] = dN_dy[i, k] * dN_dy[j, k];

                        dictdNdx_dNdxT_detJ[(IntegralPointsNames)k][i, j] = dictdNdx_dNdxT[(IntegralPointsNames)k][i, j] * detJ[k];
                        dictdNdy_dNdyT_detJ[(IntegralPointsNames)k][i, j] = dictdNdy_dNdyT[(IntegralPointsNames)k][i, j] * detJ[k];

                        dictdNdx_dNdxT_dNdy_dNdyT_detJ_K[(IntegralPointsNames)k][i, j] =
                            K * (dictdNdx_dNdxT_detJ[(IntegralPointsNames)k][i, j] + dictdNdy_dNdyT_detJ[(IntegralPointsNames)k][i, j]);

                        dMatrixH[IntegralPointsNames.All][i, j] += dictdNdx_dNdxT_dNdy_dNdyT_detJ_K[(IntegralPointsNames)k][i, j];
                    }
                }

            var dictH = new Dictionary<MatrixHVariables, Dictionary<IntegralPointsNames, MatrixFEM>>();

            dictH[MatrixHVariables.dN_dx_T] = dictdNdx_dNdxT;
            dictH[MatrixHVariables.dN_dy_T] = dictdNdy_dNdyT;
            dictH[MatrixHVariables.dN_dx_T_detJ] = dictdNdx_dNdxT_detJ;
            dictH[MatrixHVariables.dN_dy_T_detJ] = dictdNdy_dNdyT_detJ;
            dictH[MatrixHVariables.dN_dx_dy_detJ_K] = dictdNdx_dNdxT_dNdy_dNdyT_detJ_K;
            dictH[MatrixHVariables.H] = dMatrixH;

            return dictH;
        }

        public static Dictionary<IntegralPointsNames, MatrixFEM> Count_matrixC(MatrixFEM matrixN, Double c, Double ro, Double[] detJ = null)
        {
            if (detJ == null)
                detJ = UniversalElementService.detJ;

            var dictN = new Dictionary<IntegralPointsNames, MatrixFEM>();
            dictN[IntegralPointsNames.All] = new MatrixFEM();
            for (var k = 0; k < 4; k++)
            {
                dictN[(IntegralPointsNames)k] = new MatrixFEM();
                for (var i = 0; i < 4; i++)
                    for (var j = 0; j < 4; j++)
                    {
                        dictN[(IntegralPointsNames)k][i, j] = matrixN[i, k] * matrixN[j, k] * detJ[k] * c * ro;
                        dictN[IntegralPointsNames.All][i, j] += dictN[(IntegralPointsNames)k][i, j];
                    }
            }

            return dictN;
        }

        public static MatrixFEM Count_Hbc_for_edge(Grid grid, Element element, Double alpha, BCForHbc pow)
        {
            MatrixFEM matrixHbc = new MatrixFEM();
            Node[] nodes = GetElementNodes(grid, element);
            Double[][] matrixN = new Double[2][];
            var points = UniversalElement.GetBCPoints(pow);
            Double L = 0;

            switch (pow)
            {
                case BCForHbc.Pow1:
                    L = (Double)Math.Sqrt(Math.Pow(nodes[1].X - nodes[0].X, 2) + Math.Pow(nodes[1].Y - nodes[0].Y, 2));
                    break;
                case BCForHbc.Pow2:
                    L = (Double)Math.Sqrt(Math.Pow(nodes[2].X - nodes[1].X, 2) + Math.Pow(nodes[2].Y - nodes[1].Y, 2));
                    break;
                case BCForHbc.Pow3:
                    L = (Double)Math.Sqrt(Math.Pow(nodes[3].X - nodes[2].X, 2) + Math.Pow(nodes[3].Y - nodes[2].Y, 2));
                    break;
                case BCForHbc.Pow4:
                    L = (Double)Math.Sqrt(Math.Pow(nodes[0].X - nodes[3].X, 2) + Math.Pow(nodes[0].Y - nodes[3].Y, 2));
                    break;
            }
            Double detJ = L / 2f;
            for (var i = 0; i < 2; i++)
            {
                matrixN[i] = new Double[4];
                matrixN[i][0] = func_N1(points[i].Xi, points[i].Eta);
                matrixN[i][1] = func_N2(points[i].Xi, points[i].Eta);
                matrixN[i][2] = func_N3(points[i].Xi, points[i].Eta);
                matrixN[i][3] = func_N4(points[i].Xi, points[i].Eta);
            }

            for (var i = 0; i < 4; i++)
                for (var j = 0; j < 4; j++)
                {
                    matrixHbc[i, j] = matrixN[0][j] * matrixN[0][i];
                    matrixHbc[i, j] += matrixN[1][j] * matrixN[1][i];
                    matrixHbc[i, j] *= alpha * detJ;
                }

            return matrixHbc;
        }

        static Dictionary<BCForHbc, Double[]> dictvP = new Dictionary<BCForHbc, Double[]>();

        public static Double[] Count_PforPow(Double L, Double T0, Double alpha, BCForHbc pow, bool count = false)
        {
            Double[] vP;
            if (!dictvP.TryGetValue(pow, out vP) || count)
            {
                dictvP[pow] = Count_vP(UniversalElement.GetBCPoints(pow), L, T0, alpha);
                vP = dictvP[pow];
            }
            return vP;
        }

        private static Double[] Count_vP(Point[] points, Double L, Double T0, Double alpha)
        {
            var vp = new Double[4];
            vp[0] = (func_N1(points[0].Xi, points[0].Eta) + func_N1(points[1].Xi, points[1].Eta)) * L / 2f * T0 * alpha;
            vp[1] = (func_N2(points[0].Xi, points[0].Eta) + func_N2(points[1].Xi, points[1].Eta)) * L / 2f * T0 * alpha;
            vp[2] = (func_N3(points[0].Xi, points[0].Eta) + func_N3(points[1].Xi, points[1].Eta)) * L / 2f * T0 * alpha;
            vp[3] = (func_N4(points[0].Xi, points[0].Eta) + func_N4(points[1].Xi, points[1].Eta)) * L / 2f * T0 * alpha;
            return vp;
        }

        private static Node[] GetElementNodes(Grid grid, Element element)
        {
            var nodes = new Node[element.Nodes.Length];
            for (var i = 0; i < element.Nodes.Length; ++i)
            {
                nodes[i] = grid.Nodes[element.Nodes[i]];
            }
            return nodes;
        }
    }
}
