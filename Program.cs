using System;
using System.IO;
using System.Linq;
using System.Text;

namespace PROGRAM
{
    class Program
    {
        static void Main(string[] args)
        {
            Net _net = new Net(@"E:\STUDIA\MES\PROGRAM\PROGRAM\ProgramData.json");
            Grid _grid = _net.GetGrid();

            var matrixNormal = UniversalElementService.GenerateMatrixFEM_WithStDerative();
            var matrixXi = UniversalElementService.GenerateMatrixFEM_WithStDerative(UniversalDeriative.DerXi);
            var matrixEta = UniversalElementService.GenerateMatrixFEM_WithStDerative(UniversalDeriative.DerEta);

            var time = _net.GetGlobalData().Time;
            var step = _net.GetGlobalData().Step;
            var nE = _net.GetGlobalData().nE;

            var csv = new StringBuilder();

            //Console.WriteLine(NetService.PrintNodesTemperature(_net));

            for (var i = 0d; i < time; i += step)
            {

                for (var id = 0; id < nE; id++)
                {
                    var element = _grid.Elements[id];

                    var J = UniversalElementService.Count_J1_1(_grid, matrixEta, matrixXi, element);
                    element.detJ = UniversalElementService.Count_detJ(J, true);

                    var J1_1_1 = UniversalElementService.Count_J1_1_1(J, element.detJ, true);
                    var dN_dx = UniversalElementService.Count_dN_dx(matrixEta, matrixXi, true, J1_1_1);
                    var dN_dy = UniversalElementService.Count_dN_dy(matrixEta, matrixXi, true, J1_1_1);
                    var dictH = UniversalElementService.Count_matrixH(_net.GetGlobalData().Conductivity, dN_dx, dN_dy, element.detJ);
                    element.H = dictH[MatrixHVariables.H][IntegralPointsNames.All];

                    var dictC = UniversalElementService.Count_matrixC(matrixNormal, _net.GetGlobalData().Specific_Heat, _net.GetGlobalData().Density, element.detJ);
                    element.C = dictC[IntegralPointsNames.All];

                    var Hbc = new MatrixFEM();
                    for (var pow = 0; pow < 4; pow++)
                    {
                        if (element.AreaBC[pow])
                        {
                            var L = 0d;
                            if ((BCForHbc)pow == BCForHbc.Pow1 || (BCForHbc)pow == BCForHbc.Pow3)
                                L = _net.GetGlobalData().l_el_W;
                            else
                                L = _net.GetGlobalData().l_el_H;

                            var vP = UniversalElementService.Count_PforPow(L, _net.GetGlobalData().T_ambient, _net.GetGlobalData().Alpha, (BCForHbc)pow);
                            element.Add_Pvector(vP);
                            Hbc = Hbc + UniversalElementService.Count_Hbc_for_edge(_grid, element, _net.GetGlobalData().Alpha, (BCForHbc)pow);
                        }
                    }
                    element.Hbc = Hbc;

                    NetService.AgregateToGlobalWithElement(element);
                }

                //Console.WriteLine("------------  [C]  -----------------\n" + NetService.PrintGlobalC());
                //Console.WriteLine("\n\n------------  [H]  -----------------\n" + NetService.PrintGlobalH());

                NetService.Count_H_C_dT(_net.GetGlobalData().nN, step);
                
                //Console.WriteLine("\n\n------------  [H] + [C]/dT  -----------------\n" + NetService.PrintMatrix(H_CdT));
                var P_C_T0 = NetService.Count_P_CdT_T0(_net.GetGlobalData().nN, step, _grid);
                //Console.WriteLine("\n\n------------  P_Vector  -----------------\n" + NetService.PrintMatrix(P_C_T0));
                var T_final = NetService.Count_Final_vector(_net.GetGlobalData().nN);
                NetService.PrepareToNextStep(_net.GetGlobalData().nN, _grid, T_final);
                
                var max = NetService.FindMax(T_final);
                var min = NetService.FindMin(T_final);

                //in your loop
                var newLine = string.Format("{0},{1},{2}", i, min, max);
                csv.AppendLine(newLine);



                //Console.WriteLine(NetService.PrintNodesTemperature(_net));
                Console.WriteLine(String.Format("{2} {0} {1}", min, max, i));
            }

            //after your loop
            File.WriteAllText(@"E:\STUDIA\MES\PROGRAM\PROGRAM\ProgramResult" + DateTime.Now.Millisecond + ".csv", csv.ToString());
        }
    }
}
