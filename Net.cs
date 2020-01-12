using Newtonsoft.Json.Linq;
using System;
using System.IO;

namespace PROGRAM
{
    public class Net
    {
        Grid grid;
        GlobalData globalData;

        public static Double[][] globalH = null;
        public static Double[][] globalHbc = null;
        public static Double[][] globalC = null;
        public static Double[] globalP = null;

        public Net(string path)
        {
            globalData = new GlobalData();
            if (globalData.ReadDataFromFile(path)) GenerateNet();
            SetAllEdgesInGrid(BanderCondition.Special);
            SetNodesInitialTemperature(globalData.T_init);
            NetService.InitializeGlogalMatrixes(this);
        }

        void GenerateNet()
        {
            grid = new Grid();
            grid.Elements = new Element[globalData.nE];
            grid.Nodes = new Node[globalData.nN];

            Double dX = globalData.H / (globalData.nH-1);
            Double dY = globalData.L / (globalData.nL-1);

            Double actualX = 0f, actualY = 0f;

            for (var iL = 0; iL < globalData.nL; iL++)
            {
                for (var iH = 0; iH < globalData.nH; iH++)
                {
                    //----------Node
                    Node tmpNode = new Node((Int64)iH + iL * globalData.nH, actualX, actualY, 21f);
                    grid.Nodes[iH + iL * globalData.nH] = tmpNode;
                    actualY += dX;

                    //----------Element
                    if (iH < globalData.nH - 1 && iL < globalData.nL - 1)
                    {
                        var i = iH + iL * globalData.nH;
                        var tmpPoints = new int[4];
                        tmpPoints[0] = i;
                        tmpPoints[1] = i + globalData.nH;
                        tmpPoints[2] = tmpPoints[1] + 1;
                        tmpPoints[3] = tmpPoints[0] + 1;
                        grid.Elements[i - iL] = new Element(i - iL, tmpPoints);
                    }
                }
                actualX += dY;
                actualY = 0f;
            }
        }

        public void SetAllEdgesInGrid(BanderCondition banderCondition)
        {
            for (var iL = 0; iL < globalData.nL; iL++)
                for (var iH = 0; iH < globalData.nH; iH++)
                    if (iL == 0 || iL == globalData.nL - 1 || iH == 0 || iH == (globalData.nH - 1))
                        grid.Nodes[iH + iL * globalData.nH].BC = banderCondition;
            SetAreaStateForElements();
        }

        private void SetAreaStateForElements()
        {
            foreach (var el in grid.Elements)
            {
                if (grid.Nodes[el.Nodes[0]].BC == BanderCondition.Special && grid.Nodes[el.Nodes[1]].BC == BanderCondition.Special)
                    el.AreaBC[(int)BCForHbc.Pow1] = true;
                if (grid.Nodes[el.Nodes[1]].BC == BanderCondition.Special && grid.Nodes[el.Nodes[2]].BC == BanderCondition.Special)
                    el.AreaBC[(int)BCForHbc.Pow2] = true;
                if (grid.Nodes[el.Nodes[2]].BC == BanderCondition.Special && grid.Nodes[el.Nodes[3]].BC == BanderCondition.Special)
                    el.AreaBC[(int)BCForHbc.Pow3] = true;
                if (grid.Nodes[el.Nodes[3]].BC == BanderCondition.Special && grid.Nodes[el.Nodes[0]].BC == BanderCondition.Special)
                    el.AreaBC[(int)BCForHbc.Pow4] = true;
            }
        }

        public void SetNodesInitialTemperature(Double T)
        {
            foreach (var node in grid.Nodes)
                node.T = T;
        }

        public void SetNodesTemperatures(Double[] T)
        {
            for (var i = 0; i < globalData.nE; i++)
                grid.Nodes[i].T = T[i];
        }

        public Grid GetGrid()
        {
            return grid;
        }
        public GlobalData GetGlobalData()
        {
            return globalData;
        }
    }

    public struct GlobalData
    {
        public Double H, L;
        public int nH, nL, nN, nE;
        public Double T_init, Time, Step, T_ambient,
            Alpha, Specific_Heat, Conductivity, Density;
        public Double l_el_H, l_el_W;

        public bool ReadDataFromFile(string path)
        {
            if (!File.Exists(path)) return false;

            var json = File.ReadAllText(path);
            JObject jsonData = JObject.Parse(json);

            H = (Double)jsonData.GetValue("H");
            L = (Double)jsonData.GetValue("L");
            nH = (int)jsonData.GetValue("nH");
            nL = (int)jsonData.GetValue("nL");
            T_init = (Double)jsonData.GetValue("T_init");
            Time = (Double)jsonData.GetValue("Time");
            Step = (Double)jsonData.GetValue("Step");
            T_ambient = (Double)jsonData.GetValue("T_ambient");
            Alpha = (Double)jsonData.GetValue("Alpha");
            Specific_Heat = (Double)jsonData.GetValue("Specific_heat");
            Conductivity = (Double)jsonData.GetValue("Conductivity");
            Density = (Double)jsonData.GetValue("Density");
            nN = nH * nL;
            nE = (nH - 1) * (nL - 1);
            l_el_H = H / (nH - 1f);
            l_el_W = L / (nL - 1f);
            return true;
        }
    }

    public struct Grid
    {
        public Element[] Elements;
        public Node[] Nodes;

        public String PrintElement(int index)
        {
            var el = Elements[index];
            return "Element " + (index + 1).ToString() +
                "\n" + (el.Nodes[0] + 1).ToString() + ". " + Nodes[el.Nodes[0]].Print() +
                "\n" + (el.Nodes[1] + 1).ToString() + ". " + Nodes[el.Nodes[1]].Print() +
                "\n" + (el.Nodes[2] + 1).ToString() + ". " + Nodes[el.Nodes[2]].Print() +
                "\n" + (el.Nodes[3] + 1).ToString() + ". " + Nodes[el.Nodes[3]].Print();

        }
    }

    public enum BanderCondition { Normal, Special };
}
