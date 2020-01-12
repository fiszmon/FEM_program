using System;
using System.Collections.Generic;

namespace PROGRAM
{
    public static class UniversalElement
    {
        static Point[] tabIntegralPoints = null;
        static Point[] tabPoints = null;
        private static Dictionary<BCForHbc, Point[]> dictBCpoints = new Dictionary<BCForHbc, Point[]>()
        {
            { BCForHbc.Pow1,
                    new Point[] {new Point(-1 / Math.Sqrt(3), -1f), new Point(1 / Math.Sqrt(3), -1) } },
            { BCForHbc.Pow2,
                    new Point[] {new Point(1, -1 / Math.Sqrt(3)), new Point(1 , 1 / Math.Sqrt(3)) } },
            { BCForHbc.Pow3,
                    new Point[] {new Point(1 / Math.Sqrt(3), 1), new Point(-1 / Math.Sqrt(3), 1) } },
            { BCForHbc.Pow4,
                    new Point[] {new Point(-1 , 1 / Math.Sqrt(3)), new Point(-1, -1 / Math.Sqrt(3)) } }
        };

        public static Point[] IntegralPoints
        {
            get
            {
                if (tabIntegralPoints is null)
                    tabIntegralPoints = LocalGaussPoints2d.Points;
                return tabIntegralPoints;
            }
            set
            {
                tabIntegralPoints = value;
            }
        }
        public static Point[] Points
        {
            get
            {
                if (tabPoints is null)
                    tabPoints = generateLocalPoints();
                return tabPoints;
            }
            set
            {
                tabPoints = value;
            }
        }

        public static Point[] GetBCPoints(BCForHbc pow)
        {
            return dictBCpoints[pow];
        }

        static Point[] generateLocalPoints()
        {
            var points = new Point[4];
            points[0] = new Point(-1f, -1f, 0);
            points[1] = new Point(1f, -1f, 0);
            points[2] = new Point(1f, 1f, 0);
            points[3] = new Point(-1f, 1f, 0);

            return points;
        }

        public static class LocalGaussPoints2d
        {
            public static Double LPc = (Double)(1 / Math.Sqrt(3));
            static Point Pc1 = new Point(-LPc, -LPc);
            static Point Pc2 = new Point(LPc, -LPc);
            static Point Pc3 = new Point(LPc, LPc);
            static Point Pc4 = new Point(-LPc, LPc);

            public static Point[] Points = { Pc1, Pc2, Pc3, Pc4 };
        }
    }
}
