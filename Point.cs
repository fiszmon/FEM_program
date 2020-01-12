using System;

namespace PROGRAM
{
    public class Point
    {
        public Double Xi;
        public Double Eta;
        public Double Weight;

        public Point(Double xi, Double eta, Double weight = 1f)
        {
            Xi = xi;
            Eta = eta;
            Weight = weight;
        }
    }
}
