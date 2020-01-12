using System;

namespace PROGRAM
{
    public class Node
    {
        public Int64 ID;
        public Double X, Y, T;
        public BanderCondition BC;

        public Node(Int64 id, Double x, Double y, Double t, BanderCondition bander = BanderCondition.Normal)
        {
            ID = id;
            X = x;
            Y = y;
            T = t;
            BC = bander;
        }

        public string Print()
        {
            return X + " : " + Y + " / t: " + T + " / s: " + (BC == BanderCondition.Normal ? "0" : "1");
        }
    }
}
