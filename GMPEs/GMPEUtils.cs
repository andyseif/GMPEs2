using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Hazard;

namespace GMPEs
{

    public enum SiteType
    {
        FIRM_ROCK,
        HARD_ROCK
    }

    public enum MagnitudeType
    {
        MOMENT,
        LG_PHASE
    }

    public enum FaultStyle
    {
        NORMAL,
        REVERSE,
        STRIKE_SLIP
    }

    static class HelperMethods
    {
        // method to interpolate double yInterp from input vectors
        public static double InterpFromVector(double[] xVector, double[] yVector, double xInterp)
        {
            // if xInterp outside of xVector range, return first or last value of yVector
            if (xInterp <= xVector.First())
            {
                return yVector.First();
            }
            else if (xInterp >= xVector.Last())
            {
                return yVector.Last();
            }

            // otherwise, linearly interpolate
            double x1, x2, y1, y2;
            int ind = 0;

            // stop at first index of xVector greater than xInterp
            while ((xVector[ind] < xInterp) && (ind < xVector.Length))
            {
                ind++;
            }
            if (ind == xVector.Length)
            {
                return Double.NaN;
            }
            x1 = xVector[ind - 1];
            y1 = yVector[ind - 1];
            x2 = xVector[ind];
            y2 = yVector[ind];

            return y1 + (xInterp - x1) * (y2 - y1) / (x2 - x1);

        }
    }
    

}