using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Hazard;

namespace GMPEs
{
    public class Campbell_2003_AttenRel
    {
        public readonly  string SHORT_NAME = "Campbell2003";

        // coefficients:
        // @formatter:off
        // some coefficients are labeled differnetly than in paper
        // localCoeff(paperCoeff):
        // c5(c7) c6(c8) c7(c9) c8(c10) c9(c5) c10(c6)
        private double[] pd = { 0.0, 0.2, 1.0, 0.1, 0.3, 0.4, 0.5, 2.0, 0.03, 0.04, 0.05 };
        private double[] c1 = { 0.4492, 0.1325, -0.3177, 0.4064, -0.1483, -0.17039, -0.1333, -1.2483, 1.68, 1.28, 0.87 };
        private double[] c1h = { 0.0305, -0.4328, -0.6104, -0.1475, -0.6906, -0.67076736, -0.5907, -1.4306, 1.186, 0.72857, 0.3736 };
        private double[] c2 = { 0.633, 0.617, 0.451, 0.613, 0.609, 0.5722637, 0.534, 0.459, 0.622, 0.618622, 0.616 };
        private double[] c3 = { -0.0427, -0.0586, -0.2090, -0.0353, -0.0786, -0.10939892, -0.1379, -0.2552, -0.0362, -0.035693, -0.0353 };
        private double[] c4 = { -1.591, -1.32, -1.158, -1.369, -1.28, -1.244142, -1.216, -1.124, -1.691, -1.5660, -1.469 };
        private double[] c5 = { 0.683, 0.399, 0.299, 0.484, 0.349, 0.32603806, 0.318, 0.310, 0.922, 0.75759, 0.630 };
        private double[] c6 = { 0.416, 0.493, 0.503, 0.467, 0.502, 0.5040741, 0.503, 0.499, 0.376, 0.40246, 0.423 };
        private double[] c7 = { 1.140, 1.25, 1.067, 1.096, 1.241, 1.1833254, 1.116, 1.015, 0.759, 0.76576, 0.771 };
        private double[] c8 = { -0.873, -0.928, -0.482, -1.284, -0.753, -0.6529481, -0.606, -0.4170, -0.922, -1.1005, -1.239 };
        private double[] c9 = { -0.00428, -0.0046, -0.00255, -0.00454, -0.00414, -3.7463151E-3, -0.00341, -0.00187, -0.00367, -0.0037319, -0.00378 };
        private double[] c10 = { 0.000483, 0.000337, 0.000141, 0.00046, 0.000263, 2.1878805E-4, 0.000194, 0.000103, 0.000501, 0.00050044, 0.0005 };
        private double[] c11 = { 1.030, 1.077, 1.110, 1.059, 1.081, 1.0901983, 1.098, 1.093, 1.03, 1.037, 1.042 };
        private double[] c12 = { -0.0860, -0.0838, -0.0793, -0.0838, -0.0838, -0.083180725, -0.0824, -0.0758, -0.086, -0.0848, -0.0838 };
        private double[] c13 = { 0.414, 0.478, 0.543, 0.460, 0.482, 0.49511834, 0.508, 0.551, 0.414, 0.43, 0.443 };

        private Dictionary<double, int> indexFromPerHashMap = new Dictionary<double, int> { };

        private int iper;
        private double rRup, mag;
        private SiteType siteType;


        public void setParamDefaults()
        {

            siteType = HazardCalculation.ThisScenario.SiteType;
            mag = HazardCalculation.ThisScenario.Magnitude;
            rRup = HazardCalculation.ThisScenario.RuptureDistance;
        }

        //public enum SiteType
        //{

        //    FIRM_ROCK,
        //    HARD_ROCK

        //}

        public Campbell_2003_AttenRel()
        {

            for (int i = 0; i < pd.Length; i++)
            {
                indexFromPerHashMap.Add(pd[i], i);
            }
            setParamDefaults();

        }

        private void setCoeffIndex()
        {
            iper = indexFromPerHashMap[HazardCalculation.ThisScenario.saPeriodParam];
        }


        public double getMean()
        {
            setCoeffIndex();
            return getMean(iper, siteType, rRup, mag);
        }

        public double getStdDev()
        {
            return getStdDev(iper, mag);
        }

        private double getMean(int iper, SiteType st, double rRup, double mag)
        {
            double period = pd[iper];
            double gnd0 = (st == SiteType.HARD_ROCK) ? c1h[iper] : c1[iper];
            // if (magType == LG_PHASE) mag = Utils.mblgToMw(magConvCode, mag);
            double gndm = gnd0 + c2[iper] * mag + c3[iper] * (8.5 - mag) *
                (8.5 - mag);
            double cfac = Math.Pow((c5[iper] * Math.Exp(c6[iper] * mag)), 2);


            double arg = Math.Sqrt(rRup * rRup + cfac);
            double fac = 0.0;
            if (rRup > 70.0) fac = c7[iper] * (Math.Log(rRup) - Math.Log(70));
            if (rRup > 130.0) fac = fac + c8[iper] * (Math.Log(rRup) - Math.Log(130));
            double gnd = gndm + c4[iper] * Math.Log(arg) + fac +
                (c9[iper] + c10[iper] * mag) * rRup;

            return gnd;
        }

        private double getStdDev(int iper, double mag)
        {
            // if (magType == LG_PHASE) mag = Utils.mblgToMw(magConvCode, mag);
            return (mag < 7.16) ? c11[iper] + c12[iper] * mag : c13[iper];
        }

    }
}
