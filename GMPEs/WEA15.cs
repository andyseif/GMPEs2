using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Hazard;

namespace GMPEs
{
    public class WongEtAl_2015_AttenRel
    {
        public readonly string SHORT_NAME = "WongEtAl_2015_AttenRel";
        //Periods
        private double[] pd = { 0.0, 0.01, 0.02, 0.025, 0.03, 0.04, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 2.0, 3.0, 5.0, 10.0 };
        //Coefficients
        private double[] c1 = { 68.52187, 68.77858, 65.29077, 67.40193, 69.46166, 78.01486, 85.62312, 103.0432, 108.8136, 98.22293, 73.99161, 58.53881, 46.91267, 43.47321, 32.93593, 21.84053, 9.41267, -0.08881, -6.9155 };
        private double[] c2 = { -5.09631, -5.13883, -4.98649, -5.24299, -5.48027, -6.22686, -6.76578, -7.20814, -6.61834, -5.32783, -2.90241, -1.71001, -0.93125, -0.50426, 0.08357, 0.6191, 1.23611, 1.63345, 1.66675 };
        private double[] c3 = { 5.8, 5.8, 5.7, 5.7, 5.7, 5.8, 5.9, 6.2, 6.4, 6.4, 6.3, 6.2, 6.1, 6.2, 6.1, 6.1, 5.9, 5.7, 5.4 };
        private double[] c4 = { -12.9601, -13.00533, -12.56134, -12.88585, -13.1943, -14.39446, -15.41352, -17.48802, -17.96824, -16.41981, -13.06945, -10.95674, -9.36813, -8.89132, -7.51731, -6.22484, -4.64883, -3.56561, -3.03522 };
        private double[] c5 = { 1.03629, 1.04366, 1.02752, 1.06672, 1.10188, 1.20623, 1.27568, 1.29736, 1.18618, 1.00088, 0.66484, 0.50386, 0.40094, 0.35245, 0.28394, 0.24654, 0.18193, 0.15429, 0.18615 };
        private double[] c6 = { -0.14898, -0.14918, -0.14202, -0.13721, -0.13149, -0.12877, -0.13241, -0.16376, -0.1876, -0.20822, -0.2386, -0.25817, 0.27491, -0.30295, -0.31387, -0.30737, -0.28698, -0.23807, -0.15158 };
        private double[] csigma = { 0.7803, 0.7816, 0.7975, 0.8063, 0.8194, 0.828, 0.8327, 0.8254, 0.829, 0.8381, 0.8492, 0.8412, 0.8188, 0.8183, 0.7954, 0.9512, 1.0442, 1.1854, 1.3092 };

        private Dictionary<double, int> indexFromPerHashMap = new Dictionary<double, int> { };

        private int iper;
        private double rjb, mag, rMin = 20;

        public void setParamDefaults()
        {

            mag = HazardCalculation.ThisScenario.Magnitude;
            rjb = HazardCalculation.ThisScenario.JoynerBooreDistance;
            // zTor = HazardCalculation.ThisScenario.RuptureTop;
            //Need to incorporate constraints; this GMM is only valid for ZTOR GTE 20 km and LTE 60 km:
            //.set(MW, Range.closed(5.0, 8.0))
            //.set(RRUP, Range.closed(0.0, 300.0))
            //.set(ZTOP, Range.closed(20.0, 60.0))
            //.set(VS30, Range.singleton(760.0))

        }

        public WongEtAl_2015_AttenRel()
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
            return getMeanLocal();
        }

        public double getStdDev()
        {
            return getStdDev(iper);
        }

        private double getMeanLocal()
        {
            double period = pd[iper];
            return c1[iper] + c2[iper] * mag + (c4[iper] + c5[iper] * mag) * Math.Log(Math.Max(rMin, rjb) + Math.Exp(c3[iper])) + c6[iper] * (mag - 6.0) * (mag - 6.0);

        }

        private double getStdDev(int iper)
        {
            return csigma[iper];
        }
    }
}
