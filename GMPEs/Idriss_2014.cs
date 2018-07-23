using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Hazard;

// NGA West2 GMPE: Idriss (I) 2014
// implemented, compiles, but not yet tested
// no range checks yet implemented for parameters (mag, dist, per, etc)

namespace GMPEs
{
    class I2014_AttenRel
    {
        public readonly string SHORT_NAME = "I2014";

        // coefficients:
        private static double[] pd = { 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10, 0 };
        private static double[] a1_lo = { 7.0887,7.1157,7.2087,6.2638,5.9051,7.5791,8.019,9.2812,9.5804,9.8912,9.5342,9.2142,8.3517,7.0453,5.1307,3.361,0.1784,-2.4301,-4.357,-7.8275,-9.2857, 7.0887 };
        private static double[] a2_lo = { 0.2058,0.2058,0.2058,0.0625,0.1128,0.0848,0.1713,0.1041,0.0875,0.0003,0.0027,0.0399,0.0689,0.16,0.2429,0.3966,0.756,0.9283,1.1209,1.4016,1.5574, 0.2058 };
        private static double[] b1_lo = { 2.9935,2.9935,2.9935,2.8664,2.9406,3.019,2.7871,2.8611,2.8289,2.8423,2.83,2.856,2.7544,2.7339,2.68,2.6837,2.6907,2.5782,2.5468,2.4478,2.3922, 2.9935 };
        private static double[] b2_lo = { -0.2287,-0.2287,-0.2287,-0.2418,-0.2513,-0.2516,-0.2236,-0.2229,-0.22,-0.2284,-0.2318,-0.2337,-0.2392,-0.2398,-0.2417,-0.245,-0.2389,-0.2514,-0.2541,-0.2593,-0.2586,-0.2287 };
        private static double[] a1_hi = { 9.0138,9.0408,9.1338,7.9837,7.756,9.4252,9.6242,11.13,11.3629,11.7818,11.6097,11.4484,10.9065,9.8565,8.3363,6.8656,4.1178,1.8102,0.0977,-3.0563,-4.4387, 9.0138 };
        private static double[] a2_hi = { -0.0794,-0.0794,-0.0794,-0.1923,-0.1614,-0.1887,-0.0665,-0.1698,-0.1766,-0.2798,-0.3048,-0.2911,-0.3097,-0.2565,-0.232,-0.1226,0.1724,0.3001,0.4609,0.6948,0.8393,-0.0794 };
        private static double[] b1_hi = { 2.9935,2.9935,2.9935,2.7995,2.8143,2.8131,2.4091,2.4938,2.3773,2.3772,2.3413,2.3477,2.2042,2.1493,2.0408,2.0013,1.9408,1.7763,1.703,1.5212,1.4195, 2.9935 };
        private static double[] b2_hi = { -0.2287,-0.2287,-0.2287,-0.2319,-0.2326,-0.2211,-0.1676,-0.1685,-0.1531,-0.1595,-0.1594,-0.1584,-0.1577,-0.1532,-0.147,-0.1439,-0.1278,-0.1326,-0.1291,-0.122,-0.1145,-0.2287 };
        private static double[] a3 = { 0.0589,0.0589,0.0589,0.0417,0.0527,0.0442,0.0329,0.0188,0.0095,-0.0039,-0.0133,-0.0224,-0.0267,-0.0198,-0.0367,-0.0291,-0.0214,-0.024,-0.0202,-0.0219,-0.0035, 0.0589 };
        private static double[] xi = { -0.854,-0.854,-0.854,-0.631,-0.591,-0.757,-0.911,-0.998,-1.042,-1.03,-1.019,-1.023,-1.056,-1.009,-0.898,-0.851,-0.761,-0.675,-0.629,-0.531,-0.586,-0.854 };
        private static double[] gamma = { -0.0027,-0.0027,-0.0027,-0.0061,-0.0056,-0.0042,-0.0046,-0.003,-0.0028,-0.0029,-0.0028,-0.0021,-0.0029,-0.0032,-0.0033,-0.0032,-0.0031,-0.0051,-0.0059,-0.0057,-0.0061,-0.0027 };
        private static double[] phi = { 0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.06,0.04,0.02,0.02,0,0,0,0, 0.08 };

        // initialize dictionary of period indices. double is period, int is index.
        private Dictionary<double, int> indexFromPerHashMap = new Dictionary<double, int> { };

        private int iper;
        private double period, rRup, mag, vs30;
        private FaultStyle style;

        public void setParamDefaults()
        {
            mag = HazardCalculation.ThisScenario.Magnitude;
            rRup = HazardCalculation.ThisScenario.RuptureDistance;
            vs30 = HazardCalculation.ThisScenario.VsThirty;
            style = HazardCalculation.ThisScenario.FaultStyle;
        }

        public I2014_AttenRel()
        {
            // populate dictionary of period indices
            for (int i = 0; i < pd.Length; i++)
            {
                indexFromPerHashMap.Add(pd[i], i);
            }
            setParamDefaults();

        }

        public GroundMotion GetGroundMotion(GroundMotion newGm)
        {
            double mu, sig;

            // Check if key (period) is directly available from GMPE
            period = HazardCalculation.ThisScenario.saPeriodParam;
            if (indexFromPerHashMap.ContainsKey(period))
            {
                // Get median directly from GMPE
                mu = getMeanHere();
                sig = getStdDevHere();
            }
            // If key (period) is not directly available from GMPE, interpolate
            // using median of next highest and lowest periods that are directly available
            else
            {
                // created ordered list of keys (ascending order)
                var perList = indexFromPerHashMap.Keys.ToList();
                perList.Sort(); // ascending is default

                // find two indicies for keys (periods) just above and below period
                if (period < perList.First())
                {
                    // desired period is not in range of GMPE
                    return newGm;
                }
                else if (period > perList.Last())
                {
                    // desired period is not in range of GMPE
                    return newGm;
                }
                else
                {
                    int ind = 0;
                    while ((perList[ind] < period) && (ind < perList.Count()))
                    {
                        ind++;
                    }

                    // create log period and log mean vectors for interpolation
                    double[] logPerVect = { Math.Log(perList[ind - 1]), Math.Log(perList[ind]) };
                    double[] muVect = { 0, 0 };
                    double[] sigVect = { 0, 0 };

                    period = perList[ind - 1];
                    muVect[0] = getMeanHere();
                    sigVect[0] = getStdDevHere();

                    period = perList[ind];
                    muVect[1] = getMeanHere();
                    sigVect[1] = getStdDevHere();

                    period = HazardCalculation.ThisScenario.saPeriodParam;

                    //interpolate in log-log space
                    mu = HelperMethods.InterpFromVector(logPerVect, muVect, Math.Log(period));
                    sig = HelperMethods.InterpFromVector(logPerVect, sigVect, Math.Log(period));
                }

            }
            newGm.SetLogMean(mu);
            newGm.SetLogStd(sig);

            return newGm;
        }

        private double getMeanHere()
        {
            // get index coefficients for current period
            setCoeffIndex(period);

            // ****** Mean ground motion model ******
            double a1 = a1_lo[iper], a2 = a2_lo[iper];
            double b1 = b1_lo[iper], b2 = b2_lo[iper];
            if (mag > 6.75)
            {
                a1 = a1_hi[iper];
                a2 = a2_hi[iper];
                b1 = b1_hi[iper];
                b2 = b2_hi[iper];
            }

            return a1 + a2 * mag + a3[iper] * (8.5 - mag) * (8.5 - mag) - (b1 + b2 * mag) * Math.Log(rRup + 10.0) +
                xi[iper] * Math.Log(Math.Min(vs30, 1200.0)) + gamma[iper] * rRup + (style == FaultStyle.REVERSE ? phi[iper] : 0.0);
        }

        // Aleatory uncertainty model
        private double getStdDevHere()
        {
            // get index coefficients for current period
            setCoeffIndex(period);

            double s1 = 0.035;
            s1 *= (period <= 0.05) ? Math.Log(0.05) :      // omitted if T == null
                (period < 3.0) ? Math.Log(period) : Math.Log(3d);
            double s2 = 0.06;
            s2 *= (mag <= 5.0) ? 5.0 : (mag < 7.5) ? mag : 7.5;
            return 1.18 + s1 - s2;
        }

        // Get index for period in coefficient arrays 
        private void setCoeffIndex(double perKey)
        {
            iper = indexFromPerHashMap[perKey];
        }

    }

}
