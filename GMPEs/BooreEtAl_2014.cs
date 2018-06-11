using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Hazard;

// NGA West2 GMPE: Boore, Stewart, Seyhan, Atkinson (BSSA) 2014
// implemented, compiles, but not yet tested
// no range checks yet implemented for parameters (mag, dist, per, etc)

namespace GMPEs
{
    class BSSA2014_AttenRel
    {
        public readonly string SHORT_NAME = "BSSA2014";

        // coefficients and constants:
        private static double A = Math.Pow(570.94, 4);
        private static double B = Math.Pow(1360, 4) + A;
        private static double M_REF = 4.5;
        private static double R_REF = 1.0;
        private static double DC3_CA_TW = 0.0;
        private static double V_REF = 760.0;
        private static double F1 = 0.0;
        private static double F3 = 0.1;
        private static double V1 = 225;
        private static double V2 = 300;

        private static double[] pd = { 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10, 0, -1 };
        private static double[] e0 = { 0.4534, 0.48598, 0.56916, 0.75436, 0.96447, 1.1268, 1.3095, 1.3255, 1.2766, 1.2217, 1.1046, 0.96991, 0.66903, 0.3932, -0.14954, -0.58669, -1.1898, -1.6388, -1.966, -2.5865, -3.0702, 0.4473, 5.037 };
        private static double[] e1 = { 0.4916, 0.52359, 0.6092, 0.79905, 1.0077, 1.1669, 1.3481, 1.359, 1.3017, 1.2401, 1.1214, 0.99106, 0.69737, 0.4218, -0.11866, -0.55003, -1.142, -1.5748, -1.8882, -2.4874, -2.9537, 0.4856, 5.078 };
        private static double[] e2 = { 0.2519, 0.29707, 0.40391, 0.60652, 0.77678, 0.8871, 1.0648, 1.122, 1.0828, 1.0246, 0.89765, 0.7615, 0.47523, 0.207, -0.3138, -0.71466, -1.23, -1.6673, -2.0245, -2.8176, -3.3776, 0.2459, 4.849 };
        private static double[] e3 = { 0.4599, 0.48875, 0.55783, 0.72726, 0.9563, 1.1454, 1.3324, 1.3414, 1.3052, 1.2653, 1.1552, 1.012, 0.69173, 0.4124, -0.1437, -0.60658, -1.2664, -1.7516, -2.0928, -2.6854, -3.1726, 0.4539, 5.033 };
        private static double[] e4 = { 1.421, 1.4331, 1.4261, 1.3974, 1.4174, 1.4293, 1.2844, 1.1349, 1.0166, 0.95676, 0.96766, 1.0384, 1.2871, 1.5004, 1.7622, 1.9152, 2.1323, 2.204, 2.2299, 2.1187, 1.8837, 1.431, 1.073 };
        private static double[] e5 = { 0.04932, 0.053388, 0.061444, 0.067357, 0.073549, 0.055231, -0.042065, -0.11096, -0.16213, -0.1959, -0.22608, -0.23522, -0.21591, -0.18983, -0.1467, -0.11237, -0.04332, -0.014642, -0.014855, -0.081606, -0.15096, 0.05053, -0.1536 };
        private static double[] e6 = { -0.1659, -0.16561, -0.1669, -0.18082, -0.19665, -0.19838, -0.18234, -0.15852, -0.12784, -0.092855, -0.023189, 0.029119, 0.10829, 0.17895, 0.33896, 0.44788, 0.62694, 0.76303, 0.87314, 1.0121, 1.0651, -0.1662, 0.2252 };
        private static double[] Mh = { 5.5, 5.5, 5.5, 5.5, 5.5, 5.54, 5.74, 5.92, 6.05, 6.14, 6.2, 6.2, 6.2, 6.2, 6.2, 6.2, 6.2, 6.2, 6.2, 6.2, 6.2, 5.5, 6.2 };
        private static double[] c1 = { -1.134, -1.1394, -1.1421, -1.1159, -1.0831, -1.0652, -1.0532, -1.0607, -1.0773, -1.0948, -1.1243, -1.1459, -1.1777, -1.193, -1.2063, -1.2159, -1.2179, -1.2162, -1.2189, -1.2543, -1.3253, -1.134, -1.243 };
        private static double[] c2 = { 0.1916, 0.18962, 0.18842, 0.18709, 0.18225, 0.17203, 0.15401, 0.14489, 0.13925, 0.13388, 0.12512, 0.12015, 0.11054, 0.10248, 0.09645, 0.09636, 0.09764, 0.10218, 0.10353, 0.12507, 0.15183, 0.1917, 0.1489 };
        private static double[] c3 = { -0.008088, -0.008074, -0.008336, -0.009819, -0.01058, -0.0102, -0.008977, -0.007717, -0.006517, -0.005475, -0.004053, -0.00322, -0.001931, -0.00121, -0.000365, 0, 0, -0.000052, 0, 0, 0, -0.008088, -0.00344 };
        private static double[] h = { 4.5, 4.5, 4.49, 4.2, 4.04, 4.13, 4.39, 4.61, 4.78, 4.93, 5.16, 5.34, 5.6, 5.74, 6.18, 6.54, 6.93, 7.32, 7.78, 9.48, 9.66, 4.5, 5.3 };
        private static double[] c = { -0.6037, -0.5739, -0.5341, -0.458, -0.4441, -0.4872, -0.5796, -0.6876, -0.7718, -0.8417, -0.9109, -0.9693, -1.0154, -1.05, -1.0454, -1.0392, -1.0112, -0.9694, -0.9195, -0.7766, -0.6558, -0.6, -0.84 };
        private static double[] Vc = { 1500.2, 1500.36, 1502.95, 1501.42, 1494, 1479.12, 1442.85, 1392.61, 1356.21, 1308.47, 1252.66, 1203.91, 1147.59, 1109.95, 1072.39, 1009.49, 922.43, 844.48, 793.13, 771.01, 775, 1500, 1300 };
        private static double[] f4 = { -0.1483, -0.1471, -0.1549, -0.1963, -0.2287, -0.2492, -0.2571, -0.2466, -0.2357, -0.2191, -0.1958, -0.1704, -0.1387, -0.1052, -0.0679, -0.0361, -0.0136, -0.0032, -0.0003, -0.0001, 0, -0.15, -0.1 };
        private static double[] f5 = { -0.00701, -0.00728, -0.00735, -0.00647, -0.00573, -0.0056, -0.00585, -0.00614, -0.00644, -0.0067, -0.00713, -0.00744, -0.00812, -0.00844, -0.00771, -0.00479, -0.00183, -0.00152, -0.00144, -0.00137, -0.00136, -0.00701, -0.00844 };
        private static double[] f6 = { -9.9, -9.9, -9.9, -9.9, -9.9, -9.9, -9.9, -9.9, -9.9, -9.9, -9.9, -9.9, 0.092, 0.367, 0.638, 0.871, 1.135, 1.271, 1.329, 1.329, 1.183, -9.9, -9.9 };
        private static double[] f7 = { -9.9, -9.9, -9.9, -9.9, -9.9, -9.9, -9.9, -9.9, -9.9, -9.9, -9.9, -9.9, 0.059, 0.208, 0.309, 0.382, 0.516, 0.629, 0.738, 0.809, 0.703, -9.9, -9.9 };
        private static double[] r1 = { 111.67, 113.1, 112.13, 97.93, 85.99, 79.59, 81.33, 90.91, 97.04, 103.15, 106.02, 105.54, 108.39, 116.39, 125.38, 130.37, 130.36, 129.49, 130.22, 130.72, 130, 110, 105 };
        private static double[] r2 = { 270, 270, 270, 270, 270.04, 270.09, 270.16, 270, 269.45, 268.59, 266.54, 265, 266.51, 270, 262.41, 240.14, 195, 199.45, 230, 250.39, 210, 270, 272 };
        private static double[] dPhiR = { 0.096, 0.092, 0.081, 0.063, 0.064, 0.087, 0.12, 0.136, 0.141, 0.138, 0.122, 0.109, 0.1, 0.098, 0.104, 0.105, 0.088, 0.07, 0.061, 0.058, 0.06, 0.1, 0.082 };
        private static double[] dPhiV = { 0.07, 0.03, 0.029, 0.03, 0.022, 0.014, 0.015, 0.045, 0.055, 0.05, 0.049, 0.06, 0.07, 0.02, 0.01, 0.008, 0, 0, 0, 0, 0, 0.07, 0.08 };
        private static double[] phi1 = { 0.698, 0.702, 0.721, 0.753, 0.745, 0.728, 0.72, 0.711, 0.698, 0.675, 0.643, 0.615, 0.581, 0.553, 0.532, 0.526, 0.534, 0.536, 0.528, 0.512, 0.51, 0.695, 0.644 };
        private static double[] phi2 = { 0.499, 0.502, 0.514, 0.532, 0.542, 0.541, 0.537, 0.539, 0.547, 0.561, 0.58, 0.599, 0.622, 0.625, 0.619, 0.618, 0.619, 0.616, 0.622, 0.634, 0.604, 0.495, 0.552 };
        private static double[] tau1 = { 0.402, 0.409, 0.445, 0.503, 0.474, 0.415, 0.354, 0.344, 0.35, 0.363, 0.381, 0.41, 0.457, 0.498, 0.525, 0.532, 0.537, 0.543, 0.532, 0.511, 0.487, 0.398, 0.401 };
        private static double[] tau2 = { 0.345, 0.346, 0.364, 0.426, 0.466, 0.458, 0.388, 0.309, 0.266, 0.229, 0.21, 0.224, 0.266, 0.298, 0.315, 0.329, 0.344, 0.349, 0.335, 0.27, 0.239, 0.348, 0.346 };

        // initialize dictionary of period indices. double is period, int is index.
        private Dictionary<double, int> indexFromPerHashMap = new Dictionary<double, int> { };

        private int iper;
        private double period, rJB, mag, z1p0, vs30;
        private FaultStyle style;

        public void setParamDefaults()
        {
            period = HazardCalculation.ThisScenario.saPeriodParam;
            mag = HazardCalculation.ThisScenario.Magnitude;
            rJB = HazardCalculation.ThisScenario.JoynerBooreDistance;
            vs30 = HazardCalculation.ThisScenario.VsThirty;
            z1p0 = HazardCalculation.ThisScenario.Z1p0;
            style = HazardCalculation.ThisScenario.FaultStyle;
        }

        public BSSA2014_AttenRel()
        {
            // populate dictionary of period indices
            for (int i = 0; i < pd.Length; i++)
            {
                indexFromPerHashMap.Add(pd[i], i);
            }
            setParamDefaults();

        }

        public double getMean()
        {
            // Check if key (period) is directly available from GMPE
            period = HazardCalculation.ThisScenario.saPeriodParam;
            if (indexFromPerHashMap.ContainsKey(period))
            {
                // Get median directly from GMPE
                //setCoeffIndex(period);
                return getMeanHere();
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
                    return 0;
                }
                else if (period > perList.Last())
                {
                    // desired period is not in range of GMPE
                    return 0;
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
                    double[] meanVect = { 0, 0 };

                    period = perList[ind - 1];
                    meanVect[0] = getMeanHere();

                    period = perList[ind];
                    meanVect[1] = getMeanHere();

                    period = HazardCalculation.ThisScenario.saPeriodParam;

                    //interpolate in log-log space
                    return HelperMethods.InterpFromVector(logPerVect, meanVect, Math.Log(period));
                }

            }

        }

        public double getStdDev()
        {
            // Check if key (period) is directly available from GMPE
            period = HazardCalculation.ThisScenario.saPeriodParam;
            if (indexFromPerHashMap.ContainsKey(period))
            {
                // Get median directly from GMPE
                //setCoeffIndex(period);
                return getStdDevHere();
            }
            // If key (period) is not directly available from GMPE, interpolate
            // using median of next highest and lowest periods that are directly available
            else
            {
                // created ordered list of keys
                var perList = indexFromPerHashMap.Keys.ToList();
                perList.Sort();

                // find two indicies for keys (periods) just above and below saPeriodParam
                if (period < perList.First())
                {
                    // desired period is not in range of GMPE
                    return 0;
                }
                else if (period > perList.Last())
                {
                    // desired period is not in range of GMPE
                    return 0;
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
                    double[] sigVect = { 0, 0 };

                    period = perList[ind - 1];
                    sigVect[0] = getStdDevHere();

                    period = perList[ind];
                    sigVect[1] = getStdDevHere();

                    period = HazardCalculation.ThisScenario.saPeriodParam;

                    //interpolate in log-log space
                    return HelperMethods.InterpFromVector(logPerVect, sigVect, Math.Log(period));
                }

            }

        }

        // Get index 
        private void setCoeffIndex(double perKey)
        {
            iper = indexFromPerHashMap[perKey];
        }

        private double getMeanHere()
        {
            // ****** Mean ground motion model ******

            // use PGA coefficients to get pgaRock
            setCoeffIndex(0.0);
            double pgaRock = calcPGArock();

            // use current period coefficients for remaining calculations
            setCoeffIndex(period);

            // Source/Event Term -- Equation 2
            double Fe = calcSourceTerm();

            // Path Term -- Equations 3, 4
            double R = Math.Sqrt(rJB * rJB + h[iper] * h[iper]);
            double Fp = calcPathTerm(R);

            // Site Linear Term -- Equation 6
            double vsLin = (vs30 <= Vc[iper]) ? vs30 : Vc[iper];
            double lnFlin = c[iper] * Math.Log(vsLin / V_REF);

            // Site Nonlinear Term -- Equations 7, 8
            double f2 = f4[iper] * (Math.Exp(f5[iper] * (Math.Min(vs30, 760.0) - 360.0)) - Math.Exp(f5[iper] * (760.0 - 360.0)));
            double lnFnl = F1 + f2 * Math.Log((pgaRock + F3) / F3);

            // Basin depth term -- Equations 9, 10 , 11
            double DZ1 = calcDeltaZ1();
            double Fdz1 = (period >= 0.65) ? (DZ1 <= f7[iper] / f6[iper]) ? f6[iper] * DZ1 : f7[iper] : 0.0;

            // Total site term -- Equation 5
            double Fs = lnFlin + lnFnl + Fdz1;

            // Total model -- Equation 1
            return Fe + Fp + Fs;
        }

        // Median PGA for ref rock (Vs30=760m/s); always called with PGA coeffs
        private double calcPGArock()
        {
            // Source/Event Term -- Equation 2
            double FePGA = calcSourceTerm();

            // Path Term -- Equation 3
            double R = Math.Sqrt(rJB * rJB + h[iper] * h[iper]);
            double FpPGA = calcPathTerm(R);

            // No Site term -- [Vs30rk==760] < [Vc(PGA)=1500] &&
            // ln(Vs30rk / V_REF) = ln(760/760) = 0

            // Total PGA model -- Equation 1
            return Math.Exp(FePGA + FpPGA);
        }

        // Aleatory uncertainty model
        private double getStdDevHere()
        {
            // get index coefficients for current period
            setCoeffIndex(period);

            // Inter-event Term -- Equation 14
            double tau = (mag >= 5.5) ? tau2[iper] : (mag <= 4.5) ? tau1[iper] : tau1[iper] + (tau2[iper] - tau1[iper]) * (mag - 4.5);

            // Intra-event Term -- Equations 15, 16, 17
            double phi_m = (mag >= 5.5) ? phi2[iper] : (mag <= 4.5) ? phi1[iper] : phi1[iper] + (phi2[iper] - phi1[iper]) * (mag - 4.5);

            double phi_mr = phi_m;
            if (rJB > r2[iper])
            {
                phi_mr += dPhiR[iper];
            }
            else if (rJB > r1[iper])
            {
                phi_mr += dPhiR[iper] * (Math.Log(rJB / r1[iper]) / Math.Log(r2[iper] / r1[iper]));
            }

            double phi_mrv = phi_mr;
            if (vs30 <= V1)
            {
                phi_mrv -= dPhiV[iper];
            }
            else if (vs30 < V2)
            {
                phi_mrv -= dPhiV[iper] * (Math.Log(V2 / vs30) / Math.Log(V2 / V1));
            }

            // Total model -- Equation 13
            return Math.Sqrt(phi_mrv * phi_mrv + tau * tau);
        }

        // Source/Event Term -- Equation 2
        private double calcSourceTerm()
        {
            double Fe = (style == FaultStyle.STRIKE_SLIP) ? e1[iper]
                : (style == FaultStyle.REVERSE) ? e3[iper] : (style == FaultStyle.NORMAL) ? e2[iper] : e0[iper]; // UNKNOWN
            double MwMh = mag - Mh[iper];
            Fe += (mag <= Mh[iper]) ? e4[iper] * MwMh + e5[iper] * MwMh * MwMh : e6[iper] * MwMh;
            return Fe;
        }

        // Path Term, base model -- Equation 3
        private double calcPathTerm(double R)
        {
            return (c1[iper] + c2[iper] * (mag - M_REF)) * Math.Log(R / R_REF) +
                (c3[iper] + DC3_CA_TW) * (R - R_REF);
        }

        // Calculate delta Z1 in km as a function of vs30 and using the default
        // model of ChiouYoungs_2013 -- Equations 10, 11
        private double calcDeltaZ1()
        {
            if (Double.IsNaN(z1p0))
            {
                return 0.0;
            }
            double vsPow4 = vs30 * vs30 * vs30 * vs30;
            return z1p0 - Math.Exp(-7.15 / 4.0 * Math.Log((vsPow4 + A) / B)) / 1000.0;
        }

    }

}
