using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Hazard;

// NGA West2 GMPE: Chiou Youngs (CY) 2014
// not fully implemented
// no range checks yet implemented for parameters (mag, dist, per, etc)

namespace GMPEs
{
    class CY2014_AttenRel
    {
        public readonly string SHORT_NAME = "CY2014";

        // coefficients and constants:
        private static double C2 = 1.06;
        private static double C4 = -2.1;
        private static double C4A = -0.5;
        private static double dC4 = C4A - C4;
        private static double CRB = 50.0;
        private static double CRBsq = CRB * CRB;
        private static double C11 = 0.0;
        private static double PHI6 = 300.0;
        private static double A = Math.Pow(571, 4);
        private static double B = Math.Pow(1360, 4) + A;

        private static double[] pd = { 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10, 0, -1 };

        private static double[] c1 = { -1.5065,-1.4798,-1.2972,-1.1007,-0.9292,-0.658,-0.5613,-0.5462,-0.6798,-0.8663,-1.0514,-1.3794,-1.6508,-2.1511,-2.5365,-3.0686,-3.4148,-3.9013,-4.2466,-4.5143,-5.0009,-5.3461,-1.5065,2.3549 };
        private static double[] c1a = { 0.165,0.165,0.165,0.165,0.165,0.165,0.165,0.165,0.165,0.165,0.165,0.165,0.165,0.165,0.165,0.165,0.1645,0.1168,0.0732,0.0484,0.022,0.0124,0.165,0.165 };
        private static double[] c1b = { -0.255,-0.255,-0.255,-0.255,-0.255,-0.254,-0.253,-0.25,-0.2449,-0.2382,-0.2313,-0.2146,-0.1972,-0.162,-0.14,-0.1184,-0.11,-0.104,-0.102,-0.101,-0.101,-0.1,-0.255,-0.0626 };
        private static double[] c1c = { -0.165,-0.165,-0.165,-0.165,-0.165,-0.165,-0.165,-0.165,-0.165,-0.165,-0.165,-0.165,-0.165,-0.165,-0.165,-0.165,-0.1645,-0.1168,-0.0732,-0.0484,-0.022,-0.0124,-0.165,-0.165 };
        private static double[] c1d = { 0.255,0.255,0.255,0.255,0.255,0.254,0.253,0.25,0.2449,0.2382,0.2313,0.2146,0.1972,0.162,0.14,0.1184,0.11,0.104,0.102,0.101,0.101,0.1,0.255,0.0626 };
        private static double[] cn = { 16.0875,15.7118,15.8819,16.4556,17.6453,20.1772,19.9992,16.6246,13.7012,11.2667,9.1908,6.5459,5.2305,3.7896,3.3024,2.8498,2.5417,2.1488,1.8957,1.7228,1.5737,1.5265,16.0875,3.3024 };
        private static double[] cM = { 4.9993,4.9993,4.9993,4.9993,4.9993,5.0031,5.0172,5.0547,5.0939,5.1315,5.167,5.2317,5.2893,5.4109,5.5106,5.6705,5.7981,5.9983,6.1552,6.2856,6.5428,6.7415,4.9993,5.423 };
        private static double[] c3 = { 1.9636,1.9636,1.9636,1.9636,1.9636,1.9636,1.9636,2.0362,2.1521,2.2574,2.344,2.4709,2.5567,2.6812,2.7474,2.8161,2.8514,2.8875,2.9058,2.9169,2.932,2.9396,1.9636,2.3152 };
        private static double[] c5 = { 6.4551,6.4551,6.4551,6.4551,6.4551,6.4551,6.8305,7.3621,7.4972,7.5416,7.56,7.5735,7.5778,7.5808,7.5814,7.5817,7.5818,7.5818,7.5818,7.5818,7.5818,7.5818,6.4551,5.8096 };
        private static double[] cHM = { 3.0956,3.0963,3.0974,3.0988,3.1011,3.1094,3.2381,3.43,3.5146,3.5746,3.6232,3.6945,3.7401,3.7941,3.8144,3.8284,3.833,3.8361,3.8369,3.8376,3.838,3.838,3.0956,3.0514 };
        private static double[] c6 = { 0.4908,0.4925,0.4992,0.5037,0.5048,0.5048,0.5048,0.5045,0.5016,0.4971,0.4919,0.4807,0.4707,0.4575,0.4522,0.4501,0.45,0.45,0.45,0.45,0.45,0.45,0.4908,0.4407 };
        private static double[] c7 = { 0.0352,0.0352,0.0352,0.0352,0.0352,0.0352,0.0352,0.0352,0.0352,0.0352,0.0352,0.0352,0.0352,0.0352,0.0352,0.0352,0.0352,0.016,0.0062,0.0029,0.0007,0.0003,0.0352,0.0324 };
        private static double[] c7b = { 0.0462,0.0472,0.0533,0.0596,0.0639,0.063,0.0532,0.0345,0.0202,0.009,-0.0004,-0.0155,-0.0278,-0.0477,-0.0559,-0.063,-0.0665,-0.0516,-0.0448,-0.0424,-0.0348,-0.0253,0.0462,0.0097 };
        private static double[] c9 = { 0.9228,0.9296,0.9396,0.9661,0.9794,1.026,1.0177,0.9801,0.9459,0.9196,0.8829,0.8302,0.7884,0.6754,0.6196,0.5101,0.3917,0.1244,0.0086,0,0,0,0.9228,0.3079 };
        private static double[] c9a = { 0.1202,0.1217,0.1194,0.1166,0.1176,0.1171,0.1146,0.1106,0.1208,0.1208,0.1175,0.106,0.1061,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1202,0.1 };
        private static double[] c9b = { 6.8607,6.8697,6.9113,7.0271,7.0959,7.3298,7.2588,7.2109,7.2988,7.3691,6.8789,6.5334,6.526,6.5,6.5,6.5,6.5,6.5,6.5,6.5,6.5,6.5,6.8607,6.5 };
        private static double[] c11b = { -0.4536,-0.4536,-0.4536,-0.4536,-0.4536,-0.4536,-0.4536,-0.4536,-0.444,-0.3539,-0.2688,-0.1793,-0.1428,-0.1138,-0.1062,-0.102,-0.1009,-0.1003,-0.1001,-0.1001,-0.1,-0.1,-0.4536,-0.3834 };
        private static double[] gamma1 = { -0.007146,-0.007249,-0.007869,-0.008316,-0.008743,-0.009537,-0.00983,-0.009896,-0.009505,-0.008918,-0.008251,-0.007267,-0.006492,-0.005147,-0.004277,-0.002979,-0.002301,-0.001344,-0.001084,-0.00101,-0.000964,-0.00095,-0.007146,-0.001852 };
        private static double[] gamma2 = { -0.006758,-0.006758,-0.006758,-0.006758,-0.006758,-0.00619,-0.005332,-0.003806,-0.00269,-0.002128,-0.001812,-0.001274,-0.001074,-0.001115,-0.001197,-0.001675,-0.002349,-0.003306,-0.003566,-0.00364,-0.003686,-0.0037,-0.006758,-0.007403 };
        private static double[] gamma3 = { 4.2542,4.2386,4.2519,4.296,4.3578,4.5455,4.7603,5.0644,5.188,5.2164,5.1954,5.0899,4.7854,4.3304,4.1667,4.0029,3.8949,3.7928,3.7443,3.709,3.6632,3.623,4.2542,4.3439 };
        private static double[] phi1 = { -0.521,-0.5055,-0.4368,-0.3752,-0.3469,-0.3747,-0.444,-0.5477,-0.6693,-0.7766,-0.8501,-0.9431,-1.0044,-1.0602,-1.0941,-1.1142,-1.1154,-1.1081,-1.0603,-0.9872,-0.8274,-0.7053,-0.521,-0.7936 };
        private static double[] phi2 = { -0.1417,-0.1364,-0.1403,-0.1591,-0.1862,-0.2538,-0.2943,-0.3113,-0.2927,-0.2662,-0.2405,-0.1975,-0.1633,-0.1028,-0.0699,-0.0425,-0.0302,-0.0129,-0.0016,0,0,0,-0.1417,-0.0699 };
        private static double[] phi3 = { -0.00701,-0.007279,-0.007354,-0.006977,-0.006467,-0.005734,-0.005604,-0.005845,-0.006141,-0.006439,-0.006704,-0.007125,-0.007435,-0.00812,-0.008444,-0.007707,-0.004792,-0.001828,-0.001523,-0.00144,-0.001369,-0.001361,-0.00701,-0.008444 };
        private static double[] phi4 = { 0.102151,0.10836,0.119888,0.133641,0.148927,0.190596,0.230662,0.266468,0.255253,0.231541,0.207277,0.165464,0.133828,0.085153,0.058595,0.031787,0.019716,0.009643,0.005379,0.003223,0.001134,0.000515,0.102151,5.41 };
        private static double[] phi5 = { 0,0,0,0,0,0,0,0,0,0,0.001,0.004,0.01,0.034,0.067,0.143,0.203,0.277,0.309,0.321,0.329,0.33,0,0.0202 };
        private static double[] tau1 = { 0.4,0.4026,0.4063,0.4095,0.4124,0.4179,0.4219,0.4275,0.4313,0.4341,0.4363,0.4396,0.4419,0.4459,0.4484,0.4515,0.4534,0.4558,0.4574,0.4584,0.4601,0.4612,0.4,0.3894 };
        private static double[] tau2 = { 0.26,0.2637,0.2689,0.2736,0.2777,0.2855,0.2913,0.2993,0.3047,0.3087,0.3119,0.3165,0.3199,0.3255,0.3291,0.3335,0.3363,0.3398,0.3419,0.3435,0.3459,0.3474,0.26,0.2578 };
        private static double[] sig1 = { 0.4912,0.4904,0.4988,0.5049,0.5096,0.5179,0.5236,0.5308,0.5351,0.5377,0.5395,0.5422,0.5433,0.5294,0.5105,0.4783,0.4681,0.4617,0.4571,0.4535,0.4471,0.4426,0.4912,0.4785 };
        private static double[] sig2 = { 0.3762,0.3762,0.3849,0.391,0.3957,0.4043,0.4104,0.4191,0.4252,0.4299,0.4338,0.4399,0.4446,0.4533,0.4594,0.468,0.4681,0.4617,0.4571,0.4535,0.4471,0.4426,0.3762,0.3629 };
        private static double[] sig3 = { 0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.7999,0.7997,0.7988,0.7966,0.7792,0.7504,0.7136,0.7035,0.7006,0.7001,0.7,0.7,0.7,0.8,0.7504 };

        // initialize dictionary of period indices. double is period, int is index.
        private Dictionary<double, int> indexFromPerHashMap = new Dictionary<double, int> { };

        private int iper;
        private double period, rRup, rJB, rX, mag, zTop, z1p0, dip, vs30;
        private FaultStyle style;
        private bool vsInf;
        
        public void setParamDefaults()
        {
            period = HazardCalculation.ThisScenario.saPeriodParam;
            mag = HazardCalculation.ThisScenario.Magnitude;
            rX = HazardCalculation.ThisScenario.RxDistance;
            rRup = HazardCalculation.ThisScenario.RuptureDistance;
            rJB = HazardCalculation.ThisScenario.JoynerBooreDistance;
            vs30 = HazardCalculation.ThisScenario.VsThirty;
            vsInf = HazardCalculation.ThisScenario.IsInferred;
            zTop = HazardCalculation.ThisScenario.Ztor;
            z1p0 = HazardCalculation.ThisScenario.Z1p0;
            style = HazardCalculation.ThisScenario.FaultStyle;
            dip = HazardCalculation.ThisScenario.Dip;
        }

        public CY2014_AttenRel()
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

            // Check if key (period) is directly available from GMPE
            period = HazardCalculation.ThisScenario.saPeriodParam;
            if (indexFromPerHashMap.ContainsKey(period))
            {
                // Get median directly from GMPE
                newGm = getGmHere(newGm);
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
                    GroundMotion gm = new GroundMotion();

                    period = perList[ind - 1];
                    gm = getGmHere(gm);
                    muVect[0] = gm.GetLogMean();
                    sigVect[0] = gm.GetLogStd();

                    // zero-out gm
                    newGm.SetLogMean(0.0);
                    newGm.SetLogStd(0.0);

                    period = perList[ind];
                    gm = getGmHere(gm);
                    muVect[1] = gm.GetLogMean();
                    sigVect[1] = gm.GetLogStd();

                    period = HazardCalculation.ThisScenario.saPeriodParam;

                    //interpolate in log-log space
                    double mu = HelperMethods.InterpFromVector(logPerVect, muVect, Math.Log(period));
                    double sig = HelperMethods.InterpFromVector(logPerVect, sigVect, Math.Log(period));

                    newGm.SetLogMean(mu);
                    newGm.SetLogStd(sig);
                }

            }

            return newGm;
        }

        // Get index 
        private void setCoeffIndex(double perKey)
        {
            iper = indexFromPerHashMap[perKey];
        }

        private GroundMotion getGmHere(GroundMotion gm)
        {
            setCoeffIndex(period);

            // terms used by both mean and stdDev
            double saRef = CalcSAref();
            double soilNonLin = CalcSoilNonLin();

            double mu = GetMeanForThisT(soilNonLin, saRef);
            double std = GetStdForThisT(soilNonLin, saRef);

            gm.SetLogMean(mu);
            gm.SetLogStd(std);

            return gm;
        }
        
        // Seismic Source Scaling -- Equation 11
        private double CalcSAref()
        {
            // Magnitude scaling
            double r1 = c1[iper] + C2 * (mag - 6.0) + ((C2 - c3[iper]) / cn[iper]) *
                Math.Log(1.0 + Math.Exp(cn[iper] * (cM[iper] - mag)));

            // Near-field magnitude and distance scaling
            double r2 = C4 * Math.Log(rRup + c5[iper] * Math.Cosh(c6[iper] * Math.Max(mag - cHM[iper], 0.0)));

            // Far-field distance scaling
            double gamma = (gamma1[iper] + gamma2[iper] / Math.Cosh(Math.Max(mag - gamma3[iper], 0.0)));
            double r3 = dC4 * Math.Log(Math.Sqrt(rRup * rRup + CRBsq)) + rRup * gamma;

            // Scaling with other source variables
            double coshM = Math.Cosh(2 * Math.Max(mag - 4.5, 0.0));
            double cosdelta = Math.Cos(dip * Math.PI / 180.0);

            // Center zTop on the zTop-M relation
            double deltaZtop = zTop - calcMwZtop();
            double r4 = (c7[iper] + c7b[iper] / coshM) * deltaZtop + (C11 + c11b[iper] / coshM) * cosdelta * cosdelta;
            r4 += (style == FaultStyle.REVERSE) ? (c1a[iper] + c1c[iper] / coshM)
                : (style == FaultStyle.NORMAL) ? (c1b[iper] + c1d[iper] / coshM) : 0.0;

            // Hanging-wall effect
            double r5 = 0.0;
            if (rX >= 0.0) {
                r5 = c9[iper] * Math.Cos(dip * Math.PI / 180.0) *
                    (c9a[iper] + (1.0 - c9a[iper]) * Math.Tanh(rX / c9b[iper])) *
                    (1 - Math.Sqrt(rJB * rJB + zTop * zTop) / (rRup + 1.0));
            }

            // Directivity effect (not implemented)
            // cDPP = centered DPP (direct point directivity parameter)
            // double c8 = 0.2154; // corrected from 2.154 12/3/13 per email from Sanaz
            // double c8a = 0.2695;
            // double Mc8 = mag-c.c8b;
            // double r6 = c8 * exp(-c8a * Mc8 * Mc8) *
            // max(0.0, 1.0 - max(0, rRup - 40.0) / 30.0) *
            // min(max(0, mag - 5.5) / 0.8, 1.0) * cDPP;

            return Math.Exp(r1 + r2 + r3 + r4 + r5);
        }

        private double CalcSoilNonLin()
        {
            double exp1 = Math.Exp(phi3[iper] * (Math.Min(vs30, 1130.0) - 360.0));
            double exp2 = Math.Exp(phi3[iper] * (1130.0 - 360.0));
            return phi2[iper] * (exp1 - exp2);
        }

        // Mean ground motion model -- Equation 12
        private double GetMeanForThisT(double snl, double saRef)
        {

            // Soil effect: linear response
            double sl = phi1[iper] * Math.Min(Math.Log(vs30 / 1130.0), 0.0);

            // Soil effect: nonlinear response (base passed in)
            double snl_mod = snl * Math.Log((saRef + phi4[iper]) / phi4[iper]);

            // Soil effect: sediment thickness
            double dZ1 = calcDeltaZ1();
            double rkdepth = phi5[iper] * (1.0 - Math.Exp(-dZ1 / PHI6));

            // total model
            return Math.Log(saRef) + sl + snl_mod + rkdepth;
        }

        // Center zTop on the zTop-M relation -- Equations 4, 5
        private double calcMwZtop()
        {
            double mzTop = 0.0;
            if (style == FaultStyle.REVERSE)
            {
                mzTop = (mag <= 5.849) ? 2.704 : Math.Max(2.704 - 1.226 * (mag - 5.849), 0);
            }
            else
            {
                mzTop = (mag <= 4.970) ? 2.673 : Math.Max(2.673 - 1.136 * (mag - 4.970), 0);
            }
            return mzTop * mzTop;
        }

        // -- Equation 1
        private double calcDeltaZ1()
        {
            if (Double.IsNaN(z1p0))
            {
                return 0.0;
            }
            double vsPow4 = vs30 * vs30 * vs30 * vs30;
            return z1p0 * 1000.0 - Math.Exp(-7.15 / 4.0 * Math.Log((vsPow4 + A) / B));
        }

        // Aleatory uncertainty model -- Equation 3.9
        private double GetStdForThisT(double snl, double saRef)
        {

            // Response Term - linear vs. non-linear
            double NL0 = snl * saRef / (saRef + phi4[iper]);

            // Magnitude thresholds
            double mTest = Math.Min(Math.Max(mag, 5.0), 6.5) - 5.0;

            // Inter-event Term
            double tau = tau1[iper] + (tau2[iper] - tau1[iper]) / 1.5 * mTest;

            // Intra-event term
            double sigNL0 = sig1[iper] + (sig2[iper] - sig1[iper]) / 1.5 * mTest;
            double vsTerm = vsInf ? sig3[iper] : 0.7;
            double NL0sq = (1 + NL0) * (1 + NL0);
            sigNL0 *= Math.Sqrt(vsTerm + NL0sq);

            return Math.Sqrt(tau * tau * NL0sq + sigNL0 * sigNL0);
        }
    }

}
