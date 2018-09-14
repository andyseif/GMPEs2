using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Hazard;

// NGA West2 GMPE: Campbell and Bozorgnia (CB) 2014
// fully implemented but not fully tested
// no range checks yet implemented for parameters (mag, dist, per, etc)

namespace GMPEs
{
    class CB2014_AttenRel
    {
        public readonly string SHORT_NAME = "CB2014";

        // coefficients and constants:
        private static double H4 = 1.0;
        private static double C = 1.88;
        private static double N = 1.18;
        private static double PHI_LNAF_SQ = 0.09;

        private static double[] pd = { 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10, 0, -1 };

        private static double[] c0 = { -4.365,-4.348,-4.024,-3.479,-3.293,-3.666,-4.866,-5.411,-5.962,-6.403,-7.566,-8.379,-9.841,-11.011,-12.469,-12.969,-13.306,-14.02,-14.558,-15.509,-15.975,-4.416,-2.895 };
        private static double[] c1 = { 0.977,0.976,0.931,0.887,0.902,0.993,1.267,1.366,1.458,1.528,1.739,1.872,2.021,2.18,2.27,2.271,2.15,2.132,2.116,2.223,2.132,0.984,1.51 };
        private static double[] c2 = { 0.533,0.549,0.628,0.674,0.726,0.698,0.51,0.447,0.274,0.193,-0.02,-0.121,-0.042,-0.069,0.047,0.149,0.368,0.726,1.027,0.169,0.367,0.537,0.27 };
        private static double[] c3 = { -1.485,-1.488,-1.494,-1.388,-1.469,-1.572,-1.669,-1.75,-1.711,-1.77,-1.594,-1.577,-1.757,-1.707,-1.621,-1.512,-1.315,-1.506,-1.721,-0.756,-0.8,-1.499,-1.299 };
        private static double[] c4 = { -0.499,-0.501,-0.517,-0.615,-0.596,-0.536,-0.49,-0.451,-0.404,-0.321,-0.426,-0.44,-0.443,-0.527,-0.63,-0.768,-0.89,-0.885,-0.878,-1.077,-1.282,-0.496,-0.453 };
        private static double[] c5 = { -2.773,-2.772,-2.782,-2.791,-2.745,-2.633,-2.458,-2.421,-2.392,-2.376,-2.303,-2.296,-2.232,-2.158,-2.063,-2.104,-2.051,-1.986,-2.021,-2.179,-2.244,-2.773,-2.466 };
        private static double[] c6 = { 0.248,0.247,0.246,0.24,0.227,0.21,0.183,0.182,0.189,0.195,0.185,0.186,0.186,0.169,0.158,0.158,0.148,0.135,0.135,0.165,0.18,0.248,0.204 };
        private static double[] c7 = { 6.753,6.502,6.291,6.317,6.861,7.294,8.031,8.385,7.534,6.99,7.012,6.902,5.522,5.65,5.795,6.632,6.759,7.978,8.538,8.468,6.564,6.768,5.837 };
        private static double[] c9 = { -0.214,-0.208,-0.213,-0.244,-0.266,-0.229,-0.211,-0.163,-0.15,-0.131,-0.159,-0.153,-0.09,-0.105,-0.058,-0.028,0,0,0,0,0,-0.212,-0.168 };
        private static double[] c10 = { 0.72,0.73,0.759,0.826,0.815,0.831,0.749,0.764,0.716,0.737,0.738,0.718,0.795,0.556,0.48,0.401,0.206,0.105,0,0,0,0.72,0.305 };
        private static double[] c11 = { 1.094,1.149,1.29,1.449,1.535,1.615,1.877,2.069,2.205,2.306,2.398,2.355,1.995,1.447,0.33,-0.514,-0.848,-0.793,-0.748,-0.664,-0.576,1.09,1.713 };
        private static double[] c14 = { -0.007,-0.0167,-0.0422,-0.0663,-0.0794,-0.0294,0.0642,0.0968,0.1441,0.1597,0.141,0.1474,0.1764,0.2593,0.2881,0.3112,0.3478,0.3747,0.3382,0.3754,0.3506,-0.0064,0.106 };
        private static double[] c16 = { 0.39,0.387,0.378,0.295,0.322,0.384,0.417,0.404,0.466,0.528,0.54,0.638,0.776,0.771,0.748,0.763,0.686,0.691,0.67,0.757,0.621,0.393,0.585 };
        private static double[] c17 = { 0.0981,0.1009,0.1095,0.1226,0.1165,0.0998,0.076,0.0571,0.0437,0.0323,0.0209,0.0092,-0.0082,-0.0131,-0.0187,-0.0258,-0.0311,-0.0413,-0.0281,-0.0205,0.0009,0.0977,0.0517 };
        private static double[] c18 = { 0.0334,0.0327,0.0331,0.027,0.0288,0.0325,0.0388,0.0437,0.0463,0.0508,0.0432,0.0405,0.042,0.0426,0.038,0.0252,0.0236,0.0102,0.0034,0.005,0.0099,0.0333,0.0327 };
        private static double[] c19 = { 0.00755,0.00759,0.0079,0.00803,0.00811,0.00744,0.00716,0.00688,0.00556,0.00458,0.00401,0.00388,0.0042,0.00409,0.00424,0.00448,0.00345,0.00603,0.00805,0.0028,0.00458,0.00757,0.00613 };
        private static double[] c20 = { -0.0055,-0.0055,-0.0057,-0.0063,-0.007,-0.0073,-0.0069,-0.006,-0.0055,-0.0049,-0.0037,-0.0027,-0.0016,-0.0006,0,0,0,0,0,0,0,-0.0055,-0.0017 };
        //private static double[] Dc20_CA = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
        //private static double[] Dc20_JP = { -0.0035,-0.0035,-0.0034,-0.0037,-0.0037,-0.0034,-0.003,-0.0031,-0.0033,-0.0035,-0.0034,-0.0034,-0.0032,-0.003,-0.0019,-0.0005,0,0,0,0,0,-0.0035,-0.0006 };
        //private static double[] Dc20_CH = { 0.0036,0.0036,0.0037,0.004,0.0039,0.0042,0.0042,0.0041,0.0036,0.0031,0.0028,0.0025,0.0016,0.0006,0,0,0,0,0,0,0,0.0036,0.0017 };
        private static double[] a2 = { 0.168,0.166,0.167,0.173,0.198,0.174,0.198,0.204,0.185,0.164,0.16,0.184,0.216,0.596,0.596,0.596,0.596,0.596,0.596,0.596,0.596,0.167,0.596 };
        private static double[] h1 = { 0.242,0.244,0.246,0.251,0.26,0.259,0.254,0.237,0.206,0.21,0.226,0.217,0.154,0.117,0.117,0.117,0.117,0.117,0.117,0.117,0.117,0.241,0.117 };
        private static double[] h2 = { 1.471,1.467,1.467,1.449,1.435,1.449,1.461,1.484,1.581,1.586,1.544,1.554,1.626,1.616,1.616,1.616,1.616,1.616,1.616,1.616,1.616,1.474,1.616 };
        private static double[] h3 = { -0.714,-0.711,-0.713,-0.701,-0.695,-0.708,-0.715,-0.721,-0.787,-0.795,-0.77,-0.77,-0.78,-0.733,-0.733,-0.733,-0.733,-0.733,-0.733,-0.733,-0.733,-0.715,-0.733 };
        private static double[] h5 = { -0.336,-0.339,-0.338,-0.338,-0.347,-0.391,-0.449,-0.393,-0.339,-0.447,-0.525,-0.407,-0.371,-0.128,-0.128,-0.128,-0.128,-0.128,-0.128,-0.128,-0.128,-0.337,-0.128 };
        private static double[] h6 = { -0.27,-0.263,-0.259,-0.263,-0.219,-0.201,-0.099,-0.198,-0.21,-0.121,-0.086,-0.281,-0.285,-0.756,-0.756,-0.756,-0.756,-0.756,-0.756,-0.756,-0.756,-0.27,-0.756 };
        private static double[] k1 = { 865,865,908,1054,1086,1032,878,748,654,587,503,457,410,400,400,400,400,400,400,400,400,865,400 };
        private static double[] k2 = { -1.186,-1.219,-1.273,-1.346,-1.471,-1.624,-1.931,-2.188,-2.381,-2.518,-2.657,-2.669,-2.401,-1.955,-1.025,-0.299,0,0,0,0,0,-1.186,-1.955 };
        private static double[] k3 = { 1.839,1.84,1.841,1.843,1.845,1.847,1.852,1.856,1.861,1.865,1.874,1.883,1.906,1.929,1.974,2.019,2.11,2.2,2.291,2.517,2.744,1.839,1.929 };
        private static double[] phi1 = { 0.734,0.738,0.747,0.777,0.782,0.769,0.769,0.761,0.744,0.727,0.69,0.663,0.606,0.579,0.541,0.529,0.527,0.521,0.502,0.457,0.441,0.734,0.655 };
        private static double[] phi2 = { 0.492,0.496,0.503,0.52,0.535,0.543,0.543,0.552,0.545,0.568,0.593,0.611,0.633,0.628,0.603,0.588,0.578,0.559,0.551,0.546,0.543,0.492,0.494 };
        private static double[] tau1 = { 0.404,0.417,0.446,0.508,0.504,0.445,0.382,0.339,0.34,0.34,0.356,0.379,0.43,0.47,0.497,0.499,0.5,0.543,0.534,0.523,0.466,0.409,0.317 };
        private static double[] tau2 = { 0.325,0.326,0.344,0.377,0.418,0.426,0.387,0.338,0.316,0.3,0.264,0.263,0.326,0.353,0.399,0.4,0.417,0.393,0.421,0.438,0.438,0.322,0.297 };
        //private static double[] phiC = { 0.166,0.166,0.165,0.162,0.158,0.17,0.18,0.186,0.191,0.198,0.206,0.208,0.221,0.225,0.222,0.226,0.229,0.237,0.237,0.271,0.29,0.166,0.19 };
        private static double[] rho = { 1,0.998,0.986,0.938,0.887,0.87,0.876,0.87,0.85,0.819,0.743,0.684,0.562,0.467,0.364,0.298,0.234,0.202,0.184,0.176,0.154,1,0.684 };

        // initialize dictionary of period indices. double is period, int is index.
        private Dictionary<double, int> indexFromPerHashMap = new Dictionary<double, int> { };

        private int iper;
        private double period, rX, rRup, rJB, mag, zHypo, zTop, z2p5, vs30, dip, width;
        private FaultStyle style;

        public void setParamDefaults()
        {
            period = HazardCalculation.ThisScenario.saPeriodParam;
            mag = HazardCalculation.ThisScenario.Magnitude;
            rX = HazardCalculation.ThisScenario.RxDistance;
            rRup = HazardCalculation.ThisScenario.RuptureDistance;
            rJB = HazardCalculation.ThisScenario.JoynerBooreDistance;
            vs30 = HazardCalculation.ThisScenario.VsThirty;
            zTop = HazardCalculation.ThisScenario.Ztor;
            z2p5 = HazardCalculation.ThisScenario.Z2p5;
            zHypo = HazardCalculation.ThisScenario.HypoDepth;
            style = HazardCalculation.ThisScenario.FaultStyle;
            dip = HazardCalculation.ThisScenario.Dip;
            width = HazardCalculation.ThisScenario.Width;
        }

        public CB2014_AttenRel()
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
            // GroundMotion newGm = new GroundMotion(); // initializes logmean and logstd to zero

            // Check if key (period) is directly available from GMPE
            period = HazardCalculation.ThisScenario.saPeriodParam;
            if (indexFromPerHashMap.ContainsKey(period))
            {
                // Get median directly from GMPE
                //setCoeffIndex(period);
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
       

        // Get index 
        private void setCoeffIndex(double perKey)
        {
            iper = indexFromPerHashMap[perKey];
        }

        private double getMeanHere()
        {
            // log mean model
            setCoeffIndex(period);

            // Magnitude term -- Equation 2
            double Fmag = c0[iper] + c1[iper] * mag;
            if (mag > 6.5)
            {
                Fmag += c2[iper] * (mag - 4.5) + c3[iper] * (mag - 5.5) + c4[iper] * (mag - 6.5);
            }
            else if (mag > 5.5)
            {
                Fmag += c2[iper] * (mag - 4.5) + c3[iper] * (mag - 5.5);
            }
            else if (mag > 4.5)
            {
                Fmag += c2[iper] * (mag - 4.5);
            }

            // Distance term -- Equation 3
            double r = Math.Sqrt(rRup * rRup + c7[iper] * c7[iper]); 
            double Fr = (c5[iper] + c6[iper] * mag) * Math.Log(r);

            // Style-of-Faulting term -- Equations 4, 5, 6
            // c8 is always 0 so REVERSE switch has been removed
            double Fflt = 0.0;
            if (style == FaultStyle.NORMAL && mag > 4.5)
            {
                Fflt = c9[iper];
                if (mag <= 5.5)
                {
                    Fflt *= (mag - 4.5);
                }
            }

            // Hanging-Wall term
            double Fhw = 0.0;
            // short-circuit: f4 is 0 if rX < 0, mag <= 5.5, zTop > 16.66
            // these switches have been removed below

            if (rX >= 0.0 && mag > 5.5 && zTop <= 16.66) { // short-circuit

                // Jennifer Donahue's HW Model plus CB08 distance taper
                // -- Equations 9, 10, 11 & 12
                double r1 = width * Math.Cos(dip* Math.PI/180.0);
                double r2 = 62.0 * mag - 350.0;
                double rXr1 = rX / r1;
                double rXr2r1 = (rX - r1) / (r2 - r1);
                double f1_rX = h1[iper] + h2[iper] * rXr1 + h3[iper] * (rXr1 * rXr1);
                double f2_rX = H4 + h5[iper] * (rXr2r1) + h6[iper] * rXr2r1 * rXr2r1;

                // ... rX -- Equation 8
                double Fhw_rX = (rX >= r1) ? Math.Max(f2_rX, 0.0) : f1_rX;

                // ... rRup -- Equation 13
                double Fhw_rRup = (rRup == 0.0) ? 1.0 : (rRup - rJB) / rRup;

                // ... magnitude -- Equation 14
                double Fhw_m = 1.0 + a2[iper] * (mag - 6.5);
                if (mag <= 6.5)
                {
                    Fhw_m *= (mag - 5.5);
                }

                // ... depth -- Equation 15
                double Fhw_z = 1.0 - 0.06 * zTop;

                // ... dip -- Equation 16
                double Fhw_d = (90.0 - dip) / 45.0;

                // ... total -- Equation 7
                Fhw = c10[iper] * Fhw_rX * Fhw_rRup * Fhw_m * Fhw_z * Fhw_d;
            }

            // Equation 18
            double pgaRock = Math.Exp(getPgaRock());
            double vsk1 = vs30 / k1[iper];
            double Fsite = (vs30 <= k1[iper]) ? c11[iper] * Math.Log(vsk1) + k2[iper] * (Math.Log(pgaRock + C * Math.Pow(vsk1, N)) - Math.Log(pgaRock + C))
                : (c11[iper] + k2[iper] * N) * Math.Log(vsk1);

            // Basin Response term -- Equation 20
            double Fsed = basinResponseTerm(iper, vs30, z2p5);

            // Hypocentral Depth term -- Equations 21, 22, 23
            double Fhyp = (zHypo <= 7.0) ? 0.0 : (zHypo <= 20.0) ? zHypo - 7.0 : 13.0;
            if (mag <= 5.5)
            {
                Fhyp *= c17[iper];
            }
            else if (mag <= 6.5)
            {
                Fhyp *= (c17[iper] + (c18[iper] - c17[iper]) * (mag - 5.5));
            }
            else
            {
                Fhyp *= c18[iper];
            }

            // Fault Dip term -- Equation 24
            double Fdip = (mag > 5.5) ? 0.0 : (mag > 4.5) ? c19[iper] * (5.5 - mag) * dip : c19[iper] * dip;

            // Anelastic Attenuation term -- Equation 25
            double Fatn = (rRup > 80.0) ? c20[iper] * (rRup - 80.0) : 0.0;

            // total model -- Equation 1
            return Fmag + Fr + Fflt + Fhw + Fsite + Fsed + Fhyp + Fdip + Fatn;
        }

        private double basinResponseTerm(int thisInd, double thisVs30, double thisZ2p5)
        {
            if (Double.IsNaN(thisZ2p5))
            {
                thisZ2p5 = Math.Exp(7.089 - 1.144 * Math.Log(thisVs30));
            }
            if (thisZ2p5 <= 1.0)
            {
                return c14[thisInd] * (thisZ2p5 - 1.0);
            }
            else if (thisZ2p5 > 3.0)
            {
                return c16[thisInd] * k3[thisInd] * Math.Exp(-0.75) * (1.0 - Math.Exp(-0.25 * (thisZ2p5 - 3.0)));
            }
            else
            {
                return 0.0;
            }
        }

        // Aleatory uncertainty model
        private double getStdDevHere()
        {
            setCoeffIndex(period);

            double pgaRock = Math.Exp(getPgaRock());
            int ipga = indexFromPerHashMap[0.0];

            // -- Equation 31
            double vsk1 = vs30 / k1[iper];
            double alpha = (vs30 < k1[iper]) ? k2[iper] * pgaRock *
                (1.0 / (pgaRock + C * Math.Pow(vsk1, N)) - 1.0 / (pgaRock + C)) : 0.0;

            // Magnitude dependence -- Equations 29 & 30
            double tau_lnYB, tau_lnPGAB, phi_lnY, phi_lnPGAB;
            if (mag <= 4.5)
            {
                tau_lnYB = tau1[iper];
                phi_lnY = phi1[iper];
                tau_lnPGAB = tau1[ipga];
                phi_lnPGAB = phi1[ipga];
            }
            else if (mag < 5.5)
            {
                tau_lnYB = stdMagDep(tau1[iper], tau2[iper]);
                phi_lnY = stdMagDep(phi1[iper], phi2[iper]);
                tau_lnPGAB = stdMagDep(tau1[ipga], tau2[ipga]);
                phi_lnPGAB = stdMagDep(phi1[ipga], phi2[ipga]);
            }
            else
            {
                tau_lnYB = tau2[iper];
                phi_lnY = phi2[iper];
                tau_lnPGAB = tau2[ipga];
                phi_lnPGAB = phi2[ipga];
            }

            // intra-event std dev -- Equation 27
            double alphaTau = alpha * tau_lnPGAB;
            double tauSq = tau_lnYB * tau_lnYB + alphaTau * alphaTau +
                2.0 * alpha * rho[iper] * tau_lnYB * tau_lnPGAB;

            // inter-event std dev -- Equation 28
            double phi_lnYB = Math.Sqrt(phi_lnY * phi_lnY - PHI_LNAF_SQ);
            phi_lnPGAB = Math.Sqrt(phi_lnPGAB * phi_lnPGAB - PHI_LNAF_SQ);
            double aPhi_lnPGAB = alpha * phi_lnPGAB;

            // phi_lnaf terms in eqn. 30 cancel when expanded leaving phi_lnY only
            double phiSq = phi_lnY * phi_lnY + aPhi_lnPGAB * aPhi_lnPGAB +
                2.0 * rho[iper] * phi_lnYB * aPhi_lnPGAB;

            // total model -- Equation 32
            return Math.Sqrt(phiSq + tauSq);
        }

        // Source/Event Term -- Equation 2
        private double stdMagDep(double lo, double hi)
        {
            return hi + (lo - hi) * (5.5 - mag);
        }

        private double getPgaRock()
        {
            // PGA log mean model
            int indPga = indexFromPerHashMap[0.0];
            double vs30_Rock = 1100.0;
            double z2p5_Rock = 0.398;

            // Magnitude term -- Equation 2
            double Fmag = c0[indPga] + c1[indPga] * mag;
            if (mag > 6.5)
            {
                Fmag += c2[indPga] * (mag - 4.5) + c3[indPga] * (mag - 5.5) + c4[indPga] * (mag - 6.5);
            }
            else if (mag > 5.5)
            {
                Fmag += c2[indPga] * (mag - 4.5) + c3[indPga] * (mag - 5.5);
            }
            else if (mag > 4.5)
            {
                Fmag += c2[indPga] * (mag - 4.5);
            }

            // Distance term -- Equation 3
            double r = Math.Sqrt(rRup * rRup + c7[indPga] * c7[indPga]);
            double Fr = (c5[indPga] + c6[indPga] * mag) * Math.Log(r);

            // Style-of-Faulting term -- Equations 4, 5, 6
            // c8 is always 0 so REVERSE switch has been removed
            double Fflt = 0.0;
            if (style == FaultStyle.NORMAL && mag > 4.5)
            {
                Fflt = c9[indPga];
                if (mag <= 5.5)
                {
                    Fflt *= (mag - 4.5);
                }
            }

            // Hanging-Wall term
            double Fhw = 0.0;
            // short-circuit: f4 is 0 if rX < 0, mag <= 5.5, zTop > 16.66
            // these switches have been removed below

            if (rX >= 0.0 && mag > 5.5 && zTop <= 16.66)
            { // short-circuit

                // Jennifer Donahue's HW Model plus CB08 distance taper
                // -- Equations 9, 10, 11 & 12
                double r1 = width * Math.Cos(dip * Math.PI / 180.0);
                double r2 = 62.0 * mag - 350.0;
                double rXr1 = rX / r1;
                double rXr2r1 = (rX - r1) / (r2 - r1);
                double f1_rX = h1[indPga] + h2[indPga] * rXr1 + h3[indPga] * (rXr1 * rXr1);
                double f2_rX = H4 + h5[indPga] * (rXr2r1) + h6[indPga] * rXr2r1 * rXr2r1;

                // ... rX -- Equation 8
                double Fhw_rX = (rX >= r1) ? Math.Max(f2_rX, 0.0) : f1_rX;

                // ... rRup -- Equation 13
                double Fhw_rRup = (rRup == 0.0) ? 1.0 : (rRup - rJB) / rRup;

                // ... magnitude -- Equation 14
                double Fhw_m = 1.0 + a2[indPga] * (mag - 6.5);
                if (mag <= 6.5)
                {
                    Fhw_m *= (mag - 5.5);
                }

                // ... depth -- Equation 15
                double Fhw_z = 1.0 - 0.06 * zTop;

                // ... dip -- Equation 16
                double Fhw_d = (90.0 - dip) / 45.0;

                // ... total -- Equation 7
                Fhw = c10[indPga] * Fhw_rX * Fhw_rRup * Fhw_m * Fhw_z * Fhw_d;
            }

            // Shallow Site Response term - pgaRock term is computed through an
            // initial call to this method with vs30=1100; 1100 is higher than any
            // k1 value so else condition always prevails -- Equation 18
            double vsk1 = vs30_Rock / k1[indPga];
            double Fsite = (vs30_Rock <= k1[indPga]) ? c11[indPga] * Math.Log(vsk1) +
                k2[indPga] * (Math.Log(0.0 + C * Math.Pow(vsk1, N)) - Math.Log(0.0 + C))
                : (c11[indPga] + k2[indPga] * N) * Math.Log(vsk1);

            // Basin Response term -- Equation 20
            double Fsed = basinResponseTerm(indPga, vs30_Rock, z2p5_Rock);

            // Hypocentral Depth term -- Equations 21, 22, 23
            double Fhyp = (zHypo <= 7.0) ? 0.0 : (zHypo <= 20.0) ? zHypo - 7.0 : 13.0;
            if (mag <= 5.5)
            {
                Fhyp *= c17[indPga];
            }
            else if (mag <= 6.5)
            {
                Fhyp *= (c17[indPga] + (c18[indPga] - c17[indPga]) * (mag - 5.5));
            }
            else
            {
                Fhyp *= c18[indPga];
            }

            // Fault Dip term -- Equation 24
            double Fdip = (mag > 5.5) ? 0.0 : (mag > 4.5) ? c19[indPga] * (5.5 - mag) * dip : c19[indPga] * dip;

            // Anelastic Attenuation term -- Equation 25
            double Fatn = (rRup > 80.0) ? c20[indPga] * (rRup - 80.0) : 0.0;

            // total model -- Equation 1
            return Fmag + Fr + Fflt + Fhw + Fsite + Fsed + Fhyp + Fdip + Fatn;
        }

    }

}
