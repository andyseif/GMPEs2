using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Hazard;

//Abrahamson, Silva & Kamai (2014) GMPE

namespace GMPEs
{
    class ASK2014_AttenRel
    {
        public readonly string SHORT_NAME = "ASK2014";
        
        // coefficients and constants:
        private static double A3 = 0.275;
        private static double A4 = -0.1;
        private static double A5 = -0.41;
        private static double M2 = 5.0;
        private static double N = 1.5;
        private static double C4 = 4.5;
        private static double A = Math.Pow(610, 4);
        private static double B = Math.Pow(1360, 4) + A; 
        private static double VS_RK = 1180.0;
        private static double A2_HW = 0.2;
        private static double H1 = 0.25;
        private static double H2 = 1.5;
        private static double H3 = -0.75;
        private static double PHI_AMP_SQ = 0.16;

        private static double[] pd = { 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 6, 7.5, 10.0 };

        private static double[] a1 = { 0.587, 0.598, 0.602, 0.707, 0.973, 1.169, 1.442, 1.637, 1.701, 1.712, 1.662, 1.571, 1.299, 1.043, 0.665, 0.329, -0.06, -0.299, -0.562, -0.875, -1.303, -1.928 };
        private static double[] a2 = { -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.765, -0.711, -0.634, -0.529 };
        private static double[] a6 = { 2.1541, 2.1461, 2.1566, 2.0845, 2.0285, 2.0408, 2.1208, 2.2241, 2.3124, 2.3383, 2.4688, 2.5586, 2.6821, 2.763, 2.8355, 2.8973, 2.9061, 2.8888, 2.8984, 2.8955, 2.87, 2.8431 };
        private static double[] a8 = { -0.015, -0.015, -0.015, -0.015, -0.015, -0.015, -0.022, -0.03, -0.038, -0.045, -0.055, -0.065, -0.095, -0.11, -0.124, -0.138, -0.172, -0.197, -0.218, -0.235, -0.255, -0.285 };
        private static double[] a10 = { 1.735, 1.718, 1.615, 1.358, 1.258, 1.31, 1.66, 2.22, 2.77, 3.25, 3.99, 4.45, 4.75, 4.3, 2.6, 0.55, -0.95, -0.95, -0.93, -0.91, -0.87, -0.8 };
        private static double[] a12 = { -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.2, -0.2, -0.2 };
        private static double[] a13 = { 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.58, 0.56, 0.53, 0.5, 0.42, 0.35, 0.2, 0, 0, 0, 0, 0 };
        private static double[] a15 = { 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.03, 0.92, 0.84, 0.68, 0.57, 0.42, 0.31, 0.16, 0.05, -0.04, -0.11, -0.19, -0.3 };
        private static double[] a17 = { -0.0072, -0.0073, -0.0075, -0.008, -0.0089, -0.0095, -0.0095, -0.0086, -0.0074, -0.0064, -0.0043, -0.0032, -0.0025, -0.0025, -0.0022, -0.0019, -0.0015, -0.001, -0.001, -0.001, -0.001, -0.001 };
        private static double[] a43 = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.14, 0.17, 0.22, 0.26, 0.34, 0.41, 0.51, 0.55, 0.49, 0.42 };
        private static double[] a44 = { 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.07, 0.1, 0.14, 0.17, 0.21, 0.25, 0.3, 0.32, 0.32, 0.32, 0.275, 0.22 };
        private static double[] a45 = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.03, 0.06, 0.1, 0.14, 0.17, 0.2, 0.22, 0.23, 0.23, 0.22, 0.2, 0.17, 0.14 };
        private static double[] a46 = { -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.03, 0, 0.03, 0.06, 0.09, 0.13, 0.14, 0.16, 0.16, 0.16, 0.14, 0.13, 0.1, 0.09, 0.08 };
        private static double[] b = { -1.47, -1.459, -1.39, -1.219, -1.152, -1.23, -1.587, -2.012, -2.411, -2.757, -3.278, -3.599, -3.8, -3.5, -2.4, -1, 0, 0, 0, 0, 0, 0 };
        private static double[] c = { 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4 };
        private static double[] s1e = { 0.754, 0.76, 0.781, 0.81, 0.81, 0.81, 0.801, 0.789, 0.77, 0.74, 0.699, 0.676, 0.631, 0.609, 0.578, 0.555, 0.548, 0.527, 0.505, 0.477, 0.457, 0.429 };
        private static double[] s2e = { 0.52, 0.52, 0.52, 0.53, 0.54, 0.55, 0.56, 0.565, 0.57, 0.58, 0.59, 0.6, 0.615, 0.63, 0.64, 0.65, 0.64, 0.63, 0.63, 0.63, 0.63, 0.63 };
        private static double[] s3 = { 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47 };
        private static double[] s4 = { 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36 };
        private static double[] s1m = { 0.741, 0.747, 0.769, 0.798, 0.798, 0.795, 0.773, 0.753, 0.729, 0.693, 0.644, 0.616, 0.566, 0.541, 0.506, 0.48, 0.472, 0.447, 0.425, 0.395, 0.378, 0.359 };
        private static double[] s2m = { 0.501, 0.501, 0.501, 0.512, 0.522, 0.527, 0.519, 0.514, 0.513, 0.519, 0.524, 0.532, 0.548, 0.565, 0.576, 0.587, 0.576, 0.565, 0.568, 0.571, 0.575, 0.585 };
        private static double[] M1 = { 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.82, 6.92, 7.0, 7.06, 7.145, 7.25 };
        private static double[] Vlin = { 660, 680, 770, 915, 960, 910, 740, 590, 495, 430, 360, 340, 330, 330, 330, 330, 330, 330, 330, 330, 330, 330 };

        // initialize dictionary of period indices. double is period, int is index.
        private Dictionary<double, int> indexFromPerHashMap = new Dictionary<double, int> { };

        private int iper;
        private double period, rRup, rJB, rX, mag, dip, zTop, z1p0, vs30, width;
        private bool isInferred;
        private SiteType siteType;
        private FaultStyle style;

        public void setParamDefaults()
        {
            period = HazardCalculation.ThisScenario.saPeriodParam;
            mag = HazardCalculation.ThisScenario.Magnitude;
            rRup = HazardCalculation.ThisScenario.RuptureDistance;
            rJB = HazardCalculation.ThisScenario.JoynerBooreDistance;
            rX = HazardCalculation.ThisScenario.RxDistance;
            siteType = HazardCalculation.ThisScenario.SiteType;
            style = HazardCalculation.ThisScenario.FaultStyle;
            dip = HazardCalculation.ThisScenario.Dip;
            width = HazardCalculation.ThisScenario.Width;
            zTop = HazardCalculation.ThisScenario.Ztop;
            z1p0 = HazardCalculation.ThisScenario.Z1p0;
            vs30 = HazardCalculation.ThisScenario.VsThirty;
            isInferred = HazardCalculation.ThisScenario.IsInferred;
        }


        public ASK2014_AttenRel()
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
            if (indexFromPerHashMap.ContainsKey(period))
            {
                // Get median directly from GMPE
                setCoeffIndex(period);
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
                    while ((perList[ind] < period) && (ind < perList.Count())) {
                        ind++;
                    }

                    // create log period and log mean vectors for interpolation
                    double[] logPerVect = { Math.Log(perList[ind - 1]), Math.Log(perList[ind]) };
                    double[] meanVect = {0, 0 };
                    setCoeffIndex(perList[ind - 1]);
                    meanVect[0] = getMeanHere();
                    setCoeffIndex(perList[ind]);
                    meanVect[1] = getMeanHere();

                    //interpolate in log-log space
                    return HelperMethods.InterpFromVector(logPerVect, meanVect, Math.Log(period));
                    }

            }
            
        }

        public double getStdDev()
        {
            // Check if key (period) is directly available from GMPE
            if (indexFromPerHashMap.ContainsKey(period))
            {
                // Get median directly from GMPE
                setCoeffIndex(period);
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
                    setCoeffIndex(perList[ind - 1]);
                    sigVect[0] = getStdDevHere();
                    setCoeffIndex(perList[ind]);
                    sigVect[1] = getStdDevHere();

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
            // ****** Mean ground motion and standard deviation model ******

            // Base Model (magnitude and distance dependence for strike-slip eq)

            // Magnitude dependent taper -- Equation 4
            double c4mag = (mag > 5) ? C4 : (mag > 4) ? C4 - (C4 - 1.0) * (5.0 - mag) : 1.0;

            // -- Equation 3
            double R = Math.Sqrt(rRup * rRup + c4mag * c4mag);

            // -- Equation 2
            double MaxMwSq = (8.5 - mag) * (8.5 - mag);
            double MwM1 = mag - M1[iper];

            double f1 = a1[iper] + a17[iper] * rRup;
            if (mag > M1[iper])
            {
                f1 += A5 * MwM1 + a8[iper] * MaxMwSq + (a2[iper] + A3 * MwM1) * Math.Log(R);
            }
            else if (mag >= M2)
            {
                f1 += A4 * MwM1 + a8[iper] * MaxMwSq + (a2[iper] + A3 * MwM1) * Math.Log(R);
            }
            else
            {
                double M2M1 = M2 - M1[iper];
                double MaxM2Sq = (8.5 - M2) * (8.5 - M2);
                double MwM2 = mag - M2;
                // a7 == 0; removed a7 * MwM2 * MwM2 below
                f1 += A4 * M2M1 + a8[iper] * MaxM2Sq + a6[iper] * MwM2 + (a2[iper] + A3 * M2M1) * Math.Log(R);
            }

            // Hanging Wall Model
            double f4 = 0.0;
            // short-circuit: f4 is 0 if rJB >= 30, rX < 0, mag <= 5.5, zTop > 10
            // these switches have been removed below
            if (rJB < 30.0 && rX >= 0.0 && mag > 5.5 && zTop <= 10.0)
            {

                // ... dip taper -- Equation 11
                double T1 = (dip > 30.0) ? (90.0 - dip) / 45.0 : 1.33333333; // 60/45

                // ... mag taper -- Equation 12
                double dM = mag - 6.5;
                double T2 = (mag >= 6.5) ? 1 + A2_HW * dM : 1 + A2_HW * dM - (1 - A2_HW) * dM * dM;

                // ... rX taper -- Equation 13
                double T3 = 0.0;
                double r1 = width * Math.Cos(dip * Math.PI / 180);
                double r2 = 3.0 * r1;
                if (rX <= r1)
                {
                    double rXr1 = rX / r1;
                    T3 = H1 + H2 * rXr1 + H3 * rXr1 * rXr1;
                }
                else if (rX <= r2)
                {
                    T3 = 1 - (rX - r1) / (r2 - r1);
                }

                // ... zTop taper -- Equation 14
                double T4 = 1 - (zTop * zTop) / 100.0;

                // ... rX, rY0 taper -- Equation 15b
                double T5 = (rJB == 0.0) ? 1.0 : 1.0 - rJB / 30.0;

                // total -- Equation 10
                f4 = a13[iper] * T1 * T2 * T3 * T4 * T5;
            }

            // Depth to Rupture Top Model -- Equation 16
            double f6 = a15[iper];
            if (zTop < 20.0)
            {
                f6 *= zTop / 20.0;
            }

            // Style-of-Faulting Model -- Equations 5 & 6
            // Note: REVERSE doesn not need to be implemented as f7 always resolves
            // to 0 as a11==0; we skip f7 here
            double f78 = (style == FaultStyle.NORMAL) ? (mag > 5.0) ? a12[iper] : (mag >= 4.0) ? a12[iper] * (mag - 4.0) : 0.0
                : 0.0;

            // Soil Depth Model -- Equation 17
            double f10 = calcSoilTerm();

            // Site Response Model
            double f5 = 0.0;
            double v1 = getV1(iper); // -- Equation 9
            double vs30s = (vs30 < v1) ? vs30 : v1; // -- Equation 8

            // Site term -- Equation 7
            double saRock = 0.0; // calc Sa1180 (rock reference) if necessary
            double c_Vlin = Vlin[iper];
            double c_b = b[iper];
            double c_c = c[iper];
            if (vs30 < c_Vlin)
            {
                // soil term (f10) for Sa1180 is zero per R. Kamai's code where Z1 < 0 for Sa1180 loop
                double vs30s_rk = (VS_RK < v1) ? VS_RK : v1;
                // use this f5 form for Sa1180 Vlin is always < 1180
                double f5_rk = (a10[iper] + c_b * N) * Math.Log(vs30s_rk / c_Vlin);
                saRock = Math.Exp(f1 + f78 + f5_rk + f4 + f6);
                f5 = a10[iper] * Math.Log(vs30s / c_Vlin) - c_b * Math.Log(saRock + c_c) + c_b *
                    Math.Log(saRock + c_c * Math.Pow(vs30s / c_Vlin, N));
            }
            else
            {
                f5 = (a10[iper] + c_b * N) * Math.Log(vs30s / c_Vlin);
            }

            // total model (no aftershock f11) -- Equation 1
            return f1 + f78 + f5 + f4 + f6 + f10;
        }

        private double getStdDevHere() 
        {

            // Magnitude dependent taper -- Equation 4
            double c4mag = (mag > 5) ? C4 : (mag > 4) ? C4 - (C4 - 1.0) * (5.0 - mag) : 1.0;

            // -- Equation 3
            double R = Math.Sqrt(rRup * rRup + c4mag * c4mag);

            // -- Equation 2
            double MaxMwSq = (8.5 - mag) * (8.5 - mag);
            double MwM1 = mag - M1[iper];

            double f1 = a1[iper] + a17[iper] * rRup;
            if (mag > M1[iper])
            {
                f1 += A5 * MwM1 + a8[iper] * MaxMwSq + (a2[iper] + A3 * MwM1) * Math.Log(R);
            }
            else if (mag >= M2)
            {
                f1 += A4 * MwM1 + a8[iper] * MaxMwSq + (a2[iper] + A3 * MwM1) * Math.Log(R);
            }
            else
            {
                double M2M1 = M2 - M1[iper];
                double MaxM2Sq = (8.5 - M2) * (8.5 - M2);
                double MwM2 = mag - M2;
                // a7 == 0; removed a7 * MwM2 * MwM2 below
                f1 += A4 * M2M1 + a8[iper] * MaxM2Sq + a6[iper] * MwM2 + (a2[iper] + A3 * M2M1) * Math.Log(R);
            }

            // Hanging Wall Model
            double f4 = 0.0;
            // short-circuit: f4 is 0 if rJB >= 30, rX < 0, mag <= 5.5, zTop > 10
            // these switches have been removed below
            if (rJB < 30.0 && rX >= 0.0 && mag > 5.5 && zTop <= 10.0)
            {

                // ... dip taper -- Equation 11
                double T1 = (dip > 30.0) ? (90.0 - dip) / 45.0 : 1.33333333; // 60/45

                // ... mag taper -- Equation 12
                double dM = mag - 6.5;
                double T2 = (mag >= 6.5) ? 1 + A2_HW * dM : 1 + A2_HW * dM - (1 - A2_HW) * dM * dM;

                // ... rX taper -- Equation 13
                double T3 = 0.0;
                double r1 = width * Math.Cos(dip * Math.PI / 180);
                double r2 = 3.0 * r1;
                if (rX <= r1)
                {
                    double rXr1 = rX / r1;
                    T3 = H1 + H2 * rXr1 + H3 * rXr1 * rXr1;
                }
                else if (rX <= r2)
                {
                    T3 = 1 - (rX - r1) / (r2 - r1);
                }

                // ... zTop taper -- Equation 14
                double T4 = 1 - (zTop * zTop) / 100.0;

                // ... rX, rY0 taper -- Equation 15b
                double T5 = (rJB == 0.0) ? 1.0 : 1.0 - rJB / 30.0;

                // total -- Equation 10
                f4 = a13[iper] * T1 * T2 * T3 * T4 * T5;
            }

            // Depth to Rupture Top Model -- Equation 16
            double f6 = a15[iper];
            if (zTop < 20.0)
            {
                f6 *= zTop / 20.0;
            }

            // Style-of-Faulting Model -- Equations 5 & 6
            // Note: REVERSE doesn not need to be implemented as f7 always resolves
            // to 0 as a11==0; we skip f7 here
            double f78 = (style == FaultStyle.NORMAL) ? (mag > 5.0) ? a12[iper] : (mag >= 4.0) ? a12[iper] * (mag - 4.0) : 0.0
                : 0.0;

            // ****** Aleatory uncertainty model ******
            double saRock = 0.0;
            double c_Vlin = Vlin[iper];
            double c_b = b[iper];
            double c_c = c[iper];

            if (vs30 < c_Vlin)
            {
                // soil term (f10) for Sa1180 is zero per R. Kamai's code where Z1 < 0 for Sa1180 loop
                double v1 = getV1(iper); // -- Equation 9
                double vs30s_rk = (VS_RK < v1) ? VS_RK : v1;
                // use this f5 form for Sa1180 Vlin is always < 1180
                double f5_rk = (a10[iper] + c_b * N) * Math.Log(vs30s_rk / c_Vlin);


                saRock = Math.Exp(f1 + f78 + f5_rk + f4 + f6);
            }

            // Intra-event term -- Equation 24
            double phiAsq;
            if (isInferred) 
            {
                phiAsq = getPhiA(s1e[iper], s2e[iper]); //mag, 
            }
            else
            {
                phiAsq = getPhiA(s1m[iper], s2m[iper]); //mag, 
            }
            phiAsq *= phiAsq;

            // Inter-event term -- Equation 25
            double tauB = getTauA(); 

            // Intra-event term with site amp variability removed -- Equation 27
            double phiBsq = phiAsq - PHI_AMP_SQ;

            // Parital deriv. of ln(soil amp) w.r.t. ln(SA1180) -- Equation 30
            // saRock subject to same vs30 < Vlin test as in mean model
            double dAmp_p1 = get_dAmp(saRock) + 1.0; 

            // phi squared, with non-linear effects -- Equation 28
            double phiSq = phiBsq * dAmp_p1 * dAmp_p1 + PHI_AMP_SQ;

            // tau squared, with non-linear effects -- Equation 29
            double tau = tauB * dAmp_p1;

            // total std dev
            return Math.Sqrt(phiSq + tau * tau);
        }

        // -- Equation 9
        private double getV1(int iper)
        {
            Double T = pd[iper];
            if (Double.IsNaN(T)) // -------------------------------------------WHEN IS THIS ENCOUNTERED?
            {
                return 1500.0;
            }
            if (T >= 3.0)
            {
                return 800.0;
            }
            if (T > 0.5)
            {
                return Math.Exp(-0.35 * Math.Log(T / 0.5) + Math.Log(1500.0));
            }
            return 1500.0;
        }

        // used for interpolation in calcSoilTerm(), below
        private static double[] VS_BINS = { 150d, 250d, 400d, 700d, 1000d };

        // Soil depth model adapted from CY13 form -- Equation 17
        private double calcSoilTerm() 
        {
            // short circuit; default z1 will be the same as z1ref
            if (Double.IsNaN(z1p0))
            {
                return 0.0;
            }
            // -- Equation 18
            double vsPow4 = vs30 * vs30 * vs30 * vs30;
            double z1ref = Math.Exp(-7.67 / 4.0 * Math.Log((vsPow4 + A) / B)) / 1000.0; // km

            // double z1c = (vs30 > 500.0) ? a46 :
            // (vs30 > 300.0) ? a45 :
            // (vs30 > 200.0) ? a44 : a43;

            // new interpolation algorithm
            double[] vsCoeff = { a43[iper], a44[iper], a45[iper], a46[iper], a46[iper] };
            double z1c = HelperMethods.InterpFromVector(VS_BINS, vsCoeff, vs30);
            return z1c * Math.Log((z1p0 + 0.01) / (z1ref + 0.01));
        }

        // -- Equation 24
        private double getPhiA(double s1p, double s2p) 
        {
            return mag < 4.0 ? s1p : mag > 6.0 ? s2p : s1p + ((s2p - s1p) / 2.0) * (mag - 4.0);
        }

        // -- Equation 25
        private double getTauA() 
        {
            return mag < 5.0 ? s3[iper] : mag > 7.0 ? s4[iper] : s3[iper] + ((s4[iper] - s3[iper]) / 2.0) * (mag - 5.0);
        }

        // -- Equation 30
        private double get_dAmp(double saRockp) 
        {
            if (vs30 >= Vlin[iper])
            {
                return 0.0;
            }
            else
            {
                return (-b[iper] * saRockp) / (saRockp + c[iper]) + 
                    (b[iper] * saRockp) / (saRockp + c[iper] * Math.Pow(vs30 / Vlin[iper], N));
            }
        }

        //private double InterpFromVector(double[] xVector, double[] yVector, double xInterp)
        //{
        //    // if xInterp outside of xVector range, return first or last value of yVector
        //    if ( xInterp <= xVector.First() ) 
        //    {
        //        return yVector.First();
        //    }
        //    else if ( xInterp >= xVector.Last() )
        //    {
        //        return yVector.Last();
        //    }

        //    // otherwise, linearly interpolate
        //    double x1, x2, y1, y2;
        //    int ind = 0;

        //    // stop at first index of xVector greater than xInterp
        //    while ((xVector[ind] < xInterp) && (ind < xVector.Length))
        //    {
        //        ind++;
        //    }
        //    if (ind == xVector.Length)
        //    {
        //        return Double.NaN;
        //    }
        //    x1 = xVector[ind - 1];
        //    y1 = yVector[ind - 1];
        //    x2 = xVector[ind];
        //    y2 = yVector[ind];

        //    return y1 + (xInterp - x1) * (y2 - y1) / (x2 - x1);     
            
        //}

    }

}
