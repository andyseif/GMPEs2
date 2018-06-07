using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Hazard;
//using static GMPEs.Utils;

namespace GMPEs
{
    public class ToroEtAl_1997_AttenRel
    {
        public readonly  string SHORT_NAME = "ToroEtAl_1997_AttenRel";
        // coefficients:
        // @formatter:off
        // added 0.04 and 0.4 s coeffs july 16 2008 (NRC files)
        private double[] perx = { 0.0, 0.2, 1.0, 0.1, 0.3, 0.5, 2.0, 0.04, 0.4 };
        // tb MbLg coeffs. BC/A 2-hz Siteamp = 1.58, with BC-A coef. diff. of 0.4574.
        private double[] tb1 = { 2.489, 2.165, 0.173, 2.91, 1.7323, 1.109, -0.788, 4.0, 1.2 };
        private double[] tb1h = { 2.07, 1.6, -0.12, 2.36, 1.19, 0.652, -0.97, 3.54, 0.90 };
        private double[] tb2 = { 1.20, 1.24, 2.05, 1.23, 1.51, 1.785, 2.52, 1.19, 1.70 };
        private double[] tb3 = { 0.0, 0.0, -0.34, 0.0, -0.11, -0.2795, -0.47, 0.0, -0.26 };
        private double[] tb4 = { 1.28, 0.98, 0.90, 1.12, 0.96, 0.930, 0.93, 1.46, 0.94 };
        private double[] tb5 = { 1.23, 0.74, 0.59, 1.05, 0.6881, 0.6354, 0.6, 1.84, 0.65 };
        private double[] tb6 = { 0.0018, 0.0039, 0.0019, 0.0043, 0.0034, 0.002732, 0.0012, 0.0010, .0030 };

        // tc Mw coeffs for BC rock. 3hz BC-A is 0.5423 (BC/A siteamp is then 1.72)
        // tc Mw coeffs. 3.33 hz is log-log from the 2.5 and 5 hz values.
        private double[] tc1 = { 2.619, 2.295, 0.383, 2.92, 1.8823, 1.2887, -0.558, 4.0, 1.4 };
        private double[] tc1h = { 2.20, 1.73, 0.09, 2.37, 1.34, 0.8313, -0.740, 3.68, 1.07 };
        private double[] tc2 = { 0.81, 0.84, 1.42, 0.81, 0.964, 1.14, 1.86, 0.80, 1.05 };
        private double[] tc3 = { 0.0, 0.0, -0.2, 0.0, -0.059, -0.1244, -0.31, 0.0, -0.10 };
        private double[] tc4 = { 1.27, 0.98, 0.90, 1.1, 0.951, 0.9227, 0.92, 1.46, 0.93 };
        private double[] tc5 = { 1.16, 0.66, 0.49, 1.02, 0.601, 0.5429, 0.46, 1.77, 0.56 };
        private double[] tc6 = { 0.0021, 0.0042, 0.0023, 0.004, 0.00367, 0.00306, 0.0017, 0.0013, 0.0033 };

        private double[] tbh = { 9.3, 7.5, 6.8, 8.5, 7.35, 7.05, 7.0, 10.5, 7.2 };
        private double[] th = { 9.3, 7.5, 6.8, 8.3, 7.26, 7.027, 6.9, 10.5, 7.1 };
        private double[] clamp = { 3.0, 6.0, 0.0, 6.0, 6.0, 6.0, 0.0, 6.0, 6.0 };

        // Sigma in nat log units. Saves a divide
        // Toro : slightly larger sigma for 1 and 2 s. Toro Lg based mag has
        // larger sigma for larger M (table 3, p 50 ,srl 1997. This isn't
        // in our rendering)
        private double[] tsigma = { 0.7506, 0.7506, 0.799, 0.7506, 0.7506, 0.7506, 0.799, 0.7506, 0.7506 };
        // @formatter:on


        private Dictionary<double, int> indexFromPerHashMap = new Dictionary<double, int> { };

        private int iper;
        private double rjb, mag;
        private SiteType siteType;
        private MagnitudeType magType;


        public void setParamDefaults()
        {
            magType = HazardCalculation.ThisScenario.MagnitudeType;
            siteType = HazardCalculation.ThisScenario.SiteType;
            mag = HazardCalculation.ThisScenario.Magnitude;
            rjb = HazardCalculation.ThisScenario.RuptureDistance;
        }

        //public enum SiteType
        //{

        //    FIRM_ROCK,
        //    HARD_ROCK

        //}

        public ToroEtAl_1997_AttenRel()
        {

            for (int i = 0; i < perx.Length; i++)
            {
                indexFromPerHashMap.Add(perx[i], i);
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
            return getMean(iper, siteType, rjb, mag);
        }

        public double getStdDev()
        {
            return getStdDev(iper);
        }

        private double getMean(int iper, SiteType st, double rjb, double mag)
        {

            double period = perx[iper];

            double t1, t2, t3, t4, t5, t6;
            double thsq, t1h;

            // set coefficients based on magnitude type
            if (magType == MagnitudeType.LG_PHASE)
            {
                t1 = tb1[iper];
                t2 = tb2[iper];
                t3 = tb3[iper];
                t4 = tb4[iper];
                t5 = tb5[iper];
                t6 = tb6[iper];
                thsq = tbh[iper] * tbh[iper];
                t1h = tb1h[iper];
            }
            else
            {
                t1 = tc1[iper];
                t2 = tc2[iper];
                t3 = tc3[iper];
                t4 = tc4[iper];
                t5 = tc5[iper];
                t6 = tc6[iper];
                thsq = th[iper] * th[iper];
                t1h = tc1h[iper];
            }

            // magnitude correction
            // With toro model, you change the coefficients appropriate to the
            // magnitude.
            // New, Nov 2006: the finite-fault correction, affects the
            // fictitious depth or bending point;
            // from Toro Paducah paper. Mod. Dec 2007, mblg to Mw for the
            // correction.
            // * if(mlg) then

            /* AZ commnted out
            double mCorr;
            if (magType == MagnitudeType.LG_PHASE)
            {
                double mag1 = Utils.mblgToMw(FaultCode.M_CONV_J, mag);
                double cor1 = Math.Exp(-1.25 + 0.227 * mag1);
                double mag2 = Utils.mblgToMw(FaultCode.M_CONV_AB, mag);
                double cor2 = Math.Exp(-1.25 + 0.227 * mag2);
                mCorr = Math.Sqrt(cor1 * cor2); // geo mean
            }
            else
            {
                mCorr = Math.Exp(-1.25 + 0.227 * mag);
            }
            */

            double corsq = 1; //AZ chnaged to 1. originaly was: mCorr * mCorr;
            double dist0 = rjb;
            double dist = Math.Sqrt(dist0 * dist0 + thsq * corsq);

            // default to SOFT_ROCK values
            double gnd0 = (st == SiteType.HARD_ROCK) ? t1h : t1;
            double gndm = gnd0 + t2 * (mag - 6.0) + t3 *
                ((mag - 6.0) * (mag - 6.0));
            double gnd = gndm - t4 * Math.Log(dist) - t6 * dist;

            double factor = Math.Log(dist / 100.0);
            if (factor > 0) gnd = gnd - (t5 - t4) * factor;



            return gnd;
        }

        private double getStdDev(int iper)
        {
            return tsigma[iper];
        }



    }
}
