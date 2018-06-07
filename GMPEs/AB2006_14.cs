using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Hazard;

namespace GMPEs
{
    public class AB2006_140_AttenRel
    {
        public readonly  string SHORT_NAME = "AB2006_140";

        // coefficients and constants:
        private double[] c1, c2, c3, c4, c5, c6, c7, c8, c9, c10;
        private double[] clamp = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 3.0, 3.0, 3.0, 460.0 };
        // rounded .3968 to 0.4 s for one element of abper. SH june 30 2008
        // AZ rounded 0.0996 to 0.1
        // private double[] pd = { 5.0000, 4.0000, 3.1250, 2.5000, 2.0000, 1.5873, 1.2500, 1.0000, 0.7937, 0.6289, 0.5000, 0.40, 0.3155, 0.2506, 0.2000, 0.1580, 0.1255, 0.0996, 0.0791, 0.0629, 0.0499, 0.0396, 0.0315, 0.0250, 0.0000, -1.000 };

        private double[] pd = { 5.0000, 4.0000, 3.1250, 2.5000, 2.0000, 1.5873, 1.2500, 1.0000, 0.7937, 0.6289, 0.5000, 0.40, 0.3155, 0.2506, 0.2000, 0.1580, 0.1255, 0.1, 0.0791, 0.0629, 0.0499, 0.0396, 0.0315, 0.0250, 0.0000, -1.000 };

        private double[] Fr = { 2.00e-1, 2.50e-1, 3.20e-1, 4.00e-1, 5.00e-1, 6.30e-1, 8.00e-1, 1.00, 1.26, 1.59, 2.00, 2.52, 3.17, 3.99, 5.03, 6.33, 7.97, 1.00e1, 1.26e1, 1.59e1, 2.00e1, 2.52e1, 3.18e1, 4.00e1, 0.00, -1.00 };
        // hr coeffs from AB06 Table 6
        private double[] c1hr = { -5.41, -5.79, -6.04, -6.17, -6.18, -6.04, -5.72, -5.27, -4.60, -3.92, -3.22, -2.44, -1.72, -1.12, -6.15e-1, -1.46e-1, 2.14e-1, 4.80e-1, 6.91e-1, 9.11e-1, 1.11, 1.26, 1.44, 1.52, 9.07e-1, -1.44 };
        private double[] c2hr = { 1.71, 1.92, 2.08, 2.21, 2.30, 2.34, 2.32, 2.26, 2.13, 1.99, 1.83, 1.65, 1.48, 1.34, 1.23, 1.12, 1.05, 1.02, 9.97e-1, 9.80e-1, 9.72e-1, 9.68e-1, 9.59e-1, 9.60e-1, 9.83e-1, 9.91e-1 };
        private double[] c3hr = { -9.01e-2, -1.07e-1, -1.22e-1, -1.35e-1, -1.44e-1, -1.50e-1, -1.51e-1, -1.48e-1, -1.41e-1, -1.31e-1, -1.20e-1, -1.08e-1, -9.74e-2, -8.72e-2, -7.89e-2, -7.14e-2, -6.66e-2, -6.40e-2, -6.28e-2, -6.21e-2, -6.20e-2, -6.23e-2, -6.28e-2, -6.35e-2, -6.60e-2, -5.85e-2 };
        private double[] c4hr = { -2.54, -2.44, -2.37, -2.30, -2.22, -2.16, -2.10, -2.07, -2.06, -2.05, -2.02, -2.05, -2.08, -2.08, -2.09, -2.12, -2.15, -2.20, -2.26, -2.36, -2.47, -2.58, -2.71, -2.81, -2.70, -2.70 };
        private double[] c5hr = { 2.27e-1, 2.11e-1, 2.00e-1, 1.90e-1, 1.77e-1, 1.66e-1, 1.57e-1, 1.50e-1, 1.47e-1, 1.42e-1, 1.34e-1, 1.36e-1, 1.38e-1, 1.35e-1, 1.31e-1, 1.30e-1, 1.30e-1, 1.27e-1, 1.25e-1, 1.26e-1, 1.28e-1, 1.32e-1, 1.40e-1, 1.46e-1, 1.59e-1, 2.16e-1 };
        private double[] c6hr = { -1.27, -1.16, -1.07, -9.86e-1, -9.37e-1, -8.70e-1, -8.20e-1, -8.13e-1, -7.97e-1, -7.82e-1, -8.13e-1, -8.43e-1, -8.89e-1, -9.71e-1, -1.12, -1.30, -1.61, -2.01, -2.49, -2.97, -3.39, -3.64, -3.73, -3.65, -2.80, -2.44 };
        private double[] c7hr = { 1.16e-1, 1.02e-1, 8.95e-2, 7.86e-2, 7.07e-2, 6.05e-2, 5.19e-2, 4.67e-2, 4.35e-2, 4.30e-2, 4.44e-2, 4.48e-2, 4.87e-2, 5.63e-2, 6.79e-2, 8.31e-2, 1.05e-1, 1.33e-1, 1.64e-1, 1.91e-1, 2.14e-1, 2.28e-1, 2.34e-1, 2.36e-1, 2.12e-1, 2.66e-1 };
        private double[] c8hr = { 9.79e-1, 1.01, 1.00, 9.68e-1, 9.52e-1, 9.21e-1, 8.56e-1, 8.26e-1, 7.75e-1, 7.88e-1, 8.84e-1, 7.39e-1, 6.10e-1, 6.14e-1, 6.06e-1, 5.62e-1, 4.27e-1, 3.37e-1, 2.14e-1, 1.07e-1, -1.39e-1, -3.51e-1, -5.43e-1, -6.54e-1, -3.01e-1, 8.48e-2 };
        private double[] c9hr = { -1.77e-1, -1.82e-1, -1.80e-1, -1.77e-1, -1.77e-1, -1.73e-1, -1.66e-1, -1.62e-1, -1.56e-1, -1.59e-1, -1.75e-1, -1.56e-1, -1.39e-1, -1.43e-1, -1.46e-1, -1.44e-1, -1.30e-1, -1.27e-1, -1.21e-1, -1.17e-1, -9.84e-2, -8.13e-2, -6.45e-2, -5.50e-2, -6.53e-2, -6.93e-2 };
        private double[] c10hr = { -1.76e-04, -2.01e-04, -2.31e-04, -2.82e-04, -3.22e-04, -3.75e-04, -4.33e-04, -4.86e-04, -5.79e-04, -6.95e-04, -7.70e-04, -8.51e-04, -9.54e-04, -1.06e-03, -1.13e-03, -1.18e-03, -1.15e-03, -1.05e-03, -8.47e-04, -5.79e-04, -3.17e-04, -1.23e-04, -3.23e-05, -4.85e-05, -4.48e-04, -3.73e-04 };
        // bc coeffs from AB06 Table 9
        private double[] c1bc = { -4.85, -5.26, -5.59, -5.80, -5.85, -5.75, -5.49, -5.06, -4.45, -3.75, -3.01, -2.28, -1.56, -8.76e-1, -3.06e-1, 1.19e-1, 5.36e-1, 7.82e-1, 9.67e-1, 1.11, 1.21, 1.26, 1.19, 1.05, 5.23e-1, -1.66 };
        private double[] c2bc = { 1.58, 1.79, 1.97, 2.13, 2.23, 2.29, 2.29, 2.23, 2.12, 1.97, 1.80, 1.63, 1.46, 1.29, 1.16, 1.06, 9.65e-1, 9.24e-1, 9.03e-1, 8.88e-1, 8.83e-1, 8.79e-1, 8.88e-1, 9.03e-1, 9.69e-1, 1.05 };
        private double[] c3bc = { -8.07e-2, -9.79e-2, -1.14e-1, -1.28e-1, -1.39e-1, -1.45e-1, -1.48e-1, -1.45e-1, -1.39e-1, -1.29e-1, -1.18e-1, -1.05e-1, -9.31e-2, -8.19e-2, -7.21e-2, -6.47e-2, -5.84e-2, -5.56e-2, -5.48e-2, -5.39e-2, -5.44e-2, -5.52e-2, -5.64e-2, -5.77e-2, -6.20e-2, -6.04e-2 };
        private double[] c4bc = { -2.53, -2.44, -2.33, -2.26, -2.20, -2.13, -2.08, -2.03, -2.01, -2.00, -1.98, -1.97, -1.98, -2.01, -2.04, -2.05, -2.11, -2.17, -2.25, -2.33, -2.44, -2.54, -2.58, -2.57, -2.44, -2.50 };
        private double[] c5bc = { 2.22e-1, 2.07e-1, 1.91e-1, 1.79e-1, 1.69e-1, 1.58e-1, 1.50e-1, 1.41e-1, 1.36e-1, 1.31e-1, 1.27e-1, 1.23e-1, 1.21e-1, 1.23e-1, 1.22e-1, 1.19e-1, 1.21e-1, 1.19e-1, 1.22e-1, 1.23e-1, 1.30e-1, 1.39e-1, 1.45e-1, 1.48e-1, 1.47e-1, 1.84e-1 };
        private double[] c6bc = { -1.43, -1.31, -1.20, -1.12, -1.04, -9.57e-1, -9.00e-1, -8.74e-1, -8.58e-1, -8.42e-1, -8.47e-1, -8.88e-1, -9.47e-1, -1.03, -1.15, -1.36, -1.67, -2.10, -2.53, -2.88, -3.04, -2.99, -2.84, -2.65, -2.34, -2.30 };
        private double[] c7bc = { 1.36e-1, 1.21e-1, 1.10e-1, 9.54e-2, 8.00e-2, 6.76e-2, 5.79e-2, 5.41e-2, 4.98e-2, 4.82e-2, 4.70e-2, 5.03e-2, 5.58e-2, 6.34e-2, 7.38e-2, 9.16e-2, 1.16e-1, 1.48e-1, 1.78e-1, 2.01e-1, 2.13e-1, 2.16e-1, 2.12e-1, 2.07e-1, 1.91e-1, 2.50e-1 };
        private double[] c8bc = { 6.34e-1, 7.34e-1, 8.45e-1, 8.91e-1, 8.67e-1, 8.67e-1, 8.21e-1, 7.92e-1, 7.08e-1, 6.77e-1, 6.67e-1, 6.84e-1, 6.50e-1, 5.81e-1, 5.08e-1, 5.16e-1, 3.43e-1, 2.85e-1, 1.00e-1, -3.19e-2, -2.10e-1, -3.91e-1, -4.37e-1, -4.08e-1, -8.70e-2, 1.27e-1 };
        private double[] c9bc = { -1.41e-1, -1.56e-1, -1.72e-1, -1.80e-1, -1.79e-1, -1.79e-1, -1.72e-1, -1.70e-1, -1.59e-1, -1.56e-1, -1.55e-1, -1.58e-1, -1.56e-1, -1.49e-1, -1.43e-1, -1.50e-1, -1.32e-1, -1.32e-1, -1.15e-1, -1.07e-1, -9.00e-2, -6.75e-2, -5.87e-2, -5.77e-2, -8.29e-2, -8.70e-2 };
        private double[] c10bc = { -1.61e-04, -1.96e-04, -2.45e-04, -2.60e-04, -2.86e-04, -3.43e-04, -4.07e-04, -4.89e-04, -5.75e-04, -6.76e-04, -7.68e-04, -8.59e-04, -9.55e-04, -1.05e-03, -1.14e-03, -1.18e-03, -1.13e-03, -9.90e-04, -7.72e-04, -5.48e-04, -4.15e-04, -3.88e-04, -4.33e-04, -5.12e-04, -6.30e-04, -4.27e-04 };
        private double[] bln = { -7.52e-1, -7.45e-1, -7.40e-1, -7.35e-1, -7.30e-1, -7.26e-1, -7.16e-1, -7.00e-1, -6.90e-1, -6.70e-1, -6.00e-1, -5.00e-1, -4.45e-1, -3.90e-1, -3.06e-1, -2.80e-1, -2.60e-1, -2.50e-1, -2.32e-1, -2.49e-1, -2.86e-1, -3.14e-1, -3.22e-1, -3.30e-1, -3.61e-1, -6.00e-1 };
        private double[] b1 = { -3.00e-1, -3.10e-1, -3.30e-1, -3.52e-1, -3.75e-1, -3.95e-1, -3.40e-1, -4.40e-1, -4.65e-1, -4.80e-1, -4.95e-1, -5.08e-1, -5.13e-1, -5.18e-1, -5.21e-1, -5.28e-1, -5.60e-1, -5.95e-1, -6.37e-1, -6.42e-1, -6.43e-1, -6.09e-1, -6.18e-1, -6.24e-1, -6.41e-1, -4.95e-1 };
        private double[] b2 = { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, -2.00e-03, -3.10e-2, -6.00e-2, -9.50e-2, -1.30e-1, -1.60e-1, -1.85e-1, -1.85e-1, -1.40e-1, -1.32e-1, -1.17e-1, -1.05e-1, -1.05e-1, -1.05e-1, -1.08e-1, -1.15e-1, -1.44e-1, -6.00e-2 };
        private double[] del = { 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.11 };
        private double[] m1 = { 6.0, 5.75, 5.5, 5.25, 5.0, 4.84, 4.67, 4.5, 4.34, 4.17, 4.0, 3.65, 3.3, 2.9, 2.5, 1.85, 1.15, 0.5, 0.34, 0.17, 0.0, 0.0, 0.0, 0.0, 0.5, 2.0 };
        private double[] mh = { 8.5, 8.37, 8.25, 8.12, 8.0, 7.7, 7.45, 7.2, 6.95, 6.7, 6.5, 6.37, 6.25, 6.12, 6.0, 5.84, 5.67, 5.5, 5.34, 5.17, 5.0, 5.0, 5.0, 5.0, 5.5, 5.5 };


        private Dictionary<double, int> indexFromPerHashMap = new Dictionary<double, int> { };

        private int iper;
        private double rRup, mag;
        private SiteType siteType;
        private StressDrop stressDrop;

        public void setParamDefaults()
        {
            stressDrop = StressDrop.SD_140;
            siteType = HazardCalculation.ThisScenario.SiteType;
            mag = HazardCalculation.ThisScenario.Magnitude;
            rRup = HazardCalculation.ThisScenario.RuptureDistance;
            updateCoeffs();
        }



        public AB2006_140_AttenRel()
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
            return getMean(iper, siteType, stressDrop, rRup, mag);
        }

        public double getStdDev()
        {
            return getStdDev(iper);
        }

        // precomputed fixed-value factors and terms
        private static readonly double stressfac = 0.5146;
        private static readonly double gfac = 6.8875526; // ln (980)
        private static readonly double sfac = 2.302585; // ln(10)
        private static readonly double tfac = -0.5108256; // ln(0.6)
        private static readonly double vref = 760;
        private static readonly double v1 = 180;
        private static readonly double v2 = 300;
        private static readonly double fac70 = 1.8450980; // log10(70)
        private static readonly double fac140 = 2.1461280; // log10(140)
        private static readonly double facv1 = -0.5108256; // ln(v1/v2)
        private static readonly double facv2 = -0.9295360; // ln(v2/vref)

        private double getMean(int iper, SiteType st, StressDrop sd, double rRup, double mag)
        {
            double vs30 = (st == SiteType.HARD_ROCK) ? 2000 : 760;

            // double sigma = 0.3 * sfac;
            // double sigmaf = 1.0 / SQRT_2 / sigma;

            double period = pd[iper];

            // double H1 = Doubles.max(dtor(kk),2.0);
            // double H1sq = H1*H1;

            // if (magType == LG_PHASE) mag = Utils.mblgToMw(magConvCode, mag);

            double gndm = c1[24] + c2[24] * mag + c3[24] * mag * mag; // pga
                                                                      // reference
            double gndmp = c1[iper] + c2[iper] * mag + c3[iper] * mag * mag;

            double sf2;// , ie; ie apparantly not used

            // if (ir > 2) { // 200 bar
            if (sd == StressDrop.SD_200)
            {
                double diff = Math.Max(mag - m1[iper], 0.0);
                // c sf2 is supposed to be eqn(6) of AB06 paper. Note use of Mw.
                sf2 = stressfac *
                    Math.Min(del[iper] + 0.05, 0.05 + del[iper] * diff /
                        (mh[iper] - m1[iper]));
                // ie = 2;
            }
            else
            {
                sf2 = 0.0;
                // c default stress factor, use ie=1
                // ie = 1;
            }

            // R=sqrt(rjb**2+H1sq)
            // use rRup - NOTE min 2km for gridded, 5km for faults

            // pga calculations
            double rfac = Math.Log10(rRup);
            double f0 = Math.Max(1.0 - rfac, 0.0);
            double f1 = Math.Min(rfac, fac70);
            double f2 = Math.Max(rfac - fac140, 0.0);

            double gnd, S = 0;
            // why this check; is there a NSHMP flag for negative vs?
            if (vs30 > 0)
            {
                // compute pga on rock
                gnd = gndm + (c4[24] + c5[24] * mag) * f1 +
                    (c6[24] + c7[24] * mag) * f2 + (c8[24] + c9[24] * mag) * f0 +
                    c10[24] * rRup + sf2;
                // apply stress factor before nonlinear adjustments, which occur
                // in eqn (7) of ab paper
                double bnl;

                if (vs30 <= v1)
                {
                    bnl = b1[iper];
                }
                else if (vs30 <= v2)
                {
                    bnl = (b1[iper] - b2[iper]) * Math.Log(vs30 / v2) / facv1 +
                        b1[iper];
                }
                else if (vs30 <= vref)
                {
                    bnl = b2[iper] * Math.Log(vs30 / vref) / facv2;
                }
                else
                {
                    bnl = 0.0;
                }

                double pga_bc = Math.Pow(10, gnd);

                if (st == SiteType.HARD_ROCK)
                {
                    S = 0.0;
                }
                else if (pga_bc <= 60.0)
                {
                    S = bln[iper] * Math.Log(vs30 / vref) + bnl * tfac;
                }
                else
                {
                    S = bln[iper] * Math.Log(vs30 / vref) + bnl *
                        Math.Log(pga_bc / 100.0);
                }
                // need to take alog10(exp(S)) according to eqns. 7a and 7b
                // AB2006. p. 2200 BSSA (new nov 26 2007)
                S = Math.Log10(Math.Exp(S));
            }
            gnd = gndmp + (c4[iper] + c5[iper] * mag) * f1 +
                (c6[iper] + c7[iper] * mag) * f2 + (c8[iper] + c9[iper] * mag) *
                f0 + c10[iper] * rRup + sf2 + S;
            if (iper < 26)
            {
                gnd = gnd * sfac - gfac;
            }
            else
            {
                // pgv?
                gnd = gnd * sfac;
            }



            return gnd;
        }

        private double getStdDev(int iper)
        {
            // if (magType == LG_PHASE) mag = Utils.mblgToMw(magConvCode, mag);
            return 0.3 * sfac;
        }

        private void updateCoeffs()
        {
            if (siteType == SiteType.HARD_ROCK)
            {
                c1 = c1hr;
                c2 = c2hr;
                c3 = c3hr;
                c4 = c4hr;
                c5 = c5hr;
                c6 = c6hr;
                c7 = c7hr;
                c8 = c8hr;
                c9 = c9hr;
                c10 = c10hr;
            }
            else
            {
                c1 = c1bc;
                c2 = c2bc;
                c3 = c3bc;
                c4 = c4bc;
                c5 = c5bc;
                c6 = c6bc;
                c7 = c7bc;
                c8 = c8bc;
                c9 = c9bc;
                c10 = c10bc;
            }

        }


        public enum StressDrop
        {
            SD_140,
            SD_200
        }


    }

}

