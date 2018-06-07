using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using GMPEs;
using System.IO;

namespace Hazard
{

    class Program
    {

        static void Main(string[] args)
        {

          //  HazardCalculation.ThisScenario.saPeriodParam=1;  !!!! See below 
            HazardCalculation.ThisScenario.Magnitude = 7.1;
            HazardCalculation.ThisScenario.RuptureDistance = 25.6;
            HazardCalculation.ThisScenario.JoynerBooreDistance = 25.1;
            HazardCalculation.ThisScenario.RxDistance = 25.2;
            HazardCalculation.ThisScenario.VsThirty = 320.0;
            HazardCalculation.ThisScenario.IsInferred = false;
            HazardCalculation.ThisScenario.SiteType = SiteType.HARD_ROCK;
            HazardCalculation.ThisScenario.MagnitudeType = MagnitudeType.MOMENT;
            HazardCalculation.ThisScenario.FaultStyle = FaultStyle.NORMAL;
            HazardCalculation.ThisScenario.Dip = 90;
            HazardCalculation.ThisScenario.Width = 20;
            HazardCalculation.ThisScenario.Ztop = 1.0;
            HazardCalculation.ThisScenario.Z1p0 = 1.5;

            string GMPEName = "ASK2014_AttenRel"; // "Campbell_2003_AttenRel";  GMPEName="ToroEtAl_1997_AttenRel"; //GMPEName = "AB2006_140_AttenRel";
            // double[] periods = { 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.5, 10.0 };
            
            // alternative set of periods that require interpolation
            double[] periods = { 0.015, 0.025, 0.04, 0.0625, 0.0875, 0.125, 0.175, 0.225, 0.275, 0.35, 0.45, 0.625, 0.875, 1.25, 1.75, 2.5, 3.5, 4.5, 5.5, 6.75, 8.75 };

            double[] MedianResults = new double[periods.Length];
            double[] SigmaResults = new double[periods.Length];

            Type elementType = Type.GetType("GMPEs." + GMPEName);
            dynamic GMPE = Activator.CreateInstance(elementType); // method to instantiate a class by it's name

            try
            {

                for (int i = 0; i < periods.Length; i++)
                {
                    HazardCalculation.ThisScenario.saPeriodParam = periods[i];
                    
                    var Median = GMPE.getMean();
                    var Sigma = GMPE.getStdDev();

                    MedianResults[i] = Median;
                    SigmaResults[i] = Sigma;
                }


                using (StreamWriter writer = new StreamWriter("lnY.txt"))
                {
                    string line = "Period"+" "+ (string)GMPE.SHORT_NAME;

                    writer.WriteLine(line);

                    for (int i = 0; i < periods.Length; i++)
                    {
                        line = Convert.ToString(periods[i]);

                        line = line + " " + Convert.ToString(MedianResults[i]);
                        writer.WriteLine(line);
                    }

                

                }



                using (StreamWriter writer = new StreamWriter("Sigma.txt"))
                {
                    string line = "Period" + " " + GMPE.SHORT_NAME;

                    writer.WriteLine(line);

                    for (int i = 0; i < periods.Length; i++)
                    {
                        line = Convert.ToString(periods[i]);

                        line = line + " " + Convert.ToString(SigmaResults[i]);
                        writer.WriteLine(line);
                    }

                  

                }

            }

            catch (Exception e)
            {

                Console.WriteLine(e.Message);
                Console.WriteLine("The program will now terminate.\n");
                Console.WriteLine("Press any key to quit.");
                Console.ReadKey();
                Environment.Exit(0);
            }

        }

    }
    
}
