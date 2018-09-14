using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using GMPEs;

 

namespace Hazard
{
    public class HazardCalculation
    {

        // Instantiate the class HazardParams to save the project
        public static SetScenarioPrams ThisScenario = new SetScenarioPrams();

        public class SetScenarioPrams
        {

            //Period
            public double saPeriodParam { get; set; }

            // Others
            public double Magnitude;
            public double EpicentralDistance { get; set; }
            public double HypocenralDistance { get; set; }
            public double RuptureDistance { get; set; }
            public double JoynerBooreDistance { get; set; }
            public double RxDistance { get; set; }
            public double RyDistance { get; set; }
            public double Azimuth { get; set; }
            public double Ztor { get; set; }
            public double HypoDepth { get; set; }
            public double Dip { get; set; }
            public double Width { get; set; }
            public double Z1p0 { get; set; }
            public double Z2p5 { get; set; }
            public double VsThirty { get; set; }
            public bool IsInferred { get; set; }
            public bool HangingWallFlag { get; set; }


            public SiteType SiteType { get; set; }
            public MagnitudeType MagnitudeType { get; set; }
            public FaultStyle FaultStyle { get; set; }



        }


    } // class HazardCalculation

}// namespace Hazard


