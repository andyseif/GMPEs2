using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GMPEs
{
    public class GroundMotion
    {
        private double logMean;
        private double logStd;

        public GroundMotion()
        {
            logMean = 0;
            logStd = 0;
        }

        public GroundMotion(double mu, double sig)
        {
            logMean = mu;
            logStd = sig;
        }

        public void SetLogMean(double mean)
        {
            logMean = mean;
        }

        public void SetLogStd(double std)
        {
            logStd = std;
        }

        public double GetLogMean()
        {
            return logMean;
        }

        public double GetLogStd()
        {
            return logStd;
        }
    }
}
