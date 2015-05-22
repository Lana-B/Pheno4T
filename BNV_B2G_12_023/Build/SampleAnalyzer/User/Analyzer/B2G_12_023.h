#ifndef analysis_B2G_12_023_h
#define analysis_B2G_12_023_h

#include "SampleAnalyzer/Process/Analyzer/AnalyzerBase.h"

namespace MA5
{
	class B2G_12_023 : public AnalyzerBase
	{
		INIT_ANALYSIS(B2G_12_023,"B2G_12_023")

		public:
			virtual bool Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters);
			virtual void Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files);
			virtual bool Execute(SampleFormat& sample, const EventFormat& event);


        double dFunctionBeta(double dBeta, string lepton); //function for calculating MET efficiency
		double dFunctionEfficiencies(double dBeta, string lepton, string TorB, string channel)
        int dCounterPassedEvents;
        double dCounterSelectionEff;
        
		private:
			double dSigma_tt = 22.2; // ±1.5pb
			double dSigma_tW = 246; // ± 12pb

			double dEpsilon_X_SM_SM;
			double dEpsilon_X_BNV_SM;			
			double dEpsilon_X_BNV_BNV;
			double dEpsilon_X_SM;
			double dEpsilon_X_BNV;

			double dEpsilon_X_chan;

	};
}

#endif