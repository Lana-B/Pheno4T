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


	        double dFunctionBeta(double dBeta, std::string lepton); //function for calculating MET efficiency
			double dFunctionEfficiencies(double dBeta, std::string lepton, std::string TorB, std::string channel);
			double dJetCombiner();
	        int dCounterPassedEvents;
	        double dCounterSelectionEff;
        
		private:
			double dSigma_tt; 
			double dSigma_tW; 

			double dEpsilon_X_SM_SM;
			double dEpsilon_X_BNV_SM;			
			double dEpsilon_X_BNV_BNV;
			double dEpsilon_X_SM;
			double dEpsilon_X_BNV;

			double dEpsilon_X_chan;
			double dE_B_tt;
			double dE_T_tt;
			double dE_B_tW;
			double dE_T_tW; 

			double dN_T_bck;
			double dN_B_bck;
			double dN_B_obs;
			};
}

#endif