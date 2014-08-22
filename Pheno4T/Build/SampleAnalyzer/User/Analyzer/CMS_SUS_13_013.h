#ifndef analysis_CMS_SUS_13_013_h
#define analysis_CMS_SUS_13_013_h

#include "SampleAnalyzer/Process/Analyzer/AnalyzerBase.h"
#include <TGraph.h>
namespace MA5
{
    class CMS_SUS_13_013 : public AnalyzerBase
    {
      INIT_ANALYSIS(CMS_SUS_13_013,"CMS_SUS_13_013")

     public:
        virtual bool Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters);
        virtual void Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files);
        virtual bool Execute(SampleFormat& sample, const EventFormat& event);


        double dFunctionMET(double dGenMET, std::string SearchReg); //function for calculating MET efficiency
        double dFunctionHT(double dGenHT, std::string SearchReg); //function for calculating MET efficiency
        double dFunctionLepton(double dGenLepPt, int lepFlav, std::string SearchReg); //function for calculating Lepton efficiency
        double dFunctionBTag(double dGenJetPt, std::string SearchReg); //function for calculating BTag efficiency
        double dFunctionJetReco(double dGenJetPt, std::string SearchReg); //function for calculating BTag efficiency
        
        bool dFunctionVETO(std::vector leptons1, std::vector leptons2);

        std::vector<double> vecMET;
        std::vector<double> vecMETEff;
        std::vector<double> vecHT;
        std::vector<double> vecHTEff;
        std::vector<double> vecBTag;
        std::vector<double> vecBTagEff;
        std::vector<double> vecJetReco;
        std::vector<double> vecJetRecoEff;
        std::vector<double> vecElectron;
        std::vector<double> vecElectronEff;
        std::vector<double> vecMuon;
        std::vector<double> vecMuonEff;
        
        double dSelectionEff;
        double dLeptonEff;
        double dHTEff;
        double dMETEff;
        double dBTagEff;
        
        double dCounterSelectionEff;
        
     private:
        
        TGraph * EffMET; //Efficiency plot for MET
        TGraph * EffHT; //Efficiency plot for MET
        TGraph * EffJetReco; //Efficiency plot for MET
        TGraph * EffBTag; //Efficiency plot for MET
        TGraph * EffElectron; //Efficiency plot for MET
        TGraph * EffMuon; //Efficiency plot for MET
        
    };
}

#endif