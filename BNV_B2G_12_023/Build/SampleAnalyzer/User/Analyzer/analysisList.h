#include "SampleAnalyzer/User/Analyzer/B2G_12_023.h"
#include "SampleAnalyzer/Process/Analyzer/AnalyzerManager.h"
#include "SampleAnalyzer/Commons/Service/LogStream.h"

// -----------------------------------------------------------------------------
// BuildTable
// -----------------------------------------------------------------------------
void BuildUserTable(MA5::AnalyzerManager& manager)
{
  using namespace MA5;
  manager.Add("B2G_12_023",new B2G_12_023);
}
