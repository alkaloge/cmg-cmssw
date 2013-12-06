///
/// \class l1t::CaloMainProcessor
///
/// Description: interface for the main processor
///
/// Implementation:
///
/// \author: Jim Brooke - University of Bristol
///

//

#ifndef CaloStage1MainProcessor_h
#define CaloStage1MainProcessor_h

//#include "DataFormats/L1TCalorimeter/interface/BXVector.h"
#include <vector>
#include "DataFormats/L1TCalorimeter/interface/EGamma.h"
#include "DataFormats/L1TCalorimeter/interface/Tau.h"
#include "DataFormats/L1TCalorimeter/interface/Jet.h"
#include "DataFormats/L1TCalorimeter/interface/EtSum.h"

#include "FWCore/Framework/interface/Event.h"

namespace l1t {

  class CaloStage1MainProcessor {
  public:
    virtual void processEvent(const EcalTriggerPrimitiveDigiCollection &,
			      const HcalTriggerPrimitiveCollection &,
			      std::vector<EGamma> & egammas,
			      std::vector<Tau> & taus,
			      std::vector<Jet> & jets,
			      std::vector<EtSum> & etsums) = 0;

    virtual ~CaloStage1MainProcessor(){};
  };

}

#endif
