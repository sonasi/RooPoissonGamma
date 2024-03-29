/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooHistPoissonGamma.h 42233 2011-11-24 23:35:45Z wouter $
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifdef __MAKECINT__ 
#pragma link C++ myclass+; 
#endif

#ifndef ROO_HIST_POISSONGAMMA
#define ROO_HIST_POISSONGAMMA

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooSetProxy.h"
#include "RooAICRegistry.h"
#include "RooListProxy.h"

#include "Math/Random.h"
#include "Math/GSLRndmEngines.h"

class RooRealVar;
class RooAbsReal;
class RooDataHist;

class RooHistPoissonGamma : public RooAbsPdf {
public:
  RooHistPoissonGamma() ; 
  RooHistPoissonGamma(const char *name, const char *title, const RooArgSet& vars, const RooDataHist& dhist, const RooArgList& histList, const RooArgList& muList, Int_t intOrder=0);

  RooHistPoissonGamma(const RooHistPoissonGamma& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooHistPoissonGamma(*this,newname); }
  virtual ~RooHistPoissonGamma() ;

  RooDataHist& dataHist()  { 
    // Return RooDataHist that is represented
    return *_dataHist ; 
  }
  const RooDataHist& dataHist() const { 
    // Return RooDataHist that is represented
    return *_dataHist ; 
  }
  
  void setInterpolationOrder(Int_t order) { 
    // Set histogram interpolation order 
    _intOrder = order ; 
  }
  Int_t getInterpolationOrder() const { 
    // Return histogram interpolation order
    return _intOrder ; 
  }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

  void setCdfBoundaries(Bool_t flag) { 
    // Set use of special boundary conditions for c.d.f.s
    _cdfBoundaries = flag ; 
  }
  Bool_t getCdfBoundaries() const { 
    // If true, special boundary conditions for c.d.f.s are used
    return _cdfBoundaries ; 
  }

  void setUnitNorm(Bool_t flag) { 
    // Declare contents to have unit normalization
    _unitNorm = flag ; 
  }
  Bool_t haveUnitNorm() const { 
    // Return true if contents is declared to be unit normalized
    return _unitNorm ; 
  }

  virtual Bool_t selfNormalized() const { return _unitNorm ; }

  virtual Int_t getMaxVal(const RooArgSet& vars) const ;
  virtual Double_t maxVal(Int_t code) const ;

  virtual std::list<Double_t>* plotSamplingHint(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const ; 
  virtual std::list<Double_t>* binBoundaries(RooAbsRealLValue& /*obs*/, Double_t /*xlo*/, Double_t /*xhi*/) const ;
  virtual Bool_t isBinnedDistribution(const RooArgSet&) const { return _intOrder==0 ; }

  void genHist(Int_t seed) const ;


protected:

  Bool_t importWorkspaceHook(RooWorkspace& ws) ;
  
  Double_t evaluate() const;
  Double_t totalVolume() const ;
  friend class RooAbsCachedPdf ;
  Double_t totVolume() const ;

  RooListProxy      _muList ;
  RooListProxy      _histList ;
  RooArgSet         _histObsList ; // List of observables defining dimensions of histogram
  RooSetProxy       _pdfObsList ;  // List of observables mapped onto histogram observables
  RooDataHist*      _dataHist ;  // Unowned pointer to underlying histogram
  TIterator*        _muIter ;
  TIterator*         _histIter ; //! 
  TIterator*         _histObsIter ; //! 
  TIterator*         _pdfObsIter ; //! 
  mutable RooAICRegistry _codeReg ; //! Auxiliary class keeping tracking of analytical integration code
  Int_t             _intOrder ; // Interpolation order
  Bool_t            _cdfBoundaries ; // Use boundary conditions for CDFs.
  mutable Double_t  _totVolume ; //! Total volume of space (product of ranges of observables)
  Bool_t            _unitNorm  ; // Assume contents is unit normalized (for use as pdf cache)

  ROOT::Math::Random<ROOT::Math::GSLRngMT>* _gslRan;

  ClassDef(RooHistPoissonGamma,4) // Histogram based PDF
};

#endif
