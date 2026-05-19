#include <cstdlib>
#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>

// #define ARMA_NO_DEBUG
// #define MONOCAR_NO_DEBUG

#define ANALYTIC_GRADIENT
#define ANALYTIC_HESSIAN

#define START_NULL

#ifdef USE_NLOPT
#include <nlopt.h>
#endif

#include "wrapper.h"

using namespace arma;
using namespace std;

typedef arma::mat MyMatrix;
typedef arma::vec MyMatrixCol;
typedef arma::rowvec MyMatrixRow;
typedef arma::cube MyCube;
typedef arma::uword uint;

#define COUT Rcpp::Rcout
#define CHECK_USER_INTERRUPT R_CheckUserInterrupt()
typedef Rcpp::NumericMatrix inputMatrix;
typedef Rcpp::NumericVector inputVector;
typedef Rcpp::IntegerVector inputIntegerVector;
#define INPUTNROW nrow()
#define INPUTNCOL ncol()
#define INPUTINT R_xlen_t

#include "matrixexp.h"


#define LOGTWOPI 1.83787706640934548356

int debugLevel;

int screenWidth;
int screenUse;

MyMatrix vectorToMatrix(const vector<double>& input)
{
  MyMatrix output = conv_to<MyMatrix>::from(input);
  return(output);
}

MyMatrix RtoMatrix(inputMatrix& input)
{
  // MyMatrix output(input.begin(), input.nrow(), input.ncol());
  MyMatrix output=Rcpp::as<MyMatrix>(input);
  return(output);
}

MyMatrix RtoMatrix(inputVector& input)
{
  MyMatrix output(input.begin(), input.size(), 1);
  return(output);
}

MyMatrixCol RtoMatrixCol(inputVector& input)
{
  // MyMatrixCol output(input.begin(), input.size());
  MyMatrixCol output=Rcpp::as<MyMatrixCol>(input);
  return(output);
}

const MyMatrixCol RtoMatrixCol(const inputVector& input)
{
  //MyMatrixCol output(input.begin(), input.size(), false);
  MyMatrixCol output=Rcpp::as<MyMatrixCol>(input);
  return(output);
}

const MyMatrix RtoMatrix(const inputMatrix& input)
{
  //MyMatrix output(input.begin(), input.nrow(), input.ncol(), false);
  MyMatrix output=Rcpp::as<MyMatrix>(input);
  return(output);
}

const MyMatrix RtoMatrix(const inputVector& input)
{
  //MyMatrix output(input.begin(), input.size(), 1, false);
  MyMatrix output=Rcpp::as<MyMatrix>(input);
  return(output);
}

#include "monocar.h"

RcppExport SEXP monocarptr__new(SEXP iParam, SEXP iOutput, SEXP iEpsilon,
				SEXP iStart, SEXP iMean, SEXP iExtraMean,
				SEXP iReversion, SEXP iSimultaneous, SEXP iSigma,
				SEXP iExtraVarAdd, SEXP iExtraVarMult, SEXP iExtraVarPow,
				SEXP iExo,
                                SEXP iPriors,
				SEXP iTimeIndices, SEXP iProcessIndices,
				SEXP iSeriesIndices, SEXP iHouseIndices,
				SEXP iTimeStart, SEXP iTimeEnd,
				SEXP iExoData,
				SEXP iObservations, SEXP iVariances,
				SEXP iVarCenters, SEXP iHistTimes,
				SEXP iSizes, SEXP iScreenWidth, SEXP iGrainSize)
{
#ifdef TRY_X
  char *exceptionMesg=NULL;

  try {
#endif
    
    GetRNGstate();

    inputVector param(iParam);

    inputVector sizes(iSizes);
#ifndef MONOCAR_NO_DEBUG
    assert(sizes.size()>=4);
#endif
    int processIndexSize = sizes(0);
    //int seriesIndexSize = sizes(1);
    int houseIndexSize = sizes(2);
    bool varByHouse = sizes(3) > 0;
    
    inputIntegerVector outputInfo(iOutput);
    double epsilon = Rcpp::as<double>(iEpsilon);
    size_t grainSize = (size_t)(Rcpp::as<int>(iGrainSize));

    bool hasStart = !Rf_isNull(iStart);
    inputMatrix start;
    if (hasStart)
      start = inputMatrix(iStart);
    else
      start = inputMatrix(iMean);
    inputMatrix mean(iMean);
    inputMatrix extraMean(iExtraMean);
    inputMatrix reversion(iReversion);
    inputMatrix simultaneous(iSimultaneous);
    inputMatrix sigma(iSigma);
    inputMatrix extraVarAdd(iExtraVarAdd);
    inputMatrix extraVarMult(iExtraVarMult);
    inputMatrix extraVarPow(iExtraVarPow);
    inputMatrix exo(0,reversion.ncol());
    bool hasExo = !Rf_isNull(iExo);
    if (hasExo)
      exo = inputMatrix(iExo);

    inputMatrix priors(0,0);
    if (!Rf_isNull(iPriors))
      priors = inputMatrix(iPriors);
    
    inputVector timeIndices(iTimeIndices);    
    inputIntegerVector processIndices(iProcessIndices);
    inputIntegerVector seriesIndices(iSeriesIndices);
    inputIntegerVector houseIndices(iHouseIndices);
    inputVector timeStart(iTimeStart);
    inputVector timeEnd(iTimeEnd);
    inputVector observations(iObservations);
    inputVector variances(iVariances);
    inputVector varCenters(iVarCenters);

    MyMatrix exoData(observations.size(), exo.nrow(), arma::fill::zeros);
    if (hasExo)
      {
	if (!Rf_isNull(iExoData))
	  {
	    exoData = Rcpp::as<MyMatrix>(iExoData);
	    assert(exoData.n_cols == uint(exo.nrow()));
	  }
	else
	  hasExo = false;
      }

    vector<double> histTimes = Rcpp::as<vector<double> >(iHistTimes);
    
    MyMatrix paramMat = RtoMatrix(param);

#ifndef MONOCAR_NO_DEBUG
    assert(outputInfo.size() >= 3);
#endif

    if (outputInfo.size() >= 4)
      debugLevel = outputInfo(3) - 1;
    else
      debugLevel = 0;

    screenWidth = Rcpp::as<int>(iScreenWidth);
    screenUse = 0;

    Rcpp::XPtr<OUsuper>
      ptr(new OUsuper(timeIndices,
		      processIndices,
		      processIndexSize,
		      houseIndexSize,
		      seriesIndices,
		      houseIndices,
		      timeStart,
		      timeEnd,
		      exoData,
		      observations,
		      variances,
		      varCenters,
		      hasStart,
		      hasExo,
		      varByHouse,
		      start,
		      mean,
		      extraMean,
		      reversion,
		      simultaneous,
		      sigma,
		      extraVarAdd,
		      extraVarMult,
		      extraVarPow,
		      exo,
                      priors,
		      epsilon,
		      grainSize,
		      histTimes),
	  true );
    
    // return the external pointer to the R side
    return ptr;
    
#ifdef TRY_X
  } catch(std::exception& ex) {
    exceptionMesg = copyMessageToR(ex.what());
  } catch(...) {
    exceptionMesg = copyMessageToR("unknown reason");
  }
    
  if(exceptionMesg != NULL)
    Rf_error(exceptionMesg);
#endif

}

RcppExport SEXP monocarptr__print( SEXP xp, SEXP param_ )
{
#ifdef TRY_X
  char *exceptionMesg=NULL;
  try {
#endif

    Rcpp::XPtr<OUsuper> ou(xp);
    //MyMatrixCol param = Rcpp::as<MyMatrixCol>(param_);
    inputVector param(param_);
    MyMatrix paramMat = RtoMatrix(param);
    
    ou->printParams(paramMat);
    return R_NilValue;
}


RcppExport SEXP monocarptr__fun( SEXP xp, SEXP param_ )
{
#ifdef TRY_X
  char *exceptionMesg=NULL;
  try {
#endif

    Rcpp::XPtr<OUsuper> ou(xp);
    //MyMatrixCol param = Rcpp::as<MyMatrixCol>(param_);
    inputVector param(param_);
    MyMatrix paramMat = RtoMatrix(param);
    
    double negLogLik = - ou->getLikelihood(paramMat);

    if (debugLevel > 4)
      ou->printParams(param);
    
    return Rcpp::wrap(negLogLik);
      
#ifdef TRY_X
  } catch(std::exception& ex) {
    exceptionMesg = copyMessageToR(ex.what());
  } catch(...) {
    exceptionMesg = copyMessageToR("unknown reason");
  }    
  if(exceptionMesg != NULL)
    Rf_error(exceptionMesg);
#endif

}


double numderiv(Rcpp::XPtr<OUsuper> ou, MyMatrix& param, double eps, uint i)
{
  double logLikPlus = -ou->getLikelihood(param, eps, i);
  double logLikMinus = -ou->getLikelihood(param, -eps, i);
  return (logLikPlus - logLikMinus) / (2 * eps);
}

double numderivside(Rcpp::XPtr<OUsuper> ou, MyMatrix& param, double eps, uint i)
{
  double logLikZero = -ou->getLikelihood(param);
  double logLikPlus = -ou->getLikelihood(param, eps, i);
  return (logLikPlus - logLikZero) / (eps);
}

double numsecderiv(Rcpp::XPtr<OUsuper> ou, MyMatrix& param, double eps, uint i, uint j)
{
  double llpp = -ou->getLikelihood(param, eps, i, eps, j);
  double llpm = -ou->getLikelihood(param, eps, i, -eps, j);
  double llmm = -ou->getLikelihood(param, -eps, i, -eps, j);
  double llmp = -ou->getLikelihood(param, -eps, i, eps, j);
  return ((llpp - llpm)/(eps+eps) + (llmm - llmp)/(eps+eps))/(eps+eps);
}

double numsecderivdiag(Rcpp::XPtr<OUsuper> ou, MyMatrix& param, double eps, uint i)
{
  double llo = -ou->getLikelihood(param);
  double llpp = -ou->getLikelihood(param, eps+eps, i);
  double llmm = -ou->getLikelihood(param, -(eps+eps), i);
  return ((llpp - llo)/(eps+eps) + (llmm - llo)/(eps+eps))/(eps+eps);
}


RcppExport SEXP monocarptr__grad( SEXP xp, SEXP param_ )
{
#ifdef TRY_X
  char *exceptionMesg=NULL;
  try {
#endif

    Rcpp::XPtr<OUsuper> ou(xp);
    // MyMatrixCol param = Rcpp::as<MyMatrixCol>(param_);
    inputVector paramX(param_);
    MyMatrix param = RtoMatrix(paramX);
    // double *x = param.memptr();
    double eps = ou->getEpsilon();
    uint n = param.n_rows;
    
    MyMatrixCol gr(n);

    if (debugLevel == 1)
      printDot(":");
    if (debugLevel == 0)
      printDot(".");
    
    if (eps <= 0.0)
      {
	double logLik;
	MyMatrix score;
	ou->getScore(logLik, score, param, 0, n);
#ifndef MONOCAR_NO_DEBUG
	assert(n == score.n_rows);
#endif
	for (uint i = 0; i < n; i++)
	  gr[i] = -score(i, 0);
    }
    else
      {
    for (uint i = 0; i < n; i++)
      {
	double grvalue = (4 * numderiv(ou, param, eps/2, i) - numderiv(ou, param, eps, i))/3;
	if (!isfinite(grvalue))
	  {
	    if (debugLevel > 1)
	      COUT << "nan->" << flush;
	    grvalue = 2 * numderivside(ou, param, eps/2, i) - numderivside(ou, param, eps, i);
	    if (!isfinite(grvalue))
	      {
		if (debugLevel > 1)
		  COUT << "nan->" << flush;
		grvalue = 2 * numderivside(ou, param, -eps/2, i) - numderivside(ou, param, -eps, i);
	      }
	  }
	if (debugLevel > 1)
	  COUT << " " << grvalue << flush;

	gr[i] = grvalue;
      }
      }
    
    return Rcpp::wrap(gr);
      
#ifdef TRY_X
  } catch(std::exception& ex) {
    exceptionMesg = copyMessageToR(ex.what());
  } catch(...) {
    exceptionMesg = copyMessageToR("unknown reason");
  }    
  if(exceptionMesg != NULL)
    Rf_error(exceptionMesg);
#endif

}


struct OUdouble : public RcppParallel::Worker
{
  const Rcpp::XPtr<OUsuper>& mOUsuper;
  const MyMatrix& mParams;
  MyMatrixCol mScore;

  OUdouble(const Rcpp::XPtr<OUsuper>& iOUsuper, const MyMatrix& iParams)
    :mOUsuper(iOUsuper)
    ,mParams(iParams)
    ,mScore(iParams.n_elem, fill::zeros)
  {
  }
    
  OUdouble(const OUdouble& oud, RcppParallel::Split)
    :mOUsuper(oud.mOUsuper)
    ,mParams(oud.mParams)
    ,mScore(oud.mScore.n_elem, fill::zeros)
  {
  }

  void operator()(std::size_t begin, std::size_t end) {
    double logLik;
    MyMatrix score;
    mOUsuper->getScore(logLik, score, mParams, begin, end);
    for (std::size_t k = begin; k < end; k++)
      mScore[k] += score(k-begin,0);
  }

  void join(const OUdouble& rhs) { 
    mScore += rhs.mScore; 
  }
};


RcppExport SEXP monocarptr__parallelgrad( SEXP xp, SEXP param_ )
{
#ifdef TRY_X
  char *exceptionMesg=NULL;
  try {
#endif

    Rcpp::XPtr<OUsuper> ou(xp);
    // MyMatrixCol param = Rcpp::as<MyMatrixCol>(param_);
    inputVector paramX(param_);
    MyMatrix param = RtoMatrix(paramX);
    // double *x = param.memptr();
    double eps = ou->getEpsilon();
    uint n = param.n_rows;
    
    MyMatrixCol gr(n);

    if (debugLevel == 1)
      printDot(":");
    if (debugLevel == 0)
      printDot(".");
    
    if (eps <= 0.0)
      {
	OUdouble pgrad(ou, param);
	RcppParallel::parallelReduce(0, param.n_elem, pgrad, ou->getGrainSize());
#ifndef MONOCAR_NO_DEBUG
	assert(n == pgrad.mScore.n_rows);
#endif
	gr = -pgrad.mScore;
      }
    else
      {
    for (uint i = 0; i < n; i++)
      {
	double grvalue = (4 * numderiv(ou, param, eps/2, i) - numderiv(ou, param, eps, i))/3;
	if (!isfinite(grvalue))
	  {
	    if (debugLevel > 1)
	      COUT << "nan->" << flush;
	    grvalue = 2 * numderivside(ou, param, eps/2, i) - numderivside(ou, param, eps, i);
	    if (!isfinite(grvalue))
	      {
		if (debugLevel > 1)
		  COUT << "nan->" << flush;
		grvalue = 2 * numderivside(ou, param, -eps/2, i) - numderivside(ou, param, -eps, i);
	      }
	  }
	if (debugLevel > 1)
	  COUT << " " << grvalue << flush;

	gr[i] = grvalue;
      }
      }
    
    return Rcpp::wrap(gr);
      
#ifdef TRY_X
  } catch(std::exception& ex) {
    exceptionMesg = copyMessageToR(ex.what());
  } catch(...) {
    exceptionMesg = copyMessageToR("unknown reason");
  }    
  if(exceptionMesg != NULL)
    Rf_error(exceptionMesg);
#endif

}





RcppExport SEXP monocarptr__partial( SEXP xp, SEXP param_ )
{
#ifdef TRY_X
  char *exceptionMesg=NULL;
  try {
#endif

    Rcpp::XPtr<OUsuper> ou(xp);
    // MyMatrixCol param = Rcpp::as<MyMatrixCol>(param_);
    inputVector param(param_);
    MyMatrix paramMat = RtoMatrix(param);
    double outLogLik;
    MyMatrix outPartialScores;
    
    ou->getScore(outLogLik, outPartialScores, paramMat, 0, paramMat.n_rows, true);

    Rcpp::List ret;
    ret["grad"] = outPartialScores.t();
    ret["logLik"] = outLogLik;
    return ret;
      
#ifdef TRY_X
  } catch(std::exception& ex) {
    exceptionMesg = copyMessageToR(ex.what());
  } catch(...) {
    exceptionMesg = copyMessageToR("unknown reason");
  }    
  if(exceptionMesg != NULL)
    Rf_error(exceptionMesg);
#endif

}


RcppExport SEXP monocarptr__hess( SEXP xp, SEXP param_ )
{
#ifdef TRY_X
  char *exceptionMesg=NULL;
  try {
#endif

    Rcpp::XPtr<OUsuper> ou(xp);
    // MyMatrixCol param = Rcpp::as<MyMatrixCol>(param_);
    inputVector paramX(param_);
    MyMatrix param = RtoMatrix(paramX);
    uint n = param.n_rows;
    double eps = ou->getEpsilon();

    MyMatrix hessianMatrix(n,n);

    if (debugLevel > 2)
      ou->printParams(param);

    if (eps <= 0.0)
      {

	if (debugLevel > 1) {
	  COUT << "Computing Hessian" << endl;
	}

	double logLik;
	MyMatrix score;

	ou->getHessian(logLik, score, hessianMatrix, param, 0, n);
#ifndef MONOCAR_NO_DEBUG
	assert(n == score.n_rows);
#endif

	CHECK_USER_INTERRUPT;
	// if (debugLevel > 1)
	//   printDot(".");
	// if (debugLevel == 1 || debugLevel == 0)
	//   printDot("*");
	
      }
    else
      {
#ifndef MONOCAR_NO_DEBUG
      assert(eps != 0);
#endif

    if (debugLevel > 1) {
      COUT << "Computing Hessian" << endl;
      for (uint i = 0; i < n; i++)
    	for (uint j = i; j < n; j++)
    	  COUT << "_";
      COUT << endl;
    }


      for (uint i = 0; i < n; i++)
	{
	  hessianMatrix(i,i) = (4 * numsecderivdiag(ou, param, eps/2, i) - numsecderivdiag(ou, param, eps, i))/3;
	  for (uint j = i+1; j < n; j++)
	    {
	      CHECK_USER_INTERRUPT;
	      hessianMatrix(i,j) = (4 * numsecderiv(ou, param, eps/2, i, j) - numsecderiv(ou, param, eps, i, j))/3;
	      hessianMatrix(j,i) = hessianMatrix(i,j);
	    }
	  if (debugLevel > 1)
	    printDot(".");
	  if (debugLevel == 1 || debugLevel == 0)
	    printDot("*");
	}
      
      if (debugLevel > 1)
	COUT << endl;
      }

    hessianMatrix *= -1.0;

    return Rcpp::wrap(hessianMatrix);
      
#ifdef TRY_X
  } catch(std::exception& ex) {
    exceptionMesg = copyMessageToR(ex.what());
  } catch(...) {
    exceptionMesg = copyMessageToR("unknown reason");
  }    
  if(exceptionMesg != NULL)
    Rf_error(exceptionMesg);
#endif

}


RcppExport SEXP monocarptr__hist( SEXP xp, SEXP param_ )
{
#ifdef TRY_X
  char *exceptionMesg=NULL;
  try {
#endif

    Rcpp::XPtr<OUsuper> ou(xp);
    inputVector param(param_);
    MyMatrix paramMat = RtoMatrix(param);

    vector<MyMatrix> meanMatVector, stdDevMatVector;
    vector<vector<double> > timeVector;
    ou->getHist(paramMat, timeVector, meanMatVector, stdDevMatVector);
    
    vector<Rcpp::NumericMatrix> meanList, stdDevList;
    for (size_t i = 0; i < meanMatVector.size(); ++i)
      {
	Rcpp::NumericMatrix m = Rcpp::wrap(meanMatVector[i]);
	meanList.push_back(m);
      }

    for (size_t i = 0; i < stdDevMatVector.size(); ++i)
      {
	Rcpp::NumericMatrix m = Rcpp::wrap(stdDevMatVector[i]);
	stdDevList.push_back(m);
      }
    
    Rcpp::List ret;
    ret["hist.time"] = timeVector;
    ret["hist.mean"] = meanList;
    ret["hist.sd"] = stdDevList;
    
    return ret;
      
#ifdef TRY_X
  } catch(std::exception& ex) {
    exceptionMesg = copyMessageToR(ex.what());
  } catch(...) {
    exceptionMesg = copyMessageToR("unknown reason");
  }    
  if(exceptionMesg != NULL)
    Rf_error(exceptionMesg);
#endif

}

RcppExport SEXP monocarptr__histrand( SEXP xp, SEXP param_, SEXP num_iter_ )
{
#ifdef TRY_X
  char *exceptionMesg=NULL;
  try {
#endif

    Rcpp::XPtr<OUsuper> ou(xp);
    inputVector paramX(param_);
    MyMatrix param = RtoMatrix(paramX);
    uint n = param.n_rows;
    uint num_iter = Rcpp::as<uint>(num_iter_);

    double logLik;
    MyMatrix score, hessianMatrix, covariance, cholCovariance;

    ou->getHessian(logLik, score, hessianMatrix, param, 0, n);
    if (!inv_sympd(covariance,hessianMatrix))
      covariance = pinv(hessianMatrix);
    cholCovariance = chol(covariance, "lower");    
    
    vector<vector<double> > timeVector;
    vector<MyCube> meanMatVectors, stdDevMatVectors;
    vector<MyMatrix> meanMatVector, stdDevMatVector;

    ou->getHist(param, timeVector, meanMatVector, stdDevMatVector);
    for (size_t i = 0; i < meanMatVector.size(); ++i)
	meanMatVectors.push_back(MyCube(meanMatVector[i].n_rows, meanMatVector[i].n_cols, num_iter));
    for (size_t i = 0; i < stdDevMatVector.size(); ++i)
	stdDevMatVectors.push_back(MyCube(stdDevMatVector[i].n_rows, stdDevMatVector[i].n_cols, num_iter));

    MyMatrix epsilon(param.n_rows, param.n_cols);
    for (size_t iter = 0; iter < num_iter; ++iter)
      {
	for (size_t i = 0; i < epsilon.n_elem; ++i)
	  epsilon[i] = norm_rand();
	ou->getHist(param + cholCovariance * epsilon, timeVector, meanMatVector, stdDevMatVector);
	for (size_t i = 0; i < meanMatVector.size(); ++i)
	  {
	  }
      }
      
    vector<Rcpp::NumericMatrix> meanList, stdDevList;
    for (size_t i = 0; i < meanMatVector.size(); ++i)
      {
	Rcpp::NumericMatrix m = Rcpp::wrap(meanMatVector[i]);
	meanList.push_back(m);
      }

    for (size_t i = 0; i < stdDevMatVector.size(); ++i)
      {
	Rcpp::NumericMatrix m = Rcpp::wrap(stdDevMatVector[i]);
	stdDevList.push_back(m);
      }
    
    Rcpp::List ret;
    ret["hist.time"] = timeVector;
    ret["hist.mean"] = meanList;
    ret["hist.sd"] = stdDevList;
    
    return ret;
      
#ifdef TRY_X
  } catch(std::exception& ex) {
    exceptionMesg = copyMessageToR(ex.what());
  } catch(...) {
    exceptionMesg = copyMessageToR("unknown reason");
  }    
  if(exceptionMesg != NULL)
    Rf_error(exceptionMesg);
#endif

}


RcppExport SEXP monocarptr__endopt( SEXP success_ )
{
  int success = Rcpp::as<int>(success_);
  if (success<=0)
    {
      if ((debugLevel == 1) | (debugLevel == 0))
	printDot("-");
      else if (debugLevel > 1)
	COUT << "nlopt failed!" << endl;
    }
  else
    {
      if ((debugLevel == 1) | (debugLevel == 0))
	printDot("+");
      else if (debugLevel > 1)
	COUT << "found minimum" << endl;
    }
    Rcpp::List ret;
    return ret;
}

#ifdef USE_NLOPT
RcppExport SEXP monocarest(SEXP iParam, SEXP iOutput,
			   SEXP iOptimizers, SEXP iMaxIters,
			   SEXP iTolerances, SEXP iEpsilon,
			   SEXP iLowerBounds, SEXP iUpperBounds,
			   SEXP iStart, SEXP iMean, SEXP iExtraMean,
			   SEXP iReversion, SEXP iSimultaneous, SEXP iSigma,
			   SEXP iExtraVarAdd, SEXP iExtraVarMult, SEXP iExtraVarPow,
			   SEXP iExo,
			   SEXP iTimeIndices, SEXP iProcessIndices,
			   SEXP iSeriesIndices, SEXP iHouseIndices,
			   SEXP iTimeStart, SEXP iTimeEnd,
			   SEXP iExoData,
			   SEXP iObservations, SEXP iVariances,
			   SEXP iVarCenters, SEXP iSizes, SEXP iScreenWidth)
{
#ifdef TRY_X
  char *exceptionMesg=NULL;

  try {
#endif
    
    GetRNGstate();
    nlopt_srand((unsigned long) unif_rand()*numeric_limits<unsigned long>::max());

    inputVector param(iParam);

    inputVector sizes(iSizes);
#ifndef MONOCAR_NO_DEBUG
    assert(sizes.size()>=4);
#endif
    int processIndexSize = sizes(0);
    //int seriesIndexSize = sizes(1);
    int houseIndexSize = sizes(2);
    bool varByHouse = sizes(3) > 0;
    
    inputIntegerVector outputInfo(iOutput);
    vector<string> optimizers = Rcpp::as<vector<string> >(iOptimizers);
    vector<int> maxiters = Rcpp::as<vector<int> >(iMaxIters);
    vector<double> tolerances = Rcpp::as<vector<double> >(iTolerances);
    double epsilon = Rcpp::as<double>(iEpsilon);

    vector<double> lowerBounds(0);
    if (!Rf_isNull(iLowerBounds))
      lowerBounds = Rcpp::as<vector<double> >(iLowerBounds);

    vector<double> upperBounds(0);
    if (!Rf_isNull(iUpperBounds))
      upperBounds = Rcpp::as<vector<double> >(iUpperBounds);

    bool hasStart = !Rf_isNull(iStart);
    inputMatrix start;
    if (hasStart)
      start = inputMatrix(iStart);
    else
      start = inputMatrix(iMean);
    inputMatrix mean(iMean);
    inputMatrix extraMean(iExtraMean);
    inputMatrix reversion(iReversion);
    inputMatrix simultaneous(iSimultaneous);
    inputMatrix sigma(iSigma);
    inputMatrix extraVarAdd(iExtraVarAdd);
    inputMatrix extraVarMult(iExtraVarMult);
    inputMatrix extraVarPow(iExtraVarPow);
    inputMatrix exo(0,reversion.ncol());
    bool hasExo = !Rf_isNull(iExo);
    if (hasExo)
      exo = inputMatrix(iExo);
    
    inputVector timeIndices(iTimeIndices);    
    inputIntegerVector processIndices(iProcessIndices);
    inputIntegerVector seriesIndices(iSeriesIndices);
    inputIntegerVector houseIndices(iHouseIndices);
    inputVector timeStart(iTimeStart);
    inputVector timeEnd(iTimeEnd);
    inputVector observations(iObservations);
    inputVector variances(iVariances);
    inputVector varCenters(iVarCenters);

    MyMatrix exoData;
    if (hasExo)
      {
	if (!Rf_isNull(iExoData))
	  {
	    exoData = Rcpp::as<MyMatrix>(iExoData);
	    assert(exoData.n_cols == uint(exo.nrow()));
	  }
	else
	  hasExo = false;
      }
    
    MyMatrix paramMat = RtoMatrix(param);

#ifndef MONOCAR_NO_DEBUG
    assert(outputInfo.size() >= 3);
#endif

    if (outputInfo.size() >= 4)
      debugLevel = outputInfo(3) - 1;
    else
      debugLevel = 0;

    screenWidth = Rcpp::as<int>(iScreenWidth);
    screenUse = 0;

    //int processIndexSize = maxIndex(processIndices)+1;

    inputVector timeCumStart = getTimeCumStart(timeStart);
    OUsuper ou(timeIndices,
	       processIndices,
	       processIndexSize,
	       houseIndexSize,
	       seriesIndices,
	       houseIndices,
	       timeStart,
	       timeEnd,
	       timeCumStart,
	       exoData,
	       observations,
	       variances,
	       varCenters,
	       hasStart,
	       hasExo,
	       varByHouse,
	       start,
	       mean,
	       extraMean,
	       reversion,
	       simultaneous,
	       sigma,
	       extraVarAdd,
	       extraVarMult,
	       extraVarPow,
	       exo,
	       epsilon);


    double *x = paramMat.memptr();
    double min = 1;

    bool getGrad = false;
    bool getHess = false;
    bool getHist = false;
    if (outputInfo.size() > 0)
      getGrad = outputInfo(0) > 0;
    if (outputInfo.size() > 1)
      getHess = outputInfo(1) > 0;
    if (outputInfo.size() > 2)
      getHist = outputInfo(2) > 0;

    
    for (uint i = 0; i < optimizers.size(); i++)
      {
	bool hasOptimizer = false;
	nlopt_opt opt;
	string optimizerName = optimizers[i];
	if (optimizerName == "DIRECT" ||
	    optimizerName == "DIRECTL" ||
	    optimizerName == "DIRECT_L")
	  opt = nlopt_create(NLOPT_GN_DIRECT_L, paramMat.n_rows);
	if (optimizerName == "CRS" || optimizerName == "CRS2")
	  {
	    opt = nlopt_create(NLOPT_GN_CRS2_LM, paramMat.n_rows);
	    hasOptimizer = true;
	  }
	if (optimizerName == "MLSL")
	  {
	    opt = nlopt_create(NLOPT_G_MLSL_LDS, paramMat.n_rows);
	    hasOptimizer = true;
	  }
	if (optimizerName == "STOGO")
	  {
	    opt = nlopt_create(NLOPT_GD_STOGO, paramMat.n_rows);
	    hasOptimizer = true;
	  }
	if (optimizerName == "STOGORAND" || optimizerName == "STOGO_RAND")
	  {
	    opt = nlopt_create(NLOPT_GD_STOGO_RAND, paramMat.n_rows);
	    hasOptimizer = true;
	  }
	if (optimizerName == "ISRES")
	  {
	    opt = nlopt_create(NLOPT_GN_ISRES, paramMat.n_rows);
	    hasOptimizer = true;
	  }
	if (optimizerName == "ESCH")
	  {
	    opt = nlopt_create(NLOPT_GN_ESCH, paramMat.n_rows);
	    hasOptimizer = true;
	  }
	if (optimizerName == "COBYLA")
	  {
	    opt = nlopt_create(NLOPT_LN_COBYLA, paramMat.n_rows);
	    hasOptimizer = true;
	  }
	if (optimizerName == "BOBYQA")
	  {
	    opt = nlopt_create(NLOPT_LN_BOBYQA, paramMat.n_rows);
	    hasOptimizer = true;
	  }
	if (optimizerName == "NEWUOA")
	  {
	    opt = nlopt_create(NLOPT_LN_NEWUOA, paramMat.n_rows);
	    hasOptimizer = true;
	  }
	if (optimizerName == "PRAXIS")
	  {
	    opt = nlopt_create(NLOPT_LN_PRAXIS, paramMat.n_rows);
	    hasOptimizer = true;
	  }
	if (optimizerName == "SUBPLEX" || optimizerName == "SBPLX")
	  {
	    opt = nlopt_create(NLOPT_LN_SBPLX, paramMat.n_rows);
	    hasOptimizer = true;
	  }
	if (optimizerName == "NM" || optimizerName == "NELDERMEAD")
	  {
	    opt = nlopt_create(NLOPT_LN_NELDERMEAD, paramMat.n_rows);
	    hasOptimizer = true;
	  }
	if (optimizerName == "MMA")
	  {
	    opt = nlopt_create(NLOPT_LD_MMA, paramMat.n_rows);
	    hasOptimizer = true;
	  }
	if (optimizerName == "SLSQP")
	  {
	    opt = nlopt_create(NLOPT_LD_SLSQP, paramMat.n_rows);
	    hasOptimizer = true;
	  }
	if (optimizerName == "BFGS" || optimizerName == "LBFGS")
	  {
	    opt = nlopt_create(NLOPT_LD_LBFGS, paramMat.n_rows);
	    hasOptimizer = true;
	  }
	if (optimizerName == "TNEWTON")
	  {
	    opt = nlopt_create(NLOPT_LD_TNEWTON_RESTART, paramMat.n_rows);
	    hasOptimizer = true;
	  }
	if (optimizerName == "VAR2")
	  {
	    opt = nlopt_create(NLOPT_LD_VAR2, paramMat.n_rows);
	    hasOptimizer = true;
	  }
	if (optimizerName == "VAR1")
	  {
	    opt = nlopt_create(NLOPT_LD_VAR1, paramMat.n_rows);
	    hasOptimizer = true;
	  }
	    
	if (!hasOptimizer)
	  {
	    COUT<< "Skipping unknown optimizer: " << optimizerName << endl;
	  }
	else {
	  nlopt_set_min_objective(opt, oufg, &ou);

	  if (lowerBounds.size() == paramMat.n_rows)
	    {
	      nlopt_set_lower_bounds(opt, lowerBounds.data());
	    }
	  else if (lowerBounds.size() > 0)
	    {
	      COUT << "Invalid lower bounds: bound of size " <<
		lowerBounds.size() << " but should be size " <<
		paramMat.n_rows << endl;
	    }

	  if (upperBounds.size() == paramMat.n_rows)
	    {
	      nlopt_set_upper_bounds(opt, upperBounds.data());
	    }
	  else if (upperBounds.size() > 0)
	    {
	      COUT << "Invalid upper bounds: bound of size " <<
		upperBounds.size() << " but should be size " <<
		paramMat.n_rows << endl;
	    }

	  
	  double tol = 1e-3;
	  if (i < tolerances.size())
	    tol = tolerances[i];
	  nlopt_set_xtol_rel(opt, tol);
	  
	  if (i < maxiters.size())
	    if (maxiters[i] >= 0)
	      nlopt_set_maxeval(opt, maxiters[i]);
	  
	  if (nlopt_optimize(opt, x, &min) < 0)
	    {
	      if ((debugLevel == 1) | (debugLevel == 0))
		printDot("-");
	      else if (debugLevel > 1)
		COUT << "nlopt failed!" << endl;
	    }
	  else
	    {
	      if ((debugLevel == 1) | (debugLevel == 0))
		printDot("+");
	      else if (debugLevel > 1)
		COUT << "found minimum" << endl;
	    }
	  
	  nlopt_destroy(opt);
	  }
      }
     
    Rcpp::List ret;

    ret["min"] = min;
    //inputMatrix outParam=MatrixToR(paramMat);
    ret["par"] = paramMat;

    if (getGrad)
      {
	double outLogLik;
	MyMatrix outPartialScores;

	ou.getScore(outLogLik, outPartialScores, paramMat, 0, paramMat.n_rows, true);

	MyMatrix scoreMat(paramMat.n_rows, 1);
	double *scorePtr = scoreMat.memptr();
	ougradloud(paramMat.n_rows, x, scorePtr, &ou);
	//inputMatrix outScore=MatrixToR(scoreMat);
	ret["score"] = scoreMat;
	
	//inputMatrix outGrad=MatrixToR(outPartialScores.t());
	ret["grad"] = outPartialScores.t();
	ret["logLik"] = outLogLik;
      }
    else
      {
	ret["logLik"] = ou.getLikelihood(paramMat);
      }

    if (getHess)
      {
	MyMatrix hessMat(paramMat.n_rows, paramMat.n_rows);
	double *hess = hessMat.memptr();
	ouhess(paramMat.n_rows, x, hess, &ou);

	//inputMatrix outHess=MatrixToR(hessMat);
	ret["hess"] = hessMat;
      }

    if(getHist)
      {
	MyMatrix timeMat;
	MyMatrix meanMat;
	MyMatrix covMat;
	ou.getHistory(paramMat, timeMat, meanMat, covMat);
	// inputMatrix outTime=MatrixToR(timeMat);
	// inputMatrix outMean=MatrixToR(meanMat);
	// inputMatrix outCov=MatrixToR(covMat);
	ret["hist.time"] = timeMat;
	ret["hist.mean"] = meanMat;
	ret["hist.cov"] = covMat;
      }


    COUT << endl;
    PutRNGstate();
    return ret;
	
#ifdef TRY_X
  } catch(std::exception& ex) {
    exceptionMesg = copyMessageToR(ex.what());
  } catch(...) {
    exceptionMesg = copyMessageToR("unknown reason");
  }
    
  if(exceptionMesg != NULL)
    Rf_error(exceptionMesg);
#endif

}
#endif


RcppExport SEXP ousim(SEXP iParam,
		      SEXP iStart, SEXP iMean, SEXP iExtraMean,
		      SEXP iReversion, SEXP iSimultaneous, SEXP iSigma,
		      SEXP iExtraVarAdd, SEXP iExtraVarMult, SEXP iExtraVarPow,
		      SEXP iExo,
		      SEXP iTimeIndices, SEXP iProcessIndices,
		      SEXP iSeriesIndices, SEXP iHouseIndices,
		      SEXP iTimeStart, SEXP iTimeEnd,
		      SEXP iExoData,
		      SEXP iVariances, SEXP iVarCenters,
		      SEXP iSizes, SEXP iScreenWidth)
{    
    GetRNGstate();
    inputVector param(iParam);
    
    debugLevel = 0;

    screenWidth = Rcpp::as<int>(iScreenWidth);
    screenUse = 0;


    inputVector sizes(iSizes);
#ifndef MONOCAR_NO_DEBUG
    assert(sizes.size()>=4);
#endif
    int processIndexSize = sizes(0);
    //int seriesIndexSize = sizes(1);
    int houseIndexSize = sizes(2);
    bool varByHouse = sizes(3) > 0;

        bool hasStart = !Rf_isNull(iStart);
    inputMatrix start;
    if (hasStart)
      start = inputMatrix(iStart);
    else
      start = inputMatrix(iMean);
    inputMatrix mean(iMean);
    inputMatrix extraMean(iExtraMean);
    inputMatrix reversion(iReversion);
    inputMatrix simultaneous(iSimultaneous);
    inputMatrix sigma(iSigma);
    inputMatrix extraVarAdd(iExtraVarAdd);
    inputMatrix extraVarMult(iExtraVarMult);
    inputMatrix extraVarPow(iExtraVarPow);
    inputMatrix exo(0,reversion.ncol());
    bool hasExo = !Rf_isNull(iExo);
    if (hasExo)
      exo = inputMatrix(iExo);
    
    inputVector timeIndices(iTimeIndices);    
    inputIntegerVector processIndices(iProcessIndices);
    inputIntegerVector seriesIndices(iSeriesIndices);
    inputIntegerVector houseIndices(iHouseIndices);
    inputVector timeStart(iTimeStart);
    inputVector timeEnd(iTimeEnd);
    inputVector variances(iVariances);
    inputVector varCenters(iVarCenters);

    MyMatrix exoData;
    if (hasExo)
      {
	if (!Rf_isNull(iExoData))
	  {
	    exoData = Rcpp::as<MyMatrix>(iExoData);
	    assert(exoData.n_cols == uint(exo.nrow()));
	  }
	else
	  hasExo = false;
      }
    
    MyMatrix paramMat = RtoMatrix(param);

    try{
      OUsimulation ou(paramMat,
		      timeIndices,
		      processIndexSize,
		      houseIndexSize,
		      hasStart,
		      hasExo,
		      varByHouse,
		      varCenters,
		      start,
		      mean,
		      extraMean,
		      reversion,
		      simultaneous,
		      sigma,
		      extraVarAdd,
		      extraVarMult,
		      extraVarPow,
		      exo
		      );

      PutRNGstate();

      if (timeStart.size() > 0)
	return Rcpp::List::create(Rcpp::Named("instantaneous.sample") = ou.getSampleCurrent(),
				  Rcpp::Named("noninstantaneous.sample") = ou.getSampleHist(),
				  Rcpp::Named("obs") = 
				  ou.getSampleData(processIndices,
						   seriesIndices,
						   houseIndices,
						   timeStart,
						   timeEnd,
						   exoData,
						   variances),
				  Rcpp::Named("v") = variances,
				  Rcpp::Named("t1") = timeStart,
				  Rcpp::Named("t2") = timeEnd);
      else
	return Rcpp::List::create(Rcpp::Named("instantaneous.sample") = ou.getSampleCurrent(),
				  Rcpp::Named("noninstantaneous.sample") = ou.getSampleHist());
    }
    catch (const std::runtime_error& error)
      {
	PutRNGstate();
	return Rcpp::List::create();
      }
}
