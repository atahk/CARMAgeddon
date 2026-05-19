#include <RcppParallel.h>
#include <RcppArmadillo.h>

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
				SEXP iObservations, SEXP iVariances, SEXP iVarCenters,
				SEXP iHistTimes,
				SEXP iSizes, SEXP iScreenWidth, SEXP iGrainSize);

RcppExport SEXP monocarptr__fun( SEXP xp, SEXP param_ );
RcppExport SEXP monocarptr__grad( SEXP xp, SEXP param_ );
RcppExport SEXP monocarptr__parallelgrad( SEXP xp, SEXP param_ );
RcppExport SEXP monocarptr__partial( SEXP xp, SEXP param_ );
RcppExport SEXP monocarptr__hess( SEXP xp, SEXP param_ );
RcppExport SEXP monocarptr__hist( SEXP xp, SEXP param_ );
RcppExport SEXP monocarptr__histrand( SEXP xp, SEXP param_, SEXP num_iter_ );
RcppExport SEXP monocarptr__endopt( SEXP success_ );

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
			   SEXP iObservations, SEXP iVariances, SEXP iVarCenters,
			   SEXP iSizes, SEXP iScreenWidth);
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
		      SEXP iSizes, SEXP iScreenWidth);
