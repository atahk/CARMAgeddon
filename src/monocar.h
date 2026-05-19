
//#define RTS_SMOOTHER
//#define DOUBLE_CHECK_SMOOTHER


void trimMatrix(MyMatrix& mat, size_t a, size_t b, size_t c, size_t d)
{
    if ((c-a)*(d-b) == 0)
        mat.set_size(c-a, d-b);
    else
        mat = mat.submat(a, b, c-1, d-1);
}

void trimMatrix(MyMatrix& mat, size_t a, size_t b)
{
    return trimMatrix(mat, a, b, mat.n_rows, mat.n_cols);
}


void printDot(const string& input)
{
  screenUse += input.size();
  if ((screenWidth>0) && (screenUse > screenWidth))
    {
      screenUse = static_cast<int>(input.size());
      COUT << endl << input << flush;
    }
  else
    {
      COUT << input << flush;
    }
}

bool hasNegEigen(const MyMatrix& input)
{
  arma::cx_vec eigval = arma::eig_gen(input);
  for (arma::cx_vec::const_iterator ix = eigval.begin();
       ix != eigval.end();
       ix++)
    {
      if (real(*ix)<0)
    return true;
    }
  return false;
}


class MatrixExpFrechet
{
public:
  MatrixExpFrechet(const MyMatrix& iMatrix)
    :baseMatrix(iMatrix)
    ,n(iMatrix.n_rows)
  {
  }

  void push_back(const MyMatrix& iDMatrix)
  {
    bool iDnonzero = false;

    for(MyMatrix::const_iterator i=iDMatrix.begin(); i!=iDMatrix.end(); ++i)
      {
    if (*i != 0.0)
      {
        iDnonzero = true;
        break;
      }
      }

    if (iDnonzero)
      {
    expm.push_back(2*n);
    expm.back().getH().submat(0,0,size(n,n)) = baseMatrix;
    expm.back().getH().submat(n,n,size(n,n)) = baseMatrix;
    expm.back().getH().submat(0,n,size(n,n)) = iDMatrix;
    expm.back().init();
      }
    else
      {
    expm.push_back(0);
      }

    nonzero.push_back(iDnonzero);
  }

  MyMatrix operator()(double scale, uword i) const
  {
    if (nonzero[i])
      return expm[i](scale).submat(0,n,size(n,n));
    else
      return MyMatrix(n,n,fill::zeros);
  }

  void operator()(MyMatrix& iBaseMatrix, MyMatrix& iDMatrix, double scale, uword i) const
  {
    if (nonzero[i])
      {
    iBaseMatrix = expm[i](scale).submat(0,0,size(n,n));
    iDMatrix = expm[i](scale).submat(0,n,size(n,n));
      }
    else
      {
    iBaseMatrix = expm[i](scale).submat(0,0,size(n,n));
    iDMatrix = MyMatrix(n,n,fill::zeros);
      }
  }

private:
  vector<bool> nonzero;
  vector<MatrixExp> expm;
  MyMatrix baseMatrix;
  uword n;
};


class MatrixExpDoubleFrechet
{
public:
  MatrixExpDoubleFrechet(const MyMatrix& iMatrix)
    :baseMatrix(iMatrix)
    ,n(iMatrix.n_rows)
  {
  }

  void push_back(const MyMatrix& iD1Matrix,
         const MyMatrix& iD2Matrix,
         const MyMatrix& iDDMatrix)
  {
    bool iDnonzero = false;

    for(MyMatrix::const_iterator i=iD1Matrix.begin(); i!=iD1Matrix.end(); ++i)
      {
    if (*i != 0.0)
      {
        iDnonzero = true;
        break;
      }
      }

    if (iDnonzero)
      {
    iDnonzero = false;
    for(MyMatrix::const_iterator i=iD2Matrix.begin(); i!=iD2Matrix.end(); ++i)
      {
        if (*i != 0.0)
          {
        iDnonzero = true;
        break;
          }
      }
      }

    if (!iDnonzero)
      {
    for(MyMatrix::const_iterator i=iDDMatrix.begin(); i!=iDDMatrix.end(); ++i)
      {
        if (*i != 0.0)
          {
        iDnonzero = true;
        break;
          }
      }
      }
    
    if (iDnonzero)
      {
    expm.push_back(4*n);
    expm.back().getH().submat(0,0,size(n,n)) = baseMatrix;
    expm.back().getH().submat(n,n,size(n,n)) = baseMatrix;
    expm.back().getH().submat(2*n,2*n,size(n,n)) = baseMatrix;
    expm.back().getH().submat(3*n,3*n,size(n,n)) = baseMatrix;
    expm.back().getH().submat(0,n,size(n,n)) = iD1Matrix;
    expm.back().getH().submat(2*n,3*n,size(n,n)) = iD1Matrix;
    expm.back().getH().submat(0,2*n,size(n,n)) = iD2Matrix;
    expm.back().getH().submat(n,3*n,size(n,n)) = iD2Matrix;
    expm.back().getH().submat(0,3*n,size(n,n)) = iDDMatrix;
    expm.back().init();
      }
    else
      {
    expm.push_back(0);
      }

    nonzero.push_back(iDnonzero);
  }

  MyMatrix operator()(double scale, uword i) const
  {
    if (nonzero[i])
      return expm[i](scale).submat(0,3*n,size(n,n));
    else
      return MyMatrix(n,n,fill::zeros);
  }

  void operator()(MyMatrix& iBaseMatrix, MyMatrix& iD1D2Matrix, double scale, uword i) const
  {
    if (nonzero[i])
      {
    iBaseMatrix = expm[i](scale).submat(0,0,size(n,n));
    iD1D2Matrix = expm[i](scale).submat(0,3*n,size(n,n));
      }
    else
      {
    iBaseMatrix = expm[i](scale).submat(0,0,size(n,n));
    iD1D2Matrix = MyMatrix(n,n,fill::zeros);
      }
  }

private:
  vector<bool> nonzero;
  vector<MatrixExp> expm;
  MyMatrix baseMatrix;
  uword n;
};


MyMatrixCol rmvnorm(const MyMatrixCol mean, const MyMatrix vcov)
{
  MyMatrixCol epsilon(mean.n_rows, arma::fill::randn);
  if (arma::any(arma::vectorise(arma::abs(vcov))>1e-10))
    return mean + arma::chol(vcov) * epsilon;
  else
    return mean;
}

enum passType { timeRemoval = 1, timeAddition = 2, observation = 3, record = 4 };

class seriesPoint
{
public:
  seriesPoint(sword iProcessIndex,
          sword iSeriesIndex,
          sword iHouseIndex,
          double iTimeStart,
          double iTimeEnd,
          double iMean,
          double iVariance,
          MyMatrixRow iExoData)
    :processIndex(iProcessIndex)
    ,seriesIndex(iSeriesIndex)
    ,houseIndex(iHouseIndex)
    ,timeStart(iTimeStart)
    ,timeEnd(iTimeEnd)
    ,mean(iMean)
    ,variance(iVariance)
    ,cumTimeStart(iTimeStart)
    ,exoData(iExoData)
    ,type(observation)
  {
  }

  seriesPoint(sword iProcessIndex,
          double iTime,
          bool recorded)
    :processIndex(iProcessIndex)
    ,timeStart(iTime)
    ,timeEnd(iTime)
    ,mean(0.0)
    ,variance(numeric_limits<double>::infinity())
    ,cumTimeStart(iTime)
    ,exoData()
    ,type(recorded ? record : timeAddition)
  {
  }

  friend bool operator<(const seriesPoint& a, const seriesPoint& b)
  {
    if (a.processIndex != b.processIndex)
      return a.processIndex < b.processIndex;
    else if (a.timeEnd != b.timeEnd)
      return a.timeEnd < b.timeEnd;
    else if (a.timeStart != b.timeStart)
      return a.timeStart < b.timeStart;
    else if (a.type != b.type)
      return a.type < b.type;
    else if (a.seriesIndex != b.seriesIndex)
      return a.seriesIndex < b.seriesIndex;
    else
      return false;
  }

  bool is_record() const
  {
    return type == record;
  }

  bool is_instantaneous() const
  {
    return (type == observation) && (timeStart >= timeEnd);
  }

  bool is_noninstantaneous() const
  {
    return (type == observation) && (timeStart < timeEnd);
  }

  void print(const double ll) const
  {
    if (is_record())
      COUT << "r " << std::flush;
    else if (is_instantaneous())
      COUT << "i " << std::flush;
    else if (is_noninstantaneous())
      COUT << "n " << std::flush;
    else
      COUT << "t " << std::flush;
    COUT << processIndex << " " << seriesIndex << " " << houseIndex << " " << timeStart << " " << timeEnd << " " << cumTimeStart << " " << mean << " " << variance << " " << ll << std::endl;
  }

  double exoProduct(const MyMatrix& param) const
  {
    return as_scalar(exoData * param);
  }
  
  sword processIndex;
  sword seriesIndex;
  sword houseIndex;
  double timeStart;
  double timeEnd;
  double mean;
  double variance;
  double cumTimeStart;
private:
  MyMatrixRow exoData;
  passType type;
};

bool startTimeLessThan(const seriesPoint& a, const seriesPoint& b)
{
  return a.timeStart < b.timeStart;
}

class transformInfo
{
public:
  transformInfo(const MyMatrix& iPinvIexpPntot, const MyMatrix& iExpPnegtime)
    :PinvIexpPntot(iPinvIexpPntot)
    ,expPnegtime(iExpPnegtime)
  {
  }
  
  MyMatrix PinvIexpPntot;
  MyMatrix expPnegtime;
};

class recordInfo
{
public:
  recordInfo(const double iTime, const MyMatrix& iMeanCurrent, const MyMatrix& iMeanHist,
         const MyMatrix& iCovCurrent, const MyMatrix& iCovHist, const MyMatrix& iCovHC)
    :time(iTime)
    ,meanCurrent(iMeanCurrent)
    ,meanHist(iMeanHist)
    ,covCurrent(iCovCurrent)
    ,covHist(iCovHist)
    ,covHC(iCovHC)
  {
  }

  double time;
  MyMatrix meanCurrent, meanHist;
  MyMatrix covCurrent, covHist, covHC;
};

class observationInfo
{
public:
  observationInfo(const bool iInstantaneous,
          const MyMatrix& iHcHC, const MyMatrix& iK, const MyMatrix& iH, const double iy, const double iS)
    :instantaneous(iInstantaneous)
    ,HcHCoS(iHcHC / iS)
    ,Kt(iK.t())
    ,Ht(iH.t())
    ,y(iy)
    ,S(iS)
  {    
  }
  
  bool instantaneous;
  MyMatrix HcHCoS, Kt, Ht;
  double y, S;
};

class modifiedBrysonFrazier
{
public:
  modifiedBrysonFrazier(bool useSmoother)
    :skip(!useSmoother)
  {
  }

    void reset()
    {
    observations.clear();
    transformations.clear();
    records.clear();
    sizes.clear();
    types.clear();
#ifdef RTS_SMOOTHER
      timeBefore.clear();
      timeAfter.clear();
#endif
    }

  void addObservation(const bool instantaneous,
              const MyMatrix& HcHC, const MyMatrix& K, const MyMatrix& H, const double y, const double S)
  {
    if (skip) return;
    observationInfo newInfo(instantaneous, HcHC, K, H, y, S);
    observations.push_back(newInfo);
    types.push_back(observation);
  }

  void addTime(const MyMatrix& PinvIexpPntot, const MyMatrix& expPnegtime)
  {
    if (skip) return;
    transformInfo newInfo(PinvIexpPntot, expPnegtime);
    transformations.push_back(newInfo);
    types.push_back(timeAddition);
  }
  
  void removeTime(const size_t s)
  {
    if (skip) return;
    sizes.push_back(s);
    types.push_back(timeRemoval);
  }

  void addRecord(const double time,
         const MyMatrix& meanCurrent, const MyMatrix& meanHist,
         const MyMatrix& covCurrent, const MyMatrix& covHist, const MyMatrix& covHC)
  {
    if (skip) return;
    recordInfo newInfo(time,
               meanCurrent, meanHist,
               covCurrent, covHist, covHC);
    records.push_back(newInfo);
    types.push_back(record);
  }

#ifdef RTS_SMOOTHER
    void addPassData(const bool after, const double time,
                     const MyMatrix& meanCurrent, const MyMatrix& meanHist,
                     const MyMatrix& covCurrent, const MyMatrix& covHist, const MyMatrix& covHC)
    {
      recordInfo newInfo(time,
                         meanCurrent, meanHist,
                         covCurrent, covHist, covHC);
      if (after)
          timeAfter.push_back(newInfo);
      else
          timeBefore.push_back(newInfo);
    }

    void terminate(const MyMatrix& meanCurrent, const MyMatrix& meanHist,
                   const MyMatrix& covCurrent, const MyMatrix& covHist, const MyMatrix& covHC)

    {
      kappaLower = join_cols(meanCurrent, meanHist);
      kappaUpper = join_cols(join_rows(covCurrent, covHC.t()),
                             join_rows(covHC, covHist));
    }
#endif

  void smoothBackward(vector<double>& times, MyCube& means, MyCube& vcovs)
  {
    if (records.size() == 0)
      return;

    means.set_size(records.back().meanCurrent.n_rows,
           records.back().meanCurrent.n_cols,
           records.size());
    vcovs.set_size(records.back().covCurrent.n_rows,
           records.back().covCurrent.n_cols,
           records.size());
    times.resize(records.size());
    
    MyMatrix lambdaLowerCurrent(arma::size(records.back().meanCurrent), arma::fill::zeros);
    MyMatrix lambdaLowerHist(arma::size(records.back().meanHist), arma::fill::zeros);
    MyMatrix lambdaUpperCurrent(arma::size(records.back().covCurrent), arma::fill::zeros);
    MyMatrix lambdaUpperHist(arma::size(records.back().covHist), arma::fill::zeros);
    MyMatrix lambdaUpperHC(arma::size(records.back().covHC), arma::fill::zeros);
//    MyMatrix lambdaUpperCH(arma::size(records.back().covHC.t()), fill::zeros);

    for (; types.size() > 0; types.pop_back())
      {
        if (types.back() == observation)
          {
            observationInfo& observation = observations.back();

            if (observation.S <= 0.0)
              {
              }
            else if (observation.instantaneous)
              {
                lambdaLowerCurrent -= observation.Ht * observation.Kt * lambdaLowerCurrent;
                if (lambdaLowerHist.n_rows > 0)
                    lambdaLowerCurrent -= observation.Ht * observation.HcHCoS * lambdaLowerHist;
                lambdaLowerCurrent -= observation.Ht * observation.y / observation.S;

//                lambdaUpperCurrent -= (lambdaUpperCurrent * observation.Kt.t() + lambdaUpperHC.t() * observation.HcHCoS.t()) * observation.Ht.t();
                lambdaUpperCurrent -= lambdaUpperCurrent * observation.Kt.t() * observation.Ht.t();
                if (lambdaLowerHist.n_rows > 0)
                  {
                    lambdaUpperCurrent -= lambdaUpperHC.t() * observation.HcHCoS.t() * observation.Ht.t();
                    lambdaUpperHC -= (lambdaUpperHC * observation.Kt.t() + lambdaUpperHist * observation.HcHCoS.t()) * observation.Ht.t();
                  }

//                lambdaUpperCurrent -= observation.Ht * (observation.Kt * lambdaUpperCurrent + observation.HcHCoS * lambdaUpperHC);
                lambdaUpperCurrent -= observation.Ht * observation.Kt * lambdaUpperCurrent;
                if (lambdaLowerHist.n_rows > 0)
                    lambdaUpperCurrent -= observation.Ht * observation.HcHCoS * lambdaUpperHC;

                lambdaUpperCurrent += observation.Ht * observation.Ht.t() / observation.S;
              }
            else
              {
                lambdaLowerHist -= observation.Ht * observation.Kt * lambdaLowerHist;
                lambdaLowerHist -= observation.Ht * observation.HcHCoS * lambdaLowerCurrent;
                lambdaLowerHist -= observation.Ht * (observation.y / observation.S);

                lambdaUpperHist -= observation.Ht * (observation.Kt * lambdaUpperHist + observation.HcHCoS * lambdaUpperHC.t());
                lambdaUpperHC -= observation.Ht * (observation.Kt * lambdaUpperHC + observation.HcHCoS * lambdaUpperCurrent);

                lambdaUpperHist -= (lambdaUpperHist * observation.Kt.t() + lambdaUpperHC * observation.HcHCoS.t()) * observation.Ht.t();

                lambdaUpperHist += observation.Ht * observation.Ht.t() / observation.S;
              }

#ifdef RTS_SMOOTHER
#ifdef DOUBLE_CHECK_SMOOTHER
            COUT << "X" << (observation.instantaneous?"I":"N") << lambdaLowerHist.n_rows << std::endl;

            MyMatrix kappaMean = kappaLower.head_rows(lambdaLowerCurrent.n_rows);
            MyMatrix lambdaMean = timeBefore.back().meanCurrent - timeBefore.back().covCurrent * lambdaLowerCurrent  - timeBefore.back().covHC.t() * lambdaLowerHist;
            MyMatrix error = lambdaMean - kappaMean;

            if (accu(abs(error)) > 1e-4)
              {
                COUT << kappaMean.t() << lambdaMean.t() << std::flush;
                COUT << error.t() << std::flush;
                COUT << std::endl;
              }

            MyMatrix kappaCov = kappaUpper;

            MyMatrix P = join_cols(join_rows(timeBefore.back().covCurrent, timeBefore.back().covHC.t()),
                                   join_rows(timeBefore.back().covHC, timeBefore.back().covHist));
            MyMatrix Lhat = join_cols(join_rows(lambdaUpperCurrent, lambdaUpperHC.t()),
                                      join_rows(lambdaUpperHC, lambdaUpperHist));

            MyMatrix lambdaCov = P - P * Lhat * P;

            MyMatrix errorCov = lambdaCov - kappaCov;

            if ((errorCov.n_elem > 0) && (accu(abs(errorCov)) > 1e-4))
              {
                COUT << P.i() - P.i() * kappaUpper * P.i() << std::endl;
                COUT << Lhat << std::endl;

                COUT << kappaCov << std::endl << lambdaCov << std::endl;
                COUT << errorCov << std::flush;
                COUT << std::endl;
              }

            timeBefore.pop_back();
#endif
#endif

            observations.pop_back();
          }
        else if (types.back() == timeAddition)
          {
            transformInfo& transformation = transformations.back();

            lambdaLowerCurrent = transformation.expPnegtime.t() * lambdaLowerCurrent;
            lambdaLowerCurrent += transformation.PinvIexpPntot.t() * lambdaLowerHist.tail_rows(transformation.PinvIexpPntot.n_cols);
            lambdaLowerHist.resize(lambdaLowerHist.n_rows - transformation.PinvIexpPntot.n_cols, lambdaLowerHist.n_cols);

            lambdaUpperCurrent = transformation.expPnegtime.t() * lambdaUpperCurrent * transformation.expPnegtime;
            if (transformation.PinvIexpPntot.n_rows > 0)
              {
                lambdaUpperCurrent += transformation.PinvIexpPntot.t() * lambdaUpperHC.tail_rows(transformation.PinvIexpPntot.n_rows) * transformation.expPnegtime;
                lambdaUpperCurrent += transformation.expPnegtime.t() * lambdaUpperHC.tail_rows(transformation.PinvIexpPntot.n_rows).t() * transformation.PinvIexpPntot;
                lambdaUpperCurrent += transformation.PinvIexpPntot.t() * lambdaUpperHist.submat(lambdaUpperHist.n_rows-transformation.PinvIexpPntot.n_rows, lambdaUpperHist.n_cols-transformation.PinvIexpPntot.n_rows, lambdaUpperHist.n_rows-1, lambdaUpperHist.n_cols-1) * transformation.PinvIexpPntot;
              }

            lambdaUpperHC = lambdaUpperHC.head_rows(lambdaLowerHist.n_rows) * transformation.expPnegtime;
            if (lambdaLowerHist.n_rows > 0)
                lambdaUpperHC += lambdaUpperHist.submat(0, lambdaLowerHist.n_rows, lambdaUpperHC.n_rows - 1, lambdaUpperHist.n_rows - 1) * transformation.PinvIexpPntot;

            lambdaUpperHist.resize(lambdaLowerHist.n_rows, lambdaLowerHist.n_rows);


#ifdef RTS_SMOOTHER
            MyMatrix meanBefore = join_cols(timeBefore.back().meanCurrent, timeBefore.back().meanHist);
            MyMatrix meanAfter = join_cols(timeAfter.back().meanCurrent, timeAfter.back().meanHist);

            MyMatrix covBefore =
            join_cols(join_rows(timeBefore.back().covCurrent, timeBefore.back().covHC.t()),
                      join_rows(timeBefore.back().covHC, timeBefore.back().covHist));
            MyMatrix covAfter =
            join_cols(join_rows(timeAfter.back().covCurrent, timeAfter.back().covHC.t()),
                      join_rows(timeAfter.back().covHC, timeAfter.back().covHist));

            MyMatrix F(meanAfter.n_rows, meanBefore.n_rows, arma::fill::eye);

            F.submat(0, 0, transformation.expPnegtime.n_rows-1, transformation.expPnegtime.n_cols-1) = transformation.expPnegtime;

            F.submat(meanBefore.n_rows, 0, meanAfter.n_rows-1, transformation.PinvIexpPntot.n_cols-1) = transformation.PinvIexpPntot;

            MyMatrix C = covBefore * solve(covAfter, F).t();

            kappaLower = meanBefore + C * (kappaLower - meanAfter);

            kappaUpper = covBefore + C * (kappaUpper - covAfter) * C.t();

#ifdef DOUBLE_CHECK_SMOOTHER

            COUT << "+" << std::endl;

            MyMatrix kappaMean = kappaLower.head_rows(lambdaLowerCurrent.n_rows);
            MyMatrix lambdaMean = timeBefore.back().meanCurrent - timeBefore.back().covCurrent * lambdaLowerCurrent  - timeBefore.back().covHC.t() * lambdaLowerHist;
            MyMatrix error = lambdaMean - kappaMean;

            if (accu(abs(error)) > 1e-4)
              {
                COUT << kappaMean.t() << lambdaMean.t() << std::flush;
                COUT << error.t() << std::flush;
                COUT << std::endl;
              }

            MyMatrix kappaCov = kappaUpper;

            MyMatrix P = join_cols(join_rows(timeBefore.back().covCurrent, timeBefore.back().covHC.t()),
                                   join_rows(timeBefore.back().covHC, timeBefore.back().covHist));
            MyMatrix Lhat = join_cols(join_rows(lambdaUpperCurrent, lambdaUpperHC.t()),
                                      join_rows(lambdaUpperHC, lambdaUpperHist));

            MyMatrix lambdaCov = P - P * Lhat * P;

            MyMatrix errorCov = lambdaCov - kappaCov;

            if (accu(abs(errorCov)) > 1e-4)
              {
                COUT << P.i() - P.i() * kappaUpper * P.i() << std::endl;
                COUT << Lhat << std::endl;

                COUT << kappaCov << std::endl << lambdaCov << std::endl;
                COUT << errorCov << std::flush;
                COUT << std::endl;
              }
#endif

            timeBefore.pop_back();
            timeAfter.pop_back();
#endif

            transformations.pop_back();
          }
        else if (types.back() == timeRemoval)
          {
            size_t s = sizes.back();

            MyMatrix newLowerHist(lambdaLowerHist.n_rows + s, lambdaLowerHist.n_cols, fill::zeros);
            if (newLowerHist.n_rows > s)
                newLowerHist.submat(s, 0, newLowerHist.n_rows - 1, newLowerHist.n_cols - 1) = lambdaLowerHist;
            lambdaLowerHist = newLowerHist;

            MyMatrix newUpperHist(lambdaUpperHist.n_rows + s, lambdaUpperHist.n_cols + s, fill::zeros);
            if (newUpperHist.n_rows > s)
                newUpperHist.submat(s, s, newUpperHist.n_rows - 1, newUpperHist.n_cols - 1) = lambdaUpperHist;
            lambdaUpperHist = newUpperHist;

            MyMatrix newUpperHC(lambdaUpperHC.n_rows + s, lambdaUpperHC.n_cols, fill::zeros);
            if (newUpperHC.n_rows > s)
                newUpperHC.submat(s, 0, newUpperHC.n_rows - 1, newUpperHC.n_cols - 1) = lambdaUpperHC;
            lambdaUpperHC = newUpperHC;

#ifdef RTS_SMOOTHER
            MyMatrix meanBefore = join_cols(timeBefore.back().meanCurrent, timeBefore.back().meanHist);
            MyMatrix meanAfter = join_cols(timeAfter.back().meanCurrent, timeAfter.back().meanHist);

            MyMatrix covBefore =
            join_cols(join_rows(timeBefore.back().covCurrent, timeBefore.back().covHC.t()),
                      join_rows(timeBefore.back().covHC, timeBefore.back().covHist));
            MyMatrix covAfter =
            join_cols(join_rows(timeAfter.back().covCurrent, timeAfter.back().covHC.t()),
                      join_rows(timeAfter.back().covHC, timeAfter.back().covHist));

            MyMatrix F(meanAfter.n_rows, meanBefore.n_rows, arma::fill::zeros);
            F.submat(0, 0, timeAfter.back().meanCurrent.n_rows-1, timeAfter.back().meanCurrent.n_rows-1).eye();
            if (timeAfter.back().meanHist.n_rows > 0)
                F.submat(F.n_rows-timeAfter.back().meanHist.n_rows, F.n_cols-timeAfter.back().meanHist.n_rows, F.n_rows-1, F.n_cols-1).eye();

            MyMatrix C = covBefore * F.t() * covAfter.i();

            kappaLower = meanBefore + C * (kappaLower - meanAfter);

            kappaUpper = covBefore + C * (kappaUpper - covAfter) * C.t();

#ifdef DOUBLE_CHECK_SMOOTHER
            MyMatrix kappaMean = kappaLower.head_rows(lambdaLowerCurrent.n_rows);
            MyMatrix lambdaMean = timeBefore.back().meanCurrent - timeBefore.back().covCurrent * lambdaLowerCurrent  - timeBefore.back().covHC.t() * lambdaLowerHist;
            MyMatrix error = lambdaMean - kappaMean;

            if (accu(abs(error)) > 1e-4)
              {
                COUT << kappaMean.t() << lambdaMean.t() << std::flush;
                COUT << error.t() << std::flush;
                COUT << std::endl;
              }

            COUT << "-" << std::endl;

            MyMatrix kappaCov = kappaUpper;

            MyMatrix P = join_cols(join_rows(timeBefore.back().covCurrent, timeBefore.back().covHC.t()),
                                   join_rows(timeBefore.back().covHC, timeBefore.back().covHist));
            MyMatrix Lhat = join_cols(join_rows(lambdaUpperCurrent, lambdaUpperHC.t()),
                                      join_rows(lambdaUpperHC, lambdaUpperHist));

            MyMatrix lambdaCov = P - P * Lhat * P;

            MyMatrix errorCov = lambdaCov - kappaCov;

            if (accu(abs(errorCov)) > 1e-4)
              {
                COUT << P.i() - P.i() * kappaUpper * P.i() << std::endl;
                COUT << Lhat << std::endl;

                COUT << kappaCov << std::endl << lambdaCov << std::endl;
                COUT << errorCov << std::flush;
                COUT << std::endl;
              }
#endif

            timeBefore.pop_back();
            timeAfter.pop_back();
#endif

            sizes.pop_back();
          }
        else if (types.back() == record)
          {
            recordInfo& record = records.back();

            means.slice(records.size()-1) = record.meanCurrent - record.covCurrent * lambdaLowerCurrent  - record.covHC.t() * lambdaLowerHist;
            vcovs.slice(records.size()-1) = record.covCurrent - record.covCurrent * lambdaUpperCurrent * record.covCurrent - record.covHC.t() * lambdaUpperHist * record.covHC - record.covCurrent * lambdaUpperHC.t() * record.covHC - record.covHC.t() * lambdaUpperHC * record.covCurrent;

#ifdef RTS_SMOOTHER
#ifdef DOUBLE_CHECK_SMOOTHER
            if (accu(abs(means.slice(records.size()-1) - kappaLower.head_rows(record.meanCurrent.n_rows))) > 1e-4)
              {
                COUT << "X" << std::endl;
              }
            if (accu(abs(vcovs.slice(records.size()-1) - kappaUpper.submat(0,0,record.covCurrent.n_rows-1, record.covCurrent.n_cols-1))) > 1e-4)
              {
                MyMatrix P = join_cols(join_rows(record.covCurrent, record.covHC.t()),
                                       join_rows(record.covHC, record.covHist));
                MyMatrix Lhat = join_cols(join_rows(lambdaUpperCurrent, lambdaUpperHC.t()),
                                          join_rows(lambdaUpperHC, lambdaUpperHist));
                MyMatrix V2 = P - P * Lhat * P;

                COUT << V2;
                COUT << std::endl;
                COUT << kappaUpper;
                COUT << std::endl;
                COUT << "X" << std::endl;
              }
#else
            means.slice(records.size()-1) = kappaLower.submat(0, 0, means.n_rows-1, means.n_cols-1);
            vcovs.slice(records.size()-1) = kappaUpper.submat(0, 0, vcovs.n_rows-1, vcovs.n_cols-1);
#endif
#endif

            times[records.size()-1] = record.time;

            records.pop_back();
          }
      }
  }
    
private:
  bool skip;
  vector<observationInfo> observations;
  vector<transformInfo> transformations;
  vector<recordInfo> records;
  vector<size_t> sizes;
  vector<passType> types;
#ifdef RTS_SMOOTHER
    vector<recordInfo> timeBefore;
    vector<recordInfo> timeAfter;
    MyMatrix kappaLower;
    MyMatrix kappaUpper;
#endif
};

/*
double log_mult_gamma(const double p, const double a)
{
  double sum = p * (p-1) * 0.25;
  for (uword j = 1; j <= p; ++j)
    {
      sum += std::lgamma(a + 0.5 * (1.0 - double(j)));
    }
  return sum;
};
*/


 /*
double log_dwish(const MyMatrix& X, const double n_minus_p, const MyMatrix& V)
{
  double log_prob = 0.0;
  double p = double(X.n_rows);
  
  double log_det_X, sign;
  arma::log_det(log_det_X, sign, X);
  
  log_prob += log_det_X * 0.5 * (n_minus_p - 1.0);
  log_prob -= 0.5 * arma::trace(arma::solve(V, X));
  //log_prob -= 0.5 * (n_minus_p + p) * p * M_LN2;
  //double log_det_V;
  //arma::log_det(log_det_V, sign, V);
  //log_prob -= 0.5 * (n_minus_p + p) * log_det_V;
  //log_prob -= log_mult_gamma(p, 0.5 * (n_minus_p + p));
  
  return log_prob;
};

double D_log_dwish(const MyMatrix& X, const MyMatrix& D_X, const double n_minus_p, const MyMatrix& V)
{
  double log_prob = 0.0;
  double p = double(X.n_rows);

  double D_log_det_X = arma::trace(arma::solve(X, D_X));
  
  log_prob += D_log_det_X * 0.5 * (n_minus_p - 1.0);
  log_prob -= 0.5 * trace(solve(V, D_X));
  
  return log_prob;
};
 */


/* V=v*eye(p) */
double log_dinvwish(const MyMatrix& X, const double n_minus_p, const double v)
{
  double log_prob = 0.0;
  double p = double(X.n_rows);
  double df = n_minus_p + p;

  double trace_X_inv = 0.0;
  double log_det_X = 0.0;
  arma::vec eigval = eig_sym(X);
  for (arma::vec::const_iterator i = eigval.begin(); i != eigval.end(); ++i)
    {
      if (*i == 0.0)
        return -std::numeric_limits<double>::infinity();
      trace_X_inv += 1.0 / (*i);
      log_det_X += std::log(*i);
    }
  
  log_prob -= 0.5 * trace_X_inv * v;
  log_prob -= log_det_X * (df + p + 1) * 0.5;

  // log_prob += p * std::log(v) * df * 0.5;
  // log_prob -= 0.5 * df * p * M_LN2;
  // log_prob -= log_mult_gamma(p, 0.5 * df);
  
  return log_prob;
};

double D_log_dinvwish(const MyMatrix& X, const MyMatrix& D_X, const double n_minus_p, const double v)
{
  if (!arma::any(arma::vectorise(D_X)))
    return 0.0;
  
  double log_prob = 0.0;
  double p = double(X.n_rows);
  double df = n_minus_p + p;
  MyMatrix X_inv;  
  if (!arma::inv_sympd(X_inv, X))
    return 0.0;
    
  double D_log_det_X = arma::trace(X_inv * D_X);
  double D_trace_X_inv = arma::trace(X_inv * D_X * X_inv);
  
  log_prob += 0.5 * D_trace_X_inv * v;
  log_prob -= 0.5 * (df + p + 1) * D_log_det_X;

  // if (!arma::solve(X_inv_D_X, X / det_X, D_X, solve_opts::no_approx))
  //   return 0.0;
  
  // double D_log_det_X = arma::trace(X_inv_D_X);
  
  // log_prob += 0.5 * D_log_det_X * v;
  // log_prob -= 0.5 * (df + p + 1) * det_X * D_log_det_X;
  
  return log_prob;
};

double DD_log_dinvwish(const MyMatrix& X, const MyMatrix& D_X, const MyMatrix& DD_X, const double n_minus_p, const double v)
{
  if (!arma::any(arma::vectorise(D_X)))
    return 0.0;
  
  double log_prob = 0.0;
  double p = double(X.n_rows);
  double df = n_minus_p + p;
  MyMatrix X_inv;  
  if (!arma::inv_sympd(X_inv, X))
    return 0.0;

  MyMatrix D_X_inv = - X_inv * D_X * X_inv;
    
  double DD_log_det_X = arma::trace(D_X_inv * D_X) + arma::trace(X_inv * DD_X);
  double DD_trace_X_inv = arma::trace(D_X_inv * D_X * X_inv) +
    arma::trace(X_inv * DD_X * X_inv) +
    arma::trace(X_inv * D_X * D_X_inv);
  
  log_prob += 0.5 * DD_trace_X_inv * v;
  log_prob -= 0.5 * (df + p + 1) * DD_log_det_X;

  // if (!arma::solve(X_inv_D_X, X / det_X, D_X, solve_opts::no_approx))
  //   return 0.0;
  
  // double D_log_det_X = arma::trace(X_inv_D_X);
  
  // log_prob += 0.5 * D_log_det_X * v;
  // log_prob -= 0.5 * (df + p + 1) * det_X * D_log_det_X;
  
  return log_prob;
};


double log_dgamma(const double x, const double a, const double b)
{
  double log_prob = 0.0;

  /* log_prob += a * std::log(b); */
  /* log_prob -= std::lgamma(a); */
  log_prob -= b * x;
  if (a != 1.0)
    {
      if (x == 0.0)
	{
	  if (a > 1.0)
	    log_prob = -numeric_limits<double>::infinity();
	  else
	    log_prob = numeric_limits<double>::infinity();
	}
      else
	log_prob += (a - 1.0) * std::log(x);
    }
  
  return log_prob;
};

double D_log_dgamma(const double x, const double D_x, const double a, const double b)
{
  double log_prob = 0.0;

  log_prob -= b * D_x;
  if (a != 1.0)
    {
      if (x == 0.0)
	return 0.0;
      else
	log_prob += (a - 1.0) * D_x / x;
    }
  
  return log_prob;
};

double DD_log_dgamma(const double x, const double D_x, const double DD_x, const double a, const double b)
{
  double log_prob = 0.0;

  log_prob -= b * DD_x;
  if (a != 1.0)
    {
      if (x == 0.0)
	return 0.0;
      else
	log_prob += (a - 1.0) * (DD_x / x - D_x / (x*x));
    }
  
  return log_prob;
};

double log_dexp(const double x, const double b)
{
  if (x < 0.0)
    return -std::numeric_limits<double>::infinity();
  
  double log_prob = 0.0;

  /* log_prob += std::log(b); */
  log_prob -= b * x;
  
  return log_prob;
};

double D_log_dexp(const double x, const double D_x, const double b)
{
  if (x < 0.0)
    return 0.0;
  
  double log_prob = 0.0;

  log_prob -= b * D_x;
  
  return log_prob;
};

double DD_log_dexp(const double x, const double D_x, const double DD_x, const double b)
{
  if (x < 0.0)
    return 0.0;
  
  double log_prob = 0.0;

  log_prob -= b * D_x;
  
  return log_prob;
};

class OUparam
{
public:
  bool hasStart;
  MyMatrix start;
  MyMatrix mean;
  MyMatrix extraMean;
  MyMatrix reversion;
  MyMatrix sigma;
  MyMatrix extraVarAdd;
  MyMatrix extraVarMult;
  MyMatrix extraVarPow;
  MyMatrix exo;
  MyMatrixCol varCenters;
  bool hasExo;
  bool varByHouse;
    bool muByHouse;
    bool deltaByHouse;

    OUparam process(const size_t processIndex) const
    {
      OUparam out;
      out.hasStart = hasStart;
      if (processIndex < start.n_cols)
          out.start = start.col(processIndex);
      else
          out.start = start;
      if (processIndex < mean.n_cols)
          out.mean = mean.col(processIndex);
      else
          out.mean = mean;
      if (processIndex < extraMean.n_cols)
          out.extraMean = extraMean.col(processIndex);
      else
          out.extraMean = extraMean;
      out.reversion = reversion;
      out.sigma = sigma;
      out.extraVarAdd = extraVarAdd;
      out.extraVarMult = extraVarMult;
      out.extraVarPow = extraVarPow;
      out.exo = exo;
      out.varCenters = varCenters;
      out.hasExo = hasExo;
      out.varByHouse = varByHouse;
      out.muByHouse = muByHouse;
      out.deltaByHouse = deltaByHouse;
      return out;
    }

  double getExtraVarAdd(const size_t houseIndex, const size_t seriesIndex) const
  {
    if (varByHouse)
      return extraVarAdd(houseIndex, 0);
    else
      return extraVarAdd(seriesIndex, 0);
  }
  double getExtraVarMult(const size_t houseIndex, const size_t seriesIndex) const
  {
    if (varByHouse)
      return extraVarMult(houseIndex, 0);
    else
      return extraVarMult(seriesIndex, 0);
  }
  double getExtraVarPow(const size_t houseIndex, const size_t seriesIndex) const
  {
    if (varByHouse)
      return extraVarPow(houseIndex, 0);
    else
      return extraVarPow(seriesIndex, 0);
  }
  double getVarCenter(const size_t houseIndex, const size_t seriesIndex) const
  {
    if (varByHouse)
      return varCenters(houseIndex);
    else
      return varCenters(seriesIndex);
  }
};

class OUprior
{
 public:
 OUprior()
   :sigma_prior_df_plus(1.0)
    ,sigma_prior_scale(1.0)
    ,gamma_shape(2.0)
    {
    }

  static std::vector<OUprior> getPriorFromMatrix(const inputMatrix& x)
  {
    std::vector<OUprior> priors;
    for (INPUTINT i = 0; i < x.INPUTNROW; ++i)
      {
	OUprior prior;
	priors.push_back(prior);
      }
    return priors;
  }
  
  double getLogProb(const OUparam& param) const
    {
      double log_prob = 0.0;
      
      /* mean */
      
      /* extraMean */
      
      /* reversion */
      /*arma::cx_vec eigval = arma::eig_gen(param.reversion);
      for (arma::cx_vec::const_iterator ix = eigval.begin();
      	   ix != eigval.end();
      	   ix++)
      	{
      	  log_prob += std::log(std::real(*ix));
          log_prob -= std::real(*ix);
          }*/
      /* double val, sign; */
      /* arma::log_det(val, sign, param.reversion); */
      /* if (sign > 0.0) */
      /*   log_prob += val; */
      /* else */
      /*   return std::numeric_limits<double>::infinity(); */
      log_prob -= arma::trace(param.reversion);
      
      /* sigma */
      log_prob += log_dinvwish(param.sigma, sigma_prior_df_plus, sigma_prior_scale / (sigma_prior_df_plus + param.sigma.n_rows));
      
      for (MyMatrix::const_iterator i = param.extraVarAdd.begin();
	   i != param.extraVarAdd.end();
	   ++i)
	log_prob += log_dexp(*i, 1.0);
      for (MyMatrix::const_iterator i = param.extraVarMult.begin();
	   i != param.extraVarMult.end();
	   ++i)
	log_prob += log_dgamma(*i, gamma_shape, gamma_shape);
      for (MyMatrix::const_iterator i = param.extraVarPow.begin();
	   i != param.extraVarPow.end();
	   ++i)
	log_prob += log_dgamma(*i, gamma_shape, gamma_shape);
      
      /* exo */

      return log_prob;
    }

  double getDLogProb(const OUparam& param, const OUparam& Dparam) const
    {
      double log_prob = 0.0;
      
      /* mean */
      
      /* extraMean */
      
      /* reversion */
      /* arma::cx_vec eigval = arma::eig_gen(param.reversion); */
      /* for (arma::cx_vec::const_iterator ix = eigval.begin(); */
      /* 	   ix != eigval.end(); */
      /* 	   ix++) */
      /* 	{ */
      /* 	  log_prob += log_f(real(*ix), ...); */
      /* 	}	 */
      log_prob -= arma::trace(Dparam.reversion);
      
      /* sigma */
      log_prob += D_log_dinvwish(param.sigma, Dparam.sigma, sigma_prior_df_plus, sigma_prior_scale / (sigma_prior_df_plus + param.sigma.n_rows));
      
      for (uword i = 0; i < param.extraVarAdd.size(); ++i)
	log_prob += D_log_dexp(param.extraVarAdd[i], Dparam.extraVarAdd[i], 1.0);
      for (uword i = 0; i < param.extraVarMult.size(); ++i)
	log_prob += D_log_dgamma(param.extraVarMult[i], Dparam.extraVarMult[i], gamma_shape, gamma_shape);
      for (uword i = 0; i < param.extraVarPow.size(); ++i)
	log_prob += D_log_dgamma(param.extraVarPow[i], Dparam.extraVarPow[i], gamma_shape, gamma_shape);
      
      /* exo */

      return log_prob;
    }
  
  double getDDLogProb(const OUparam& param, const OUparam& Dparam, const OUparam& DDparam) const
    {
      double log_prob = 0.0;
      
      /* mean */
      
      /* extraMean */
      
      /* reversion */
      /* arma::cx_vec eigval = arma::eig_gen(param.reversion); */
      /* for (arma::cx_vec::const_iterator ix = eigval.begin(); */
      /* 	   ix != eigval.end(); */
      /* 	   ix++) */
      /* 	{ */
      /* 	  log_prob += log_f(real(*ix), ...); */
      /* 	}	 */
      log_prob -= arma::trace(DDparam.reversion);
      
      /* sigma */
      log_prob += DD_log_dinvwish(param.sigma, Dparam.sigma, DDparam.sigma, sigma_prior_df_plus, sigma_prior_scale / (sigma_prior_df_plus + param.sigma.n_rows));
      
      for (uword i = 0; i < param.extraVarAdd.size(); ++i)
	log_prob += DD_log_dexp(param.extraVarAdd[i], Dparam.extraVarAdd[i], DDparam.extraVarAdd[i], 1.0);
      for (uword i = 0; i < param.extraVarMult.size(); ++i)
	log_prob += DD_log_dgamma(param.extraVarMult[i], Dparam.extraVarMult[i], DDparam.extraVarMult[i], gamma_shape, gamma_shape);
      for (uword i = 0; i < param.extraVarPow.size(); ++i)
	log_prob += DD_log_dgamma(param.extraVarPow[i], Dparam.extraVarPow[i], DDparam.extraVarPow[i], gamma_shape, gamma_shape);
      
      /* exo */

      return log_prob;
    }
  
 private:
  double sigma_prior_df_plus;
  double sigma_prior_scale;
  double gamma_shape;
};


class OUparamspace{
public:
  OUparamspace(const size_t iProcessIndexSize,
           const size_t iHouseSize,
           const bool iHasStart,
           const bool iHasExo,
           const bool iVarByHouse,
           const inputVector& iVarCenters,
           const inputMatrix& iStartConstraint,
           const inputMatrix& iMeanConstraint,
           const inputMatrix& iExtraMeanConstraint,
           const inputMatrix& iReversionConstraint,
           const inputMatrix& iSimultaneousConstraint,
           const inputMatrix& iSigmaSqrtConstraint,
           const inputMatrix& iExtraVarAddConstraint,
           const inputMatrix& iExtraVarMultConstraint,
           const inputMatrix& iExtraVarPowConstraint,
           const inputMatrix& iExoConstraint)
  :processIndexSize(iProcessIndexSize)
  ,houseSize(iHouseSize)
  ,hasStart(iHasStart)
  ,hasExo(iHasExo)
  ,varByHouse(iVarByHouse)
  ,varCenters(RtoMatrixCol(iVarCenters))
  ,startConstraint(RtoMatrix(iStartConstraint))
  ,meanConstraint(RtoMatrix(iMeanConstraint))
  ,extraMeanConstraint(RtoMatrix(iExtraMeanConstraint))
  ,reversionConstraint(RtoMatrix(iReversionConstraint))
  ,simultaneousConstraint(RtoMatrix(iSimultaneousConstraint))
  ,sigmaSqrtConstraint(RtoMatrix(iSigmaSqrtConstraint))
  ,extraVarAddConstraint(RtoMatrix(iExtraVarAddConstraint))
  ,extraVarMultConstraint(RtoMatrix(iExtraVarMultConstraint))
  ,extraVarPowConstraint(RtoMatrix(iExtraVarPowConstraint))
  ,exoConstraint(RtoMatrix(iExoConstraint))
  ,paramSize(iMeanConstraint.INPUTNCOL)
  {
#ifndef MONOCAR_NO_DEBUG
    assert(iReversionConstraint.INPUTNCOL == paramSize);
    assert(iSimultaneousConstraint.INPUTNCOL == paramSize);
    assert(iSigmaSqrtConstraint.INPUTNCOL == paramSize);
    assert(iExtraVarAddConstraint.INPUTNCOL == paramSize);
    assert(iExtraVarMultConstraint.INPUTNCOL == paramSize);
    assert(iExtraVarPowConstraint.INPUTNCOL == paramSize);
    assert(uword(iReversionConstraint.INPUTNROW) ==
       uword(std::sqrt(iReversionConstraint.INPUTNROW))*
       uword(std::sqrt(iReversionConstraint.INPUTNROW)));
    assert(uword(iSimultaneousConstraint.INPUTNROW) ==
       uword(std::sqrt(iSimultaneousConstraint.INPUTNROW))*
       uword(std::sqrt(iSimultaneousConstraint.INPUTNROW)));
#endif
  }

  MyMatrix multProcess(const MyMatrix& conmat, const MyMatrix& parmat, const size_t nrow) const
  {
    MyMatrix out = conmat * parmat;
    if (out.n_rows == nrow * processIndexSize)
        out.reshape((uword) nrow, (uword) processIndexSize);
    return out;
  }
  
  OUparam getParams(const MyMatrix& origParams,
            bool& fail) const
  {
    OUparam newParams;
    
    fail = false;
    
#ifndef MONOCAR_NO_DEBUG
    assert(origParams.n_cols == 1);
#endif
    
    MyMatrix params(origParams);
    
    params.resize(params.n_rows + 1, 1);
    params(origParams.n_rows, 0) = 1.0;

    newParams.varByHouse = varByHouse;
    newParams.varCenters = varCenters;

    newParams.extraVarAdd = arma::square(extraVarAddConstraint * params);
    newParams.extraVarMult = arma::square(extraVarMultConstraint * params);
    newParams.extraVarPow = extraVarPowConstraint * params;

    newParams.reversion = reversionConstraint * params;
    newParams.reversion.reshape((uword) std::sqrt(newParams.reversion.n_rows),
                (uword) std::sqrt(newParams.reversion.n_rows));

    newParams.mean = multProcess(meanConstraint, params, newParams.reversion.n_rows);
    newParams.extraMean = multProcess(extraMeanConstraint, params, newParams.reversion.n_rows);

    if (hasStart)
      newParams.start = multProcess(startConstraint, params, newParams.reversion.n_rows);
    else
      newParams.start = multProcess(meanConstraint, params, newParams.reversion.n_rows);
    newParams.hasStart = hasStart;

//    if (newParams.mean.n_rows == newParams.reversion.n_rows * processIndexSize)
//        newParams.mean.reshape((uword) newParams.reversion.n_rows, (uword) processIndexSize);
//    if (newParams.start.n_rows == newParams.reversion.n_rows * processIndexSize)
//        newParams.start.reshape((uword) newParams.reversion.n_rows, (uword) processIndexSize);
//    if (newParams.extraMean.n_rows == ((uword) houseSize) * ((uword) processIndexSize))
//        newParams.extraMean.reshape((uword) houseSize, (uword) processIndexSize);

    newParams.hasExo = hasExo;
    if ((hasExo) && (exoConstraint.n_rows > 0))
      {
    newParams.exo = exoConstraint * params;
    newParams.hasExo = true;
      }
    else
      {
    newParams.hasExo = false;
      }

    MyMatrix sigmaSqrtUpper = sigmaSqrtConstraint * params;
    MyMatrix sigmaSqrt(newParams.reversion.n_rows, newParams.reversion.n_cols, fill::zeros);      

    uword pos = 0;
    for (uword j = 0; j < sigmaSqrt.n_rows; j++)
      {
    for (uword i = 0; i <= j; i++)
      {
        sigmaSqrt(j,i) = sigmaSqrtUpper(pos, 0);
        pos++;
      }
      }

    MyMatrix simultaneous = simultaneousConstraint * params;
    simultaneous.reshape((uword) std::sqrt(simultaneous.n_rows),
             (uword) std::sqrt(simultaneous.n_rows));

    newParams.reversion = solve(simultaneous, newParams.reversion);
    sigmaSqrt = solve(simultaneous, sigmaSqrt);
    newParams.sigma = sigmaSqrt * sigmaSqrt.t();

    if (hasNegEigen(newParams.reversion)) {
      if (debugLevel > 2)
    COUT << endl << "Matrix not valid:" << newParams.reversion << endl;
      else if (debugLevel == 2)
    COUT << "Matrix not valid" << endl;
      fail = true;
    }

    return newParams;
  }



  OUparam getDParams(const MyMatrix& origParams,
             const uword i) const
  {
    OUparam newParams;

    MyMatrix params(origParams.n_rows + 1, 1, fill::zeros);
    params(i, 0) = 1.0;

    MyMatrix realParams = origParams;
    realParams.resize(origParams.n_rows + 1, 1);
    realParams(origParams.n_rows, 0) = 1.0;

    newParams.varByHouse = varByHouse;
    newParams.varCenters = varCenters;

    MyMatrix realExtraVarAdd = extraVarAddConstraint * realParams;
    newParams.extraVarAdd = 2.0 * realExtraVarAdd % (extraVarAddConstraint * params);
    MyMatrix realExtraVarMult = extraVarMultConstraint * realParams;
    newParams.extraVarMult = 2.0 * realExtraVarMult % (extraVarMultConstraint * params);
    newParams.extraVarPow = extraVarPowConstraint * params;

    newParams.reversion = reversionConstraint * params;
    newParams.reversion.reshape((uword) std::sqrt(newParams.reversion.n_rows),
                (uword) std::sqrt(newParams.reversion.n_rows));

    newParams.mean = multProcess(meanConstraint, params, newParams.reversion.n_rows);
    newParams.extraMean = multProcess(extraMeanConstraint, params, newParams.reversion.n_rows);

    if (hasStart)
        newParams.start = multProcess(startConstraint, params, newParams.reversion.n_rows);
    else
        newParams.start = multProcess(meanConstraint, params, newParams.reversion.n_rows);
    newParams.hasStart = hasStart;

//    if (newParams.mean.n_rows == newParams.reversion.n_rows * processIndexSize)
//        newParams.mean.reshape((uword) newParams.reversion.n_rows, (uword) processIndexSize);
//    if (newParams.start.n_rows == newParams.reversion.n_rows * processIndexSize)
//        newParams.start.reshape((uword) newParams.reversion.n_rows, (uword) processIndexSize);
//    if (newParams.extraMean.n_rows == ((uword) houseSize) * ((uword) processIndexSize))
//        newParams.extraMean.reshape((uword) houseSize, (uword) processIndexSize);

    newParams.hasExo = hasExo;
    if ((hasExo) && (exoConstraint.n_rows > 0))
      {
    newParams.exo = exoConstraint * params;
    newParams.hasExo = true;
      }
    else
      {
    newParams.hasExo = false;
      }

    MyMatrix sigmaSqrtUpper = sigmaSqrtConstraint * params;
    MyMatrix sigmaSqrt(newParams.reversion.n_rows, newParams.reversion.n_cols, fill::zeros);      

    MyMatrix realSigmaSqrtUpper = sigmaSqrtConstraint * realParams;
    MyMatrix realSigmaSqrt(newParams.reversion.n_rows, newParams.reversion.n_cols, fill::zeros);      

    uword pos = 0;
    for (uword j = 0; j < sigmaSqrt.n_rows; j++)
      {
    for (uword k = 0; k <= j; k++)
      {
        sigmaSqrt(j,k) = sigmaSqrtUpper(pos, 0);
        realSigmaSqrt(j,k) = realSigmaSqrtUpper(pos, 0);
        pos++;
      }
      }

    MyMatrix simultaneous = simultaneousConstraint * params;
    simultaneous.reshape((uword) std::sqrt(simultaneous.n_rows),
             (uword) std::sqrt(simultaneous.n_rows));

    MyMatrix realSimultaneous = simultaneousConstraint * realParams;
    realSimultaneous.reshape((uword) std::sqrt(realSimultaneous.n_rows),
                 (uword) std::sqrt(realSimultaneous.n_rows));
    MyMatrix solveRealSimultaneous = inv(realSimultaneous);
    MyMatrix DsolveSimultaneous = - solveRealSimultaneous * simultaneous * solveRealSimultaneous;

    MyMatrix realReversion = reversionConstraint * realParams;
    realReversion.reshape((uword) std::sqrt(realReversion.n_rows),
              (uword) std::sqrt(realReversion.n_rows));

    newParams.reversion = DsolveSimultaneous * realReversion + solveRealSimultaneous * newParams.reversion;
    sigmaSqrt = DsolveSimultaneous * realSigmaSqrt + solveRealSimultaneous * sigmaSqrt;
    realSigmaSqrt = solve(realSimultaneous, realSigmaSqrt);
    newParams.sigma = sigmaSqrt * realSigmaSqrt.t() + realSigmaSqrt * sigmaSqrt.t();

    return newParams;

  }

  OUparam getDDParams(const MyMatrix& origParams,
              const uword i1,
              const uword i2) const
  {
    OUparam newParams;

    MyMatrix params(origParams.n_rows + 1, 1, fill::zeros);

    MyMatrix params1(origParams.n_rows + 1, 1, fill::zeros);
    params1(i1, 0) = 1.0;
    MyMatrix params2(origParams.n_rows + 1, 1, fill::zeros);
    params2(i2, 0) = 1.0;
    
    MyMatrix realParams = origParams;
    realParams.resize(origParams.n_rows + 1, 1);
    realParams(origParams.n_rows, 0) = 1.0;

    newParams.varByHouse = varByHouse;
    newParams.varCenters = varCenters;

    MyMatrix realExtraVarAdd = extraVarAddConstraint * realParams;
    newParams.extraVarAdd = 2.0 * (extraVarAddConstraint * params1) % (extraVarAddConstraint * params2) + 2.0 * realExtraVarAdd % (extraVarAddConstraint * params);
    MyMatrix realExtraVarMult = extraVarMultConstraint * realParams;
    newParams.extraVarMult = 2.0 * (extraVarMultConstraint * params1) % (extraVarMultConstraint * params2) + 2.0 * realExtraVarMult % (extraVarMultConstraint * params);
    newParams.extraVarPow = extraVarPowConstraint * params;

    newParams.reversion = reversionConstraint * params;
    newParams.reversion.reshape((uword) std::sqrt(newParams.reversion.n_rows),
                (uword) std::sqrt(newParams.reversion.n_rows));

    MyMatrix DDreversion = reversionConstraint * params;
    DDreversion.reshape((uword) std::sqrt(DDreversion.n_rows),
            (uword) std::sqrt(DDreversion.n_rows));

    MyMatrix D1reversion = reversionConstraint * params1;
    D1reversion.reshape((uword) std::sqrt(D1reversion.n_rows),
            (uword) std::sqrt(D1reversion.n_rows));
    
    MyMatrix D2reversion = reversionConstraint * params2;
    D2reversion.reshape((uword) std::sqrt(D2reversion.n_rows),
            (uword) std::sqrt(D2reversion.n_rows));

    MyMatrix realReversion = reversionConstraint * realParams;
    realReversion.reshape((uword) std::sqrt(realReversion.n_rows),
                  (uword) std::sqrt(realReversion.n_rows));

    newParams.mean = multProcess(meanConstraint, params, newParams.reversion.n_rows);
    newParams.extraMean = multProcess(extraMeanConstraint, params, newParams.reversion.n_rows);

    if (hasStart)
        newParams.start = multProcess(startConstraint, params, newParams.reversion.n_rows);
    else
        newParams.start = multProcess(meanConstraint, params, newParams.reversion.n_rows);
    newParams.hasStart = hasStart;

//    if (newParams.mean.n_rows == DDreversion.n_rows * processIndexSize)
//        newParams.mean.reshape((uword) DDreversion.n_rows, (uword) processIndexSize);
//    if (newParams.start.n_rows == DDreversion.n_rows * processIndexSize)
//        newParams.start.reshape((uword) DDreversion.n_rows, (uword) processIndexSize);
//    if (newParams.extraMean.n_rows == ((uword) houseSize) * ((uword) processIndexSize))
//        newParams.extraMean.reshape((uword) houseSize, (uword) processIndexSize);

    newParams.hasExo = hasExo;
    if ((hasExo) && (exoConstraint.n_rows > 0))
      {
    newParams.exo = exoConstraint * params;
    newParams.hasExo = true;
      }
    else
      {
    newParams.hasExo = false;
      }

    MyMatrix DDsigmaSqrtUpper = sigmaSqrtConstraint * params;
    MyMatrix DDsigmaSqrt(DDreversion.n_rows, DDreversion.n_cols, fill::zeros);

    MyMatrix D1sigmaSqrtUpper = sigmaSqrtConstraint * params1;
    MyMatrix D1sigmaSqrt(DDreversion.n_rows, DDreversion.n_cols, fill::zeros);

    MyMatrix D2sigmaSqrtUpper = sigmaSqrtConstraint * params2;
    MyMatrix D2sigmaSqrt(DDreversion.n_rows, DDreversion.n_cols, fill::zeros);
    
    MyMatrix realSigmaSqrtUpper = sigmaSqrtConstraint * realParams;
    MyMatrix realSigmaSqrt(DDreversion.n_rows, DDreversion.n_cols, fill::zeros);

    uword pos = 0;
    for (uword j = 0; j < DDsigmaSqrt.n_rows; j++)
      {
    for (uword i = 0; i <= j; i++)
      {
        DDsigmaSqrt(j,i) = DDsigmaSqrtUpper(pos, 0);
        D1sigmaSqrt(j,i) = D1sigmaSqrtUpper(pos, 0);
        D2sigmaSqrt(j,i) = D2sigmaSqrtUpper(pos, 0);
        realSigmaSqrt(j,i) = realSigmaSqrtUpper(pos, 0);
        pos++;
      }
      }

    MyMatrix DDsimultaneous = simultaneousConstraint * params;
    DDsimultaneous.reshape((uword) std::sqrt(DDsimultaneous.n_rows),
                 (uword) std::sqrt(DDsimultaneous.n_rows));

    MyMatrix D1simultaneous = simultaneousConstraint * params1;
    D1simultaneous.reshape((uword) std::sqrt(D1simultaneous.n_rows),
                 (uword) std::sqrt(D1simultaneous.n_rows));

    MyMatrix D2simultaneous = simultaneousConstraint * params2;
    D2simultaneous.reshape((uword) std::sqrt(D2simultaneous.n_rows),
                 (uword) std::sqrt(D2simultaneous.n_rows));

    MyMatrix realSimultaneous = simultaneousConstraint * realParams;
    realSimultaneous.reshape((uword) std::sqrt(realSimultaneous.n_rows),
                     (uword) std::sqrt(realSimultaneous.n_rows));
    MyMatrix solveRealSimultaneous = inv(realSimultaneous);
    MyMatrix D1solveSimultaneous = - solveRealSimultaneous * D1simultaneous * solveRealSimultaneous;
    MyMatrix D2solveSimultaneous = - solveRealSimultaneous * D2simultaneous * solveRealSimultaneous;
    MyMatrix DDsolveSimultaneous =
      + solveRealSimultaneous * D2simultaneous * solveRealSimultaneous * D1simultaneous * solveRealSimultaneous
      - solveRealSimultaneous * DDsimultaneous * solveRealSimultaneous
      + solveRealSimultaneous * D1simultaneous * solveRealSimultaneous * D2simultaneous * solveRealSimultaneous;

    DDreversion =
      + DDsolveSimultaneous * realReversion
      + D1solveSimultaneous * D2reversion
      + D2solveSimultaneous * D1reversion
      + solveRealSimultaneous * DDreversion;
    
    DDsigmaSqrt =
      + DDsolveSimultaneous * realSigmaSqrt
      + D1solveSimultaneous * D2sigmaSqrt
      + D2solveSimultaneous * D1sigmaSqrt
      + solveRealSimultaneous * DDsigmaSqrt;
    D1sigmaSqrt = D1solveSimultaneous * realSigmaSqrt + solveRealSimultaneous * D1sigmaSqrt;
    D2sigmaSqrt = D2solveSimultaneous * realSigmaSqrt + solveRealSimultaneous * D2sigmaSqrt;
    realSigmaSqrt = solveRealSimultaneous * realSigmaSqrt;
    
    //D1sigmaSqrt = D1solveSimultaneous * realSigmaSqrt + solve(realSimultaneous, D1sigmaSqrt);
    //D2sigmaSqrt = D2solveSimultaneous * realSigmaSqrt + solve(realSimultaneous, D2sigmaSqrt);
    //realSigmaSqrt = solve(realSimultaneous, realSigmaSqrt);

    newParams.reversion = DDreversion;

    newParams.sigma =
      + DDsigmaSqrt * realSigmaSqrt.t()
      + D1sigmaSqrt * D2sigmaSqrt.t()
      + D2sigmaSqrt * D1sigmaSqrt.t()
      + realSigmaSqrt * DDsigmaSqrt.t();

    return newParams;

  }

private:
  const size_t processIndexSize;
  const size_t houseSize;
  const bool hasStart;
  const bool hasExo;
  const bool varByHouse;
  const MyMatrixCol varCenters;
  const MyMatrix startConstraint;
  const MyMatrix meanConstraint;
  const MyMatrix extraMeanConstraint;
  const MyMatrix reversionConstraint;
  const MyMatrix simultaneousConstraint;
  const MyMatrix sigmaSqrtConstraint;
  const MyMatrix extraVarAddConstraint;
  const MyMatrix extraVarMultConstraint;
  const MyMatrix extraVarPowConstraint;
  const MyMatrix exoConstraint;
  
  const size_t paramSize;
};

class OUprocess
{
public:
  OUprocess(const OUparam& iParam,
        const vector<OUparam>& iDparam,
        const vector<OUparam>& iDDparam,
//        const double timeZero,
        bool iUseLong,
        const bool iPartialScores = false)
    :partialScores(iPartialScores)
    ,logLik(0.0)
    ,seriesSize(iParam.mean.n_rows)
    ,params(iParam)
    ,Dparams(iDparam)
    ,DDparams(iDDparam)
    ,mu(iParam.mean - iParam.mean)
    ,P(iParam.reversion)
    ,Sigma(iParam.sigma)
    ,Ident(eye(iParam.reversion.n_rows, iParam.reversion.n_rows))
    ,PpP(kron(iParam.reversion,eye(iParam.reversion.n_rows, iParam.reversion.n_rows)) + kron(eye(iParam.reversion.n_rows, iParam.reversion.n_rows),iParam.reversion))
    ,sPpP(inv(PpP))
    ,sP(inv(iParam.reversion))
    ,expP(iParam.reversion)
    ,covHist(0, 0, fill::zeros)
    ,covHC(0, iParam.reversion.n_rows, fill::zeros)
    ,covCurrent(iParam.reversion.n_rows, iParam.reversion.n_rows, fill::zeros)
    ,DlogLik(iDparam.size(), partialScores ? 0 : 1, fill::zeros)
    ,DexpP(iParam.reversion)
    ,DDlogLik(iDparam.size(), iDparam.size(), fill::zeros)
    ,DDexpP(iParam.reversion)
    ,lastTime(numeric_limits<double>::infinity())
    ,smoother(iUseLong)
    ,hasHess((iDDparam.size() > 0) && (iDDparam.size() == (iDparam.size() * (iDparam.size()+1) / 2)))
  {
//    reset(timeZero);
    for (vector<OUparam>::const_iterator i = Dparams.begin(); i != Dparams.end(); ++i)
      {
        DexpP.push_back(i->reversion);

        DcovHist.push_back(MyMatrix(0, 0));
        DcovHC.push_back(MyMatrix(0, params.reversion.n_rows));
        MyMatrix Dcov33(params.reversion.n_rows, params.reversion.n_rows, fill::zeros);

        const MyMatrix& DP = i->reversion;
        const MyMatrix& DSigma = i->sigma;
        MyMatrix DStart = (i->start) - (i->mean);

        DPpP.push_back(kron(DP,Ident) + kron(Ident,DP));
        DsPpP.push_back(- sPpP * DPpP.back() * sPpP);
        DsP.push_back(- sP * DP * sP);
        Dmu.push_back(i->mean - i->mean);

        DmeanHist.push_back(MyMatrix(0,DStart.n_cols));
        DmeanCurrent.push_back(DStart);

        if (!params.hasStart)
          {
            Dcov33 =  DsPpP.back() * vectorise(Sigma) + sPpP * vectorise(DSigma);
            Dcov33.reshape(P.n_rows, P.n_cols);
          }

        DcovCurrent.push_back(Dcov33);
      }

    if (hasHess)
      {
        uword k = 0;
        uword icount = 0;
        for (vector<OUparam>::const_iterator i = Dparams.begin(); i != Dparams.end(); ++i)
          {
            icount++;
            uword jcount = icount;

            const MyMatrix& D1P = i->reversion;
            const MyMatrix& D1Sigma = i->sigma;
            MyMatrix D1Start = (i->start) - (i->mean);
            MyMatrix D1PpP = kron(D1P,Ident) + kron(Ident,D1P);
            MyMatrix D1sPpP = - sPpP * D1PpP * sPpP;

            for (vector<OUparam>::const_iterator j = i; j != Dparams.end(); ++j)
              {
                jcount++;
                DDexpP.push_back(i->reversion, j->reversion, DDparams[k].reversion);

                DDcovHist.push_back(MyMatrix(0, 0));
                DDcovHC.push_back(MyMatrix(0, params.reversion.n_rows));
                MyMatrix DDcov33(params.reversion.n_rows, params.reversion.n_rows, fill::zeros);

                const MyMatrix& D2P = j->reversion;
                const MyMatrix& D1D2P = DDparams[k].reversion;
                const MyMatrix& D2Sigma = j->sigma;
                const MyMatrix& D1D2Sigma = DDparams[k].sigma;
                MyMatrix D2Start = (j->start) - (j->mean);
                MyMatrix DDStart = DDparams[k].start - DDparams[k].mean;;

                DDPpP.push_back(kron(D1D2P,Ident) + kron(Ident,D1D2P));
                MyMatrix D2PpP = kron(D2P,Ident) + kron(Ident,D2P);
                DDsPpP.push_back(sPpP * D2PpP * sPpP * D1PpP * sPpP
                                 - sPpP * DDPpP.back() * sPpP
                                 + sPpP * D1PpP * sPpP * D2PpP * sPpP);
                MyMatrix D2sPpP = - sPpP * D2PpP * sPpP;
                DDsP.push_back(+ sP * D2P * sP * D1P * sP
                               - sP * D1D2P * sP
                               + sP * D1P * sP * D2P * sP);
                DDmu.push_back(DDparams[k].mean - DDparams[k].mean);

                DDmeanHist.push_back(MyMatrix(0,DDStart.n_cols));
                DDmeanCurrent.push_back(DDStart);

                if (!params.hasStart)
                  {
                    DDcov33 =
                    + DDsPpP.back() * vectorise(Sigma)
                    + D1sPpP * vectorise(D2Sigma)
                    + D2sPpP * vectorise(D1Sigma)
                    + sPpP * vectorise(D1D2Sigma);
                    DDcov33.reshape(P.n_rows, P.n_cols);
                  }

                DDcovCurrent.push_back(DDcov33);
                
                k++;
              }
          }
      }
  }

  void reset(double timeZero)
  {
    meanHist.set_size(0,1);
    covHist.set_size(0,0);
    covHC.set_size(0, params.reversion.n_rows);

    meanCurrent = params.start - params.mean;

    if (!params.hasStart)
      {
//        MyMatrix cov33 =  sPpP * vectorise(Sigma);
//        cov33.reshape(P.n_rows, P.n_cols);
//        covCurrent = cov33;
        covCurrent = sPpP * vectorise(Sigma);
        covCurrent.reshape(P.n_rows, P.n_cols);
      }

    timeIndices.clear();
    timeIndices.push_back(timeZero);

    smoother.reset();
  }

  MyMatrixCol sampleCurrent(void)
  {
    MyMatrixCol obs = rmvnorm(meanCurrent, covCurrent);
    meanCurrent = obs;
    meanHist += covHC * solve(covCurrent, obs - meanCurrent);
    covHist -= covHC * solve(covCurrent, covHC.t());
    covCurrent.zeros();
    covHC.zeros();
    return(obs);
  }

  MyMatrixCol sampleHist(void)
  {
    MyMatrixCol obs = rmvnorm(meanHist, covHist);
    meanHist = obs;
    meanCurrent += covHC.t() * solve(covHist, obs - meanHist);
    covCurrent -= covHC.t() * solve(covHist, covHC);
    covHist.zeros();
    covHC.zeros();
    return(obs);
  }

  void smooth(vector<double>& times, MyCube& means, MyCube& vcovs)
  {
#ifdef RTS_SMOOTHER
    smoother.terminate(meanCurrent, meanHist, covCurrent, covHist, covHC);
#endif
    smoother.smoothBackward(times, means, vcovs);
  }

  void record()
  {
    smoother.addRecord(timeIndices[timeIndices.size()-1],
               meanCurrent, meanHist,
               covCurrent, covHist, covHC);
  }

    void printTimes()
    {
      COUT << std::endl << meanCurrent.t() << meanHist.t();
      for (size_t i = 0; i < timeIndices.size(); ++i)
          COUT << "\t" << timeIndices[i];
    }

  void addTime(const double timeEnd)
  {
#ifndef MONOCAR_NO_DEBUG
    assert(timeEnd > timeIndices[timeIndices.size()-1]);
#endif
    
    double time = timeEnd - timeIndices[timeIndices.size()-1];

    timeIndices.push_back(timeEnd);

    if (time != lastTime)
      {
    expPnegtime = expP(-time);
    IexpPpPnegtime = -kron(expPnegtime,expPnegtime);
    IexpPpPnegtime.diag() += 1.0;
    lastTime = time;
      }

    MyMatrix PinvIexpPntot = (sP * (Ident-expPnegtime)) / time;

    
    MyMatrix cov33 = sPpP * (IexpPpPnegtime) * vectorise(Sigma);
    cov33.reshape(P.n_rows, P.n_cols);


    MyMatrix cov12 = covHC * PinvIexpPntot.t();


    MyMatrix cov22 = sP * (Ident - PinvIexpPntot) * Sigma
      - PinvIexpPntot * Sigma * sP.t()
      + cov33 * sP.t() / time;
    cov22.reshape(P.n_rows * P.n_cols, 1);

    MyMatrix cov23 = Sigma * PinvIexpPntot.t()
      - PinvIexpPntot * Sigma * expPnegtime.t();
    cov23.reshape(P.n_rows*P.n_cols, 1);


    if (hasHess)
      {
    uword k = 0;
    for (uword i = 0; i < Dparams.size(); i++)
      {
        MyMatrix D1expPnegtime = DexpP(-time, i);
        MyMatrix D1IexpPpPnegtime =
          -kron(D1expPnegtime,expPnegtime)
          -kron(expPnegtime,D1expPnegtime);
        MyMatrix D1PinvIexpPntot = (DsP[i] * (Ident-expPnegtime) - sP * D1expPnegtime) / time;
        MyMatrix D1Sigma = Dparams[i].sigma;
        
        MyMatrix D1cov33 =
          DsPpP[i] * (IexpPpPnegtime) * vectorise(Sigma)
          + sPpP * (D1IexpPpPnegtime) * vectorise(Sigma)
          + sPpP * (IexpPpPnegtime) * vectorise(D1Sigma);
        D1cov33.reshape(P.n_rows, P.n_cols);
        
        MyMatrix D1cov22 = DsP[i] * (Ident - PinvIexpPntot) * Sigma
          - sP * (D1PinvIexpPntot) * Sigma
          + sP * (Ident - PinvIexpPntot) * D1Sigma
          - D1PinvIexpPntot * Sigma * sP.t()
          - PinvIexpPntot * D1Sigma * sP.t()
          - PinvIexpPntot * Sigma * DsP[i].t()
          + D1cov33 * sP.t() / time
          + cov33 * DsP[i].t() / time;
        D1cov22.reshape(P.n_rows * P.n_cols, 1);
        // D1cov22 = (DsPpP[i] * cov22 + sPpP * D1cov22) / time;
        // D1cov22.reshape(P.n_rows, P.n_cols);
        // D1cov22 += D1cov22.t();
        
        MyMatrix D1cov23 =
          D1Sigma * PinvIexpPntot.t()
          + Sigma * D1PinvIexpPntot.t()
          - D1PinvIexpPntot * Sigma * expPnegtime.t()
          - PinvIexpPntot * D1Sigma * expPnegtime.t()
          - PinvIexpPntot * Sigma * D1expPnegtime.t();
        D1cov23.reshape(P.n_rows*P.n_cols, 1);
        // D1cov23 = DsPpP[i] * cov23 + sPpP * D1cov23;
        // D1cov23.reshape(P.n_rows, P.n_cols);

        for (uword j = i; j < Dparams.size(); j++)
          {
        
        
        MyMatrix D2expPnegtime = DexpP(-time, j);
        MyMatrix D2IexpPpPnegtime =
          -kron(D2expPnegtime,expPnegtime)
          -kron(expPnegtime,D2expPnegtime);
        MyMatrix PinvIexpPntot = (sP * (Ident-expPnegtime)) / time;
        MyMatrix D2PinvIexpPntot = (DsP[j] * (Ident-expPnegtime) - sP * D2expPnegtime) / time;
        MyMatrix D2Sigma = Dparams[j].sigma;
        
        MyMatrix D1D2expPnegtime = DDexpP(-time, k);
        MyMatrix D1D2IexpPpPnegtime =
          -kron(D1D2expPnegtime,expPnegtime)
          -kron(D1expPnegtime,D2expPnegtime)
          -kron(D2expPnegtime,D1expPnegtime)
          -kron(expPnegtime,D1D2expPnegtime);
        MyMatrix D1D2PinvIexpPntot = (DDsP[k] * (Ident-expPnegtime)
                          - DsP[i] * D2expPnegtime
                          - DsP[j] * D1expPnegtime
                          - sP * D1D2expPnegtime) / time;
        MyMatrix D1D2Sigma = DDparams[k].sigma;

        MyMatrix D2cov33 =
          DsPpP[j] * (IexpPpPnegtime) * vectorise(Sigma)
          + sPpP * (D2IexpPpPnegtime) * vectorise(Sigma)
          + sPpP * (IexpPpPnegtime) * vectorise(D2Sigma);
        D2cov33.reshape(P.n_rows, P.n_cols);
        
        MyMatrix DDcov33 =
          + DDsPpP[k] * (IexpPpPnegtime) * vectorise(Sigma)
          + DsPpP[j] * (D1IexpPpPnegtime) * vectorise(Sigma)
          + DsPpP[j] * (IexpPpPnegtime) * vectorise(D1Sigma)
          + DsPpP[i] * (D2IexpPpPnegtime) * vectorise(Sigma)
          + sPpP * (D1D2IexpPpPnegtime) * vectorise(Sigma)
          + sPpP * (D2IexpPpPnegtime) * vectorise(D1Sigma)
          + DsPpP[i] * (IexpPpPnegtime) * vectorise(D2Sigma)
          + sPpP * (D1IexpPpPnegtime) * vectorise(D2Sigma)
          + sPpP * (IexpPpPnegtime) * vectorise(D1D2Sigma);
        DDcov33.reshape(P.n_rows, P.n_cols);

        MyMatrix D2cov22 = DsP[j] * (Ident - PinvIexpPntot) * Sigma
          - sP * (D2PinvIexpPntot) * Sigma
          + sP * (Ident - PinvIexpPntot) * D2Sigma
          - D2PinvIexpPntot * Sigma * sP.t()
          - PinvIexpPntot * D2Sigma * sP.t()
          - PinvIexpPntot * Sigma * DsP[j].t()
          + D2cov33 * sP.t() / time
          + cov33 * DsP[j].t() / time;
        D2cov22.reshape(P.n_rows * P.n_cols, 1);
        // D2cov22 = (DsPpP[j] * cov22 + sPpP * D2cov22) / time;
        // D2cov22.reshape(P.n_rows, P.n_cols);
        // D2cov22 += D2cov22.t();
        
        
        MyMatrix DDcov22 =
          + DDsP[k] * (Ident - PinvIexpPntot) * Sigma
          - DsP[i] * D2PinvIexpPntot * Sigma
          + DsP[i] * (Ident - PinvIexpPntot) * D2Sigma
          
          - DsP[j] * (D1PinvIexpPntot) * Sigma
          - sP * (D1D2PinvIexpPntot) * Sigma
          - sP * (D1PinvIexpPntot) * D2Sigma
          
          + DsP[j] * (Ident - PinvIexpPntot) * D1Sigma
          - sP * D2PinvIexpPntot * D1Sigma
          + sP * (Ident - PinvIexpPntot) * D1D2Sigma
          
          - D1D2PinvIexpPntot * Sigma * sP.t()
          - D1PinvIexpPntot * D2Sigma * sP.t()
          - D1PinvIexpPntot * Sigma * DsP[j].t()
          
          - D2PinvIexpPntot * D1Sigma * sP.t()
          - PinvIexpPntot * D1D2Sigma * sP.t()
          - PinvIexpPntot * D1Sigma * DsP[j].t()
          
          - D2PinvIexpPntot * Sigma * DsP[i].t()
          - PinvIexpPntot * D2Sigma * DsP[i].t()
          - PinvIexpPntot * Sigma * DDsP[k].t()
          
          + DDcov33 * sP.t() / time
          + D1cov33 * DsP[j].t() / time          
          + D2cov33 * DsP[i].t() / time
          + cov33 * DDsP[k].t() / time;
        
        DDcov22.reshape(P.n_rows * P.n_cols, 1);
        DDcov22 = (DDsPpP[k] * cov22 + DsPpP[i] * D2cov22 + DsPpP[j] * D1cov22 + sPpP * DDcov22) / time;
        DDcov22.reshape(P.n_rows, P.n_cols);
        DDcov22 += DDcov22.t();
                
        MyMatrix D2cov23 =
          D2Sigma * PinvIexpPntot.t()
          + Sigma * D2PinvIexpPntot.t()
          - D2PinvIexpPntot * Sigma * expPnegtime.t()
          - PinvIexpPntot * D2Sigma * expPnegtime.t()
          - PinvIexpPntot * Sigma * D2expPnegtime.t();
        D2cov23.reshape(P.n_rows*P.n_cols, 1);
        // D2cov23 = DsPpP[j] * cov23 + sPpP * D2cov23;
        // D2cov23.reshape(P.n_rows, P.n_cols);
        
        MyMatrix DDcov23 =
          + D1D2Sigma * PinvIexpPntot.t()
          + D1Sigma * D2PinvIexpPntot.t()
          
          + D2Sigma * D1PinvIexpPntot.t()
          + Sigma * D1D2PinvIexpPntot.t()
          
          - D1D2PinvIexpPntot * Sigma * expPnegtime.t()
          - D1PinvIexpPntot * D2Sigma * expPnegtime.t()
          - D1PinvIexpPntot * Sigma * D2expPnegtime.t()
          
          - D2PinvIexpPntot * D1Sigma * expPnegtime.t()
          - PinvIexpPntot * D1D2Sigma * expPnegtime.t()
          - PinvIexpPntot * D1Sigma * D2expPnegtime.t()
          
          - D2PinvIexpPntot * Sigma * D1expPnegtime.t()
          - PinvIexpPntot * D2Sigma * D1expPnegtime.t()
          - PinvIexpPntot * Sigma * D1D2expPnegtime.t();
        
        DDcov23.reshape(P.n_rows*P.n_cols, 1);
        DDcov23 = DDsPpP[k] * cov23 + DsPpP[i] * D2cov23 + DsPpP[j] * D1cov23 + sPpP * DDcov23;
        DDcov23.reshape(P.n_rows, P.n_cols);
        
        // MyMatrix D1cov12 = DcovHC[i] * PinvIexpPntot.t() + covHC * D1PinvIexpPntot.t();
        // MyMatrix D2cov12 = DcovHC[j] * PinvIexpPntot.t() + covHC * D2PinvIexpPntot.t();        
        MyMatrix DDcov12 =
          + DDcovHC[k] * PinvIexpPntot.t()
          + DcovHC[i] * D2PinvIexpPntot.t()
          + DcovHC[j] * D1PinvIexpPntot.t()
          + covHC * D1D2PinvIexpPntot.t();

        DDmeanHist[k] = join_cols(DDmeanHist[k],
                      + D1D2PinvIexpPntot * (meanCurrent - mu)
                      + D1PinvIexpPntot * (DmeanCurrent[j] - Dmu[j])
                      + D2PinvIexpPntot * (DmeanCurrent[i] - Dmu[i])
                      + PinvIexpPntot * (DDmeanCurrent[k] - DDmu[k]));
        DDmeanCurrent[k] =
          + D1D2expPnegtime * (meanCurrent - mu)
          + D1expPnegtime * (DmeanCurrent[j] - Dmu[j])
          + D2expPnegtime * (DmeanCurrent[i] - Dmu[i])
          + expPnegtime * (DDmeanCurrent[k] - DDmu[k]);
        
        DDcovHC[k] = join_cols(DDcovHC[k] * expPnegtime.t()
                       + DcovHC[i] * D2expPnegtime.t()
                       + DcovHC[j] * D1expPnegtime.t()
                       + covHC * D1D2expPnegtime.t(),
                       DDcov23
                       
                       + D1D2PinvIexpPntot * covCurrent * expPnegtime.t()
                       + D1PinvIexpPntot * DcovCurrent[j] * expPnegtime.t()
                       + D1PinvIexpPntot * covCurrent * D2expPnegtime.t()
                       
                       + D2PinvIexpPntot * DcovCurrent[i] * expPnegtime.t()
                       + PinvIexpPntot * DDcovCurrent[k] * expPnegtime.t()
                       + PinvIexpPntot * DcovCurrent[i] * D2expPnegtime.t()
                   
                       + D2PinvIexpPntot * covCurrent * D1expPnegtime.t()
                       + PinvIexpPntot * DcovCurrent[j] * D1expPnegtime.t()
                       + PinvIexpPntot * covCurrent * D1D2expPnegtime.t()
                       );
        
        DDcovHist[k] = join_cols(join_rows(DDcovHist[k], DDcov12),
                     join_rows(DDcov12.t(),
                           DDcov22

                           + D1D2PinvIexpPntot * covCurrent * PinvIexpPntot.t()
                           + D1PinvIexpPntot * DcovCurrent[j] * PinvIexpPntot.t()
                           + D1PinvIexpPntot * covCurrent * D2PinvIexpPntot.t()
                           
                           + D2PinvIexpPntot * DcovCurrent[i] * PinvIexpPntot.t()
                           + PinvIexpPntot * DDcovCurrent[k] * PinvIexpPntot.t()
                           + PinvIexpPntot * DcovCurrent[i] * D2PinvIexpPntot.t()
                           
                           + D2PinvIexpPntot * covCurrent * D1PinvIexpPntot.t()
                           + PinvIexpPntot * DcovCurrent[j] * D1PinvIexpPntot.t()
                           + PinvIexpPntot * covCurrent * D1D2PinvIexpPntot.t()
                           ));
        
        DDcovCurrent[k] = DDcov33

          + D1D2expPnegtime * covCurrent * expPnegtime.t()
          + D1expPnegtime * DcovCurrent[j] * expPnegtime.t()
          + D1expPnegtime * covCurrent * D2expPnegtime.t()
          
          + D2expPnegtime * DcovCurrent[i] * expPnegtime.t()
          + expPnegtime * DDcovCurrent[k] * expPnegtime.t()
          + expPnegtime * DcovCurrent[i] * D2expPnegtime.t()
          
          + D2expPnegtime * covCurrent * D1expPnegtime.t()
          + expPnegtime * DcovCurrent[j] * D1expPnegtime.t()
          + expPnegtime * covCurrent * D1D2expPnegtime.t();
          
        k++;
          }
      }
      }


    for (uword i = 0; i < Dparams.size(); i++)
      {

    MyMatrix DexpPnegtime = DexpP(-time, i);

    MyMatrix DIexpPpPnegtime = -kron(DexpPnegtime,expPnegtime) -kron(expPnegtime,DexpPnegtime);
    MyMatrix DPinvIexpPntot = (DsP[i] * (Ident-expPnegtime) - sP * DexpPnegtime) / time;
    MyMatrix DSigma = Dparams[i].sigma;
        
    MyMatrix Dcov33 =
      DsPpP[i] * (IexpPpPnegtime) * vectorise(Sigma)
      + sPpP * (DIexpPpPnegtime) * vectorise(Sigma)
      + sPpP * (IexpPpPnegtime) * vectorise(DSigma);
    Dcov33.reshape(P.n_rows, P.n_cols);

    MyMatrix Dcov22 = DsP[i] * (Ident - PinvIexpPntot) * Sigma
      - sP * (DPinvIexpPntot) * Sigma
      + sP * (Ident - PinvIexpPntot) * DSigma
      - DPinvIexpPntot * Sigma * sP.t()
      - PinvIexpPntot * DSigma * sP.t()
      - PinvIexpPntot * Sigma * DsP[i].t()
      + Dcov33 * sP.t() / time
      + cov33 * DsP[i].t() / time;
    Dcov22.reshape(P.n_rows * P.n_cols, 1);
    Dcov22 = (DsPpP[i] * cov22 + sPpP * Dcov22) / time;
    Dcov22.reshape(P.n_rows, P.n_cols);
    Dcov22 += Dcov22.t();

    MyMatrix Dcov23 =
      DSigma * PinvIexpPntot.t()
      + Sigma * DPinvIexpPntot.t()
      - DPinvIexpPntot * Sigma * expPnegtime.t()
      - PinvIexpPntot * DSigma * expPnegtime.t()
      - PinvIexpPntot * Sigma * DexpPnegtime.t();
    Dcov23.reshape(P.n_rows*P.n_cols, 1);
    Dcov23 = DsPpP[i] * cov23 + sPpP * Dcov23;
    Dcov23.reshape(P.n_rows, P.n_cols);
    
    MyMatrix Dcov12 = DcovHC[i] * PinvIexpPntot.t() + covHC * DPinvIexpPntot.t();

    DmeanHist[i] = join_cols(DmeanHist[i], PinvIexpPntot * (DmeanCurrent[i] - Dmu[i]) + DPinvIexpPntot * (meanCurrent - mu) + Dmu[i]);
    
    DmeanCurrent[i] = expPnegtime * (DmeanCurrent[i] - Dmu[i]) + DexpPnegtime * (meanCurrent - mu) + Dmu[i];

    DcovHC[i] = join_cols(DcovHC[i] * expPnegtime.t()
               + covHC * DexpPnegtime.t(),
               Dcov23
               + DPinvIexpPntot * covCurrent * expPnegtime.t()
               + PinvIexpPntot * DcovCurrent[i] * expPnegtime.t()
               + PinvIexpPntot * covCurrent * DexpPnegtime.t() );

    DcovHist[i] = join_cols(join_rows(DcovHist[i], Dcov12),
                join_rows(Dcov12.t(),
                      Dcov22
                      + DPinvIexpPntot * covCurrent * PinvIexpPntot.t()
                      + PinvIexpPntot * DcovCurrent[i] * PinvIexpPntot.t()
                      + PinvIexpPntot * covCurrent * DPinvIexpPntot.t()
                      ));

    DcovCurrent[i] = Dcov33
      + DexpPnegtime * covCurrent * expPnegtime.t()
      + expPnegtime * DcovCurrent[i] * expPnegtime.t()
      + expPnegtime * covCurrent * DexpPnegtime.t();
    
      }


    
    cov22 = sPpP * cov22 / time;
    cov22.reshape(P.n_rows, P.n_cols);
    cov22 += cov22.t();

    cov23 = sPpP * cov23;
    cov23.reshape(P.n_rows, P.n_cols);


#ifdef RTS_SMOOTHER
    smoother.addPassData(false, timeEnd - time,
                         meanCurrent, meanHist,
                         covCurrent, covHist, covHC);
#endif

    meanHist = join_cols(meanHist, PinvIexpPntot * (meanCurrent - mu) + mu);

    meanCurrent = expPnegtime * (meanCurrent - mu) + mu;

    covHC = join_cols(covHC * expPnegtime.t(), PinvIexpPntot * covCurrent * expPnegtime.t() + cov23);
    
    covHist = join_cols(join_rows(covHist, cov12),
            join_rows(cov12.t(),
                  cov22
                  + PinvIexpPntot * covCurrent * PinvIexpPntot.t()
                  ));

    covCurrent = expPnegtime * covCurrent * expPnegtime.t() + cov33;

#ifdef RTS_SMOOTHER
    smoother.addPassData(true, timeEnd,
                         meanCurrent, meanHist,
                         covCurrent, covHist, covHC);
#endif

    smoother.addTime(PinvIexpPntot, expPnegtime);
  }

  void addInstantaneousObservation(const seriesPoint& observation)
  {
    sword series = observation.seriesIndex;
    double timeBegin = observation.timeEnd;
    double obs = observation.mean;
    double var = observation.variance;
    const sword seriesIndex = observation.seriesIndex;
    const sword houseIndex = observation.houseIndex;
    const uword k = observation.processIndex;
    
    if (timeIndices[timeIndices.size()-1] != timeBegin)
      {
    COUT << "index["<<timeIndices.size()-1<<"]: " <<
      timeIndices[timeIndices.size()-1] << endl;
    COUT << "time: " << timeBegin << endl;
      }

#ifndef MONOCAR_NO_DEBUG
    assert(timeIndices[timeIndices.size()-1] == timeBegin);
#endif

    double varPow = 0;
    if ((houseIndex >= 0) & (var > 0))
      varPow = params.getVarCenter(houseIndex, seriesIndex) *
    std::pow(var/params.getVarCenter(houseIndex, seriesIndex),
         params.getExtraVarPow(houseIndex, seriesIndex));

    vector<double> Dobs;
    vector<double> Dvar;

    for (uword i = 0; i < Dparams.size(); i++)
      {
    Dobs.push_back(0.0);
    Dvar.push_back(0.0);

    Dobs[i] -= Dparams[i].mean(seriesIndex, k < params.mean.n_cols ? k : 0);
    if (houseIndex >= 0)
      Dobs[i] -= Dparams[i].extraMean(houseIndex, k < params.extraMean.n_cols ? k : 0);

    if (params.hasExo)
      Dobs[i] -= observation.exoProduct(Dparams[i].exo);

    if (houseIndex >= 0)
      {
        if (var > 0)
          Dvar[i] = params.getExtraVarMult(houseIndex, seriesIndex)
        * Dparams[i].getExtraVarPow(houseIndex, seriesIndex)
        * log(var/params.getVarCenter(houseIndex, seriesIndex))
        * varPow
        + Dparams[i].getExtraVarMult(houseIndex, seriesIndex)
        * varPow
        + Dparams[i].getExtraVarAdd(houseIndex, seriesIndex);
        else
          Dvar[i] = Dparams[i].getExtraVarAdd(houseIndex, seriesIndex);          
      }
      }


    vector<double> DDobs;
    vector<double> DDvar;

    if (hasHess)
      {
    uword ij = 0;
    for (uword i = 0; i < Dparams.size(); i++)
      {
        for (uword j = i; j < Dparams.size(); j++)
          {
        DDobs.push_back(0.0);
        DDvar.push_back(0.0);
        
        DDobs[ij] -= DDparams[ij].mean(seriesIndex, k < params.mean.n_cols ? k : 0);
        if (houseIndex >= 0)
          DDobs[ij] -= DDparams[ij].extraMean(houseIndex, k < params.extraMean.n_cols ? k : 0);
        
        if (params.hasExo)
          DDobs[ij] -= observation.exoProduct(DDparams[ij].exo);
        
        if (houseIndex >= 0)
          {
            if (var > 0)
              DDvar[ij] =
            + Dparams[j].getExtraVarMult(houseIndex, seriesIndex)
            * Dparams[i].getExtraVarPow(houseIndex, seriesIndex)
            * log(var/params.getVarCenter(houseIndex, seriesIndex))
            * varPow
            
            + params.getExtraVarMult(houseIndex, seriesIndex)
            * DDparams[ij].getExtraVarPow(houseIndex, seriesIndex)
            * log(var/params.getVarCenter(houseIndex, seriesIndex))
            * varPow
            
            + params.getExtraVarMult(houseIndex, seriesIndex)
            * Dparams[i].getExtraVarPow(houseIndex, seriesIndex)
            * Dparams[j].getExtraVarPow(houseIndex, seriesIndex)
            * log(var/params.getVarCenter(houseIndex, seriesIndex))
            * log(var/params.getVarCenter(houseIndex, seriesIndex))
            * varPow
        
            + DDparams[ij].getExtraVarMult(houseIndex, seriesIndex)
            * varPow
            
            + Dparams[i].getExtraVarMult(houseIndex, seriesIndex)
            * Dparams[j].getExtraVarPow(houseIndex, seriesIndex)
            * log(var/params.getVarCenter(houseIndex, seriesIndex))
            * varPow
            
            + DDparams[ij].getExtraVarAdd(houseIndex, seriesIndex);
            else
              DDvar[ij] = DDparams[ij].getExtraVarAdd(houseIndex, seriesIndex);
            
            ij++;
          }
          }
      }
      }


    obs -= params.mean(seriesIndex, k < params.mean.n_cols ? k : 0);
    if (houseIndex >= 0)
      obs -= params.extraMean(houseIndex, k < params.extraMean.n_cols ? k : 0);

    if (params.hasExo)
      obs -= observation.exoProduct(params.exo);


    if (houseIndex >= 0)
      {
    var = params.getExtraVarMult(houseIndex, seriesIndex)
      * varPow
      + params.getExtraVarAdd(houseIndex, seriesIndex);
      }

    MyMatrix H(1, covCurrent.n_cols, fill::zeros);
    H(0,series) = 1.0;

    double y = obs - as_scalar(H * meanCurrent);
    double S = as_scalar(H * covCurrent * H.t());
    if (var > 0.0)
      S += var;

    if (S > 0.0)
      {    
    MyMatrix K = covCurrent * H.t() / S;
    logLik -= (y * y / S + log(S) + LOGTWOPI) / 2.0;
    MyMatrix HcHC = H * covHC.t();

    if (hasHess)
      {
        uword ij = 0;
        for (uword i = 0; i < Dparams.size(); i++)
          {
        double D1y = Dobs[i] - as_scalar(H * DmeanCurrent[i]);
        double D1S = as_scalar(H * DcovCurrent[i] * H.t()) + Dvar[i];
        double D1Sinv = - D1S / (S*S);
        MyMatrix D1K = DcovCurrent[i] * H.t() / S + covCurrent * H.t() * D1Sinv;
        MyMatrix D1HcHC = H * DcovHC[i].t();

        for (uword j = i; j < Dparams.size(); j++)
          {
            double D2y = Dobs[j] - as_scalar(H * DmeanCurrent[j]);
            double D2S = as_scalar(H * DcovCurrent[j] * H.t()) + Dvar[j];
            double D2Sinv = - D2S / (S*S);
            MyMatrix D2K = DcovCurrent[j] * H.t() / S + covCurrent * H.t() * D2Sinv;

            double DDy = DDobs[ij] - as_scalar(H * DDmeanCurrent[ij]);
            double DDS = as_scalar(H * DDcovCurrent[ij] * H.t()) + DDvar[ij];
            double DDSinv = 2 * D1S * D2S / (S*S*S) - DDS / (S*S);
            MyMatrix DDK =
              + DDcovCurrent[ij] * H.t() / S
              + DcovCurrent[i] * H.t() * D2Sinv
              + DcovCurrent[j] * H.t() * D1Sinv
              + covCurrent * H.t() * DDSinv;

            double DDlogLikChange =
              + DDy * y / S
              + D1y * D2y / S
              + D1y * y * D2Sinv
              + D2y * y * D1Sinv
              + (y * y * DDSinv) / 2.0
              + (DDS / S) / 2.0
              - (D1S * D2S / (S*S)) / 2.0;
            
            DDlogLik(i,j) -= DDlogLikChange;
            if (i != j)
              DDlogLik(j,i) -= DDlogLikChange;
          
            MyMatrix D2HcHC = H * DcovHC[j].t();
            MyMatrix DDHcHC = H * DDcovHC[ij].t();
        
            DDmeanCurrent[ij] += DDK * y + D1K * D2y + D2K * D1y + K * DDy;
          
            DDmeanHist[ij] +=
              + DDHcHC.t() * (y / S)
              + D1HcHC.t() * (D2y / S)
              + D1HcHC.t() * (y * D2Sinv)
            
              + D2HcHC.t() * (D1y / S)
              + HcHC.t() * (DDy / S)
              + HcHC.t() * (D1y * D2Sinv)
            
              + D2HcHC.t() * (y * D1Sinv)
              + HcHC.t() * (D2y * D1Sinv)
              + HcHC.t() * (y * DDSinv);
          
            DDcovCurrent[ij] -=
              + DDK * H * covCurrent
              + D1K * H * DcovCurrent[j]
              + D2K * H * DcovCurrent[i]
              + K * H * DDcovCurrent[ij];
          
            DDcovHist[ij] -=
              + DDHcHC.t() * HcHC / S
              + D1HcHC.t() * D2HcHC / S
              + D1HcHC.t() * HcHC * D2Sinv
            
              + D2HcHC.t() * D1HcHC / S
              + HcHC.t() * DDHcHC / S
              + HcHC.t() * D1HcHC * D2Sinv
            
              + D2HcHC.t() * HcHC * D1Sinv
              + HcHC.t() * D2HcHC * D1Sinv
              + HcHC.t() * HcHC * DDSinv;
          
            DDcovHC[ij] -= trans(DDK * HcHC + D1K * D2HcHC + D2K * D1HcHC + K * DDHcHC);
          
            ij++;
          }
          }
      }

    
    if (partialScores)
      DlogLik.resize(DlogLik.n_rows, DlogLik.n_cols + 1);
    for (uword i = 0; i < Dparams.size(); i++)
      {
        double Dy = Dobs[i] - as_scalar(H * DmeanCurrent[i]);
        double DS = as_scalar(H * DcovCurrent[i] * H.t()) + Dvar[i];
        double DSinv = - DS / (S*S);
        MyMatrix DK = DcovCurrent[i] * H.t() / S + covCurrent * H.t() * DSinv;

        if (partialScores)
          DlogLik(i,DlogLik.n_cols-1) = -(Dy * y / S + (y * y * DSinv + DS / S) / 2.0);
        else
          DlogLik(i,0) -= Dy * y / S + (y * y * DSinv + DS / S) / 2.0;

        MyMatrix DHcHC = H * DcovHC[i].t();
        
        DmeanCurrent[i] += DK * y + K * Dy;
        DmeanHist[i] += DHcHC.t() * (y / S) + HcHC.t() * (Dy / S) + HcHC.t() * (y * DSinv);
        DcovCurrent[i] -= DK * H * covCurrent + K * H * DcovCurrent[i];
        DcovHist[i] -= DHcHC.t() * HcHC / S + HcHC.t() * DHcHC / S + HcHC.t() * HcHC * DSinv;
        DcovHC[i] -= trans(DK * HcHC + K * DHcHC);
      }

#ifdef RTS_SMOOTHER
#ifdef DOUBLE_CHECK_SMOOTHER
          smoother.addPassData(false, timeIndices.back(),
                               meanCurrent, meanHist,
                               covCurrent, covHist, covHC);
#endif
#endif

    meanCurrent += K * y;
    meanHist += HcHC.t() * (y / S);
    covCurrent -= K * H * covCurrent;
    covHist -= HcHC.t() * HcHC / S;
    covHC -= trans(K * HcHC);

    smoother.addObservation(true,
                HcHC, K, H, y, S);

    /*
    HcHC = H * covHC.t();

    meanCurrent
    meanHist
    covCurrent
    covHist
    covHC

    KH=
    0
    0
    HcHC.t() * H / S
    K * H

    HtSiy=
    0
    H * y / S

    HtSiH=
    0
    0
    0
    H.t() * H / S    

    */

      }
    else
      {
    logLik -= y * y / 2.0;
    if (hasHess)
      {
        uword ij = 0;
        for (uword i = 0; i < Dparams.size(); i++)
          {
        double D1y = Dobs[i] - as_scalar(H * DmeanCurrent[i]);
        for (uword j = i; j < Dparams.size(); j++)
          {
            double D2y = Dobs[j] - as_scalar(H * DmeanCurrent[j]);
            double DDy = DDobs[ij] - as_scalar(H * DDmeanCurrent[ij]);
            DDlogLik(i,j) -= DDy * y + D1y * D2y;
            if (i != j)
              DDlogLik(j,i) -= DDy * y + D1y * D2y;
            ij++;
          }
          }
      }
    if (partialScores)
      DlogLik.resize(DlogLik.n_rows, DlogLik.n_cols + 1);
    for (uword i = 0; i < Dparams.size(); i++)
      {
        double Dy = Dobs[i] - as_scalar(H * DmeanCurrent[i]);
        
        if (partialScores)
          DlogLik(i,DlogLik.n_cols-1) = -Dy * y;
        else
          DlogLik(i,0) -= Dy * y;
      }

      }
    if (debugLevel > 5)
      COUT << "LL+: " << ((S>0.0) ? (- y * y / S - log(S) - LOGTWOPI) / 2.0 : (-y*y)/2.0) << "\ty: " << y << "\tobs[" << series << "," << timeBegin << "]=" << obs << "\tS: " << S << "\tvar: " << var << "\tmean: " << trans(meanCurrent) << endl;


  }

  void addPoint(const seriesPoint& observation)
  {
    removeTime(observation.cumTimeStart, true);
    if (observation.timeEnd > timeIndices.back())
      {
        removeTime(observation.cumTimeStart, false);
        addTime(observation.timeEnd);
      }


    // observation.print(logLik);
    if (observation.is_record())
      return record();
    else if (observation.is_instantaneous())
      return addInstantaneousObservation(observation);
    else if (observation.is_noninstantaneous())
      return addNoninstantaneousObservation(observation);
  }

  void addNoninstantaneousObservation(const seriesPoint& observation)
  {
    sword series = observation.seriesIndex;
    double timeBegin = observation.timeStart;
    double timeEnd = observation.timeEnd;
    double obs = observation.mean;
    double var = observation.variance;
    const sword seriesIndex = observation.seriesIndex;
    const sword houseIndex = observation.houseIndex;
    const uword k = observation.processIndex;
    
    double varPow = 0;
    if ((houseIndex >= 0) & (var > 0))
      varPow = params.getVarCenter(houseIndex, seriesIndex) *
    std::pow(var/params.getVarCenter(houseIndex, seriesIndex),
         params.getExtraVarPow(houseIndex, seriesIndex));


    vector<double> Dobs;
    vector<double> Dvar;

    for (uword i = 0; i < Dparams.size(); i++)
      {
    Dobs.push_back(0.0);
    Dvar.push_back(0.0);

    Dobs[i] -= Dparams[i].mean(seriesIndex, k < params.mean.n_cols ? k : 0);
    if (houseIndex >= 0)
      Dobs[i] -= Dparams[i].extraMean(houseIndex, k < params.extraMean.n_cols ? k : 0);

    if (params.hasExo)
      Dobs[i] -= observation.exoProduct(Dparams[i].exo);

    if (houseIndex >= 0)
      {
        if (var > 0)
          Dvar[i] = params.getExtraVarMult(houseIndex, seriesIndex)
        * Dparams[i].getExtraVarPow(houseIndex, seriesIndex)
        * log(var/Dparams[i].getVarCenter(houseIndex, seriesIndex))
        * varPow
        + Dparams[i].getExtraVarMult(houseIndex, seriesIndex)
        * varPow
        + Dparams[i].getExtraVarAdd(houseIndex, seriesIndex);
        else
          Dvar[i] = Dparams[i].getExtraVarAdd(houseIndex, seriesIndex);          
      }
      }


    vector<double> DDobs;
    vector<double> DDvar;

    if (hasHess)
      {
    uword ij = 0;
    for (uword i = 0; i < Dparams.size(); i++)
      {
        for (uword j = i; j < Dparams.size(); j++)
          {
        DDobs.push_back(0.0);
        DDvar.push_back(0.0);
        
        DDobs[ij] -= DDparams[ij].mean(seriesIndex, k < params.mean.n_cols ? k : 0);
        if (houseIndex >= 0)
          DDobs[ij] -= DDparams[ij].extraMean(houseIndex, k < params.extraMean.n_cols ? k : 0);

        if (params.hasExo)
          DDobs[ij] -= observation.exoProduct(DDparams[ij].exo);
        
        if (houseIndex >= 0)
          {
            if (var > 0)
              DDvar[ij] =
            + Dparams[j].getExtraVarMult(houseIndex, seriesIndex)
            * Dparams[i].getExtraVarPow(houseIndex, seriesIndex)
            * log(var/params.getVarCenter(houseIndex, seriesIndex))
            * varPow
            
            + params.getExtraVarMult(houseIndex, seriesIndex)
            * DDparams[ij].getExtraVarPow(houseIndex, seriesIndex)
            * log(var/params.getVarCenter(houseIndex, seriesIndex))
            * varPow
            
            + params.getExtraVarMult(houseIndex, seriesIndex)
            * Dparams[i].getExtraVarPow(houseIndex, seriesIndex)
            * Dparams[j].getExtraVarPow(houseIndex, seriesIndex)
            * log(var/params.getVarCenter(houseIndex, seriesIndex))
            * log(var/params.getVarCenter(houseIndex, seriesIndex))
            * varPow
        
            + DDparams[ij].getExtraVarMult(houseIndex, seriesIndex)
            * varPow
            
            + Dparams[i].getExtraVarMult(houseIndex, seriesIndex)
            * Dparams[j].getExtraVarPow(houseIndex, seriesIndex)
            * log(var/params.getVarCenter(houseIndex, seriesIndex))
            * varPow
            
            + DDparams[ij].getExtraVarAdd(houseIndex, seriesIndex);
            else
              DDvar[ij] = DDparams[ij].getExtraVarAdd(houseIndex, seriesIndex);

            ij++;
          }
          }
      }
      }

    obs -= params.mean(seriesIndex, k < params.mean.n_cols ? k : 0);
    if (houseIndex >= 0)
      obs -= params.extraMean(houseIndex, k < params.extraMean.n_cols ? k : 0);

    if (params.hasExo)
      obs -= observation.exoProduct(params.exo);


    if (houseIndex >= 0)
      {
    var = params.getExtraVarMult(houseIndex, seriesIndex)
      * varPow
      + params.getExtraVarAdd(houseIndex, seriesIndex);
      }


    MyMatrix H(1, covHist.n_cols, fill::zeros);
    int j = 0;
    double Hsum = 0.0;

    for (uword i = series; i < covHist.n_cols; i+=seriesSize)
      {
    if (timeIndices[j] > timeEnd)
      {
        break;
      }
    else if (timeIndices[j+1] >= timeBegin)
      {
        H(0,i) = max(0.0, min(timeEnd, timeIndices[j+1]) - max(timeBegin, timeIndices[j]));
        Hsum += H(0,i);
      }
    j++;
      }
    H /= Hsum;

    double y = obs - as_scalar(H * meanHist);
    double S = as_scalar(H * covHist * H.t());
    if (var > 0.0)
      S += var;
    MyMatrix K = covHist * H.t() / S;

    MyMatrix HcHC = H * covHC;

    if (hasHess)
      {
    uword ij = 0;
    for (uword i = 0; i < Dparams.size(); i++)
      {
        double D1y = Dobs[i] - as_scalar(H * DmeanHist[i]);
        double D1S = as_scalar(H * DcovHist[i] * H.t()) + Dvar[i];
        double D1Sinv = - D1S / (S*S);
        MyMatrix D1K = DcovHist[i] * H.t() / S + covHist * H.t() * D1Sinv;
        MyMatrix D1HcHC = H * DcovHC[i];
        for (uword j = i; j < Dparams.size(); j++)
          {        
        double D2y = Dobs[j] - as_scalar(H * DmeanHist[j]);
        double D2S = as_scalar(H * DcovHist[j] * H.t()) + Dvar[j];
        double D2Sinv = - D2S / (S*S);
        MyMatrix D2K = DcovHist[j] * H.t() / S + covHist * H.t() * D2Sinv;

        double DDy = DDobs[ij] - as_scalar(H * DDmeanHist[ij]);
        double DDS = as_scalar(H * DDcovHist[ij] * H.t()) + DDvar[ij];
        double DDSinv = 2 * D1S * D2S / (S*S*S) - DDS / (S*S);
        MyMatrix DDK =
          + DDcovHist[ij] * H.t() / S
          + DcovHist[i] * H.t() * D2Sinv
          + DcovHist[j] * H.t() * D1Sinv
          + covHist * H.t() * DDSinv;
        
        double DDlogLikChange =
          + DDy * y / S
          + D1y * D2y / S
          + D1y * y * D2Sinv
          + D2y * y * D1Sinv
          + (y * y * DDSinv) / 2.0
          + (DDS / S) / 2.0
          - (D1S * D2S / (S*S)) / 2.0;
        
        DDlogLik(i,j) -= DDlogLikChange;
        if (i != j)
          DDlogLik(j,i) -= DDlogLikChange;

        MyMatrix D2HcHC = H * DcovHC[j];
        MyMatrix DDHcHC = H * DDcovHC[ij];

        DDmeanHist[ij] += DDK * y + D1K * D2y + D2K * D1y + K * DDy;
 
        DDmeanCurrent[ij] +=
          + DDHcHC.t() * (y / S)
          + D1HcHC.t() * (D2y / S)
          + D1HcHC.t() * (y * D2Sinv)
          
          + D2HcHC.t() * (D1y / S)
          + HcHC.t() * (DDy / S)
          + HcHC.t() * (D1y * D2Sinv)
          
          + D2HcHC.t() * (y * D1Sinv)
          + HcHC.t() * (D2y * D1Sinv)
          + HcHC.t() * (y * DDSinv);

        DDcovHist[ij] -=
          + DDK * H * covHist
          + D1K * H * DcovHist[j]
          + D2K * H * DcovHist[i]
          + K * H * DDcovHist[ij];
        
        DDcovCurrent[ij] -=
          + DDHcHC.t() * HcHC / S
          + D1HcHC.t() * D2HcHC / S
          + D1HcHC.t() * HcHC * D2Sinv
          
          + D2HcHC.t() * D1HcHC / S
          + HcHC.t() * DDHcHC / S
          + HcHC.t() * D1HcHC * D2Sinv
          
          + D2HcHC.t() * HcHC * D1Sinv
          + HcHC.t() * D2HcHC * D1Sinv
          + HcHC.t() * HcHC * DDSinv;
        
        DDcovHC[ij] -= DDK * HcHC + D1K * D2HcHC + D2K * D1HcHC + K * DDHcHC;
        
        ij++;
          }
      }
      }

    

    if (partialScores)
      DlogLik.resize(DlogLik.n_rows, DlogLik.n_cols + 1);

    for (uword i = 0; i < Dparams.size(); i++)
      {
    double Dy = Dobs[i] - as_scalar(H * DmeanHist[i]);
    double DS = as_scalar(H * DcovHist[i] * H.t()) + Dvar[i];
    double DSinv = - DS / (S*S);
    MyMatrix DK = DcovHist[i] * H.t() / S + covHist * H.t() * DSinv;

    if (partialScores)
      DlogLik(i,DlogLik.n_cols-1) = -(Dy * y / S + (y * y * DSinv + DS / S) / 2.0);
    else
      DlogLik(i,0) -= Dy * y / S + (y * y * DSinv + DS / S) / 2.0;


    MyMatrix DHcHC = H * DcovHC[i];


    DmeanHist[i] += DK * y + K * Dy;
    DmeanCurrent[i] += DHcHC.t() * (y / S) + HcHC.t() * (Dy / S) + HcHC.t() * (y * DSinv);
    DcovHist[i] -= DK * H * covHist + K * H * DcovHist[i];
    DcovCurrent[i] -= DHcHC.t() * HcHC / S + HcHC.t() * DHcHC / S + HcHC.t() * HcHC * DSinv;
    DcovHC[i] -= DK * HcHC + K * DHcHC;

      }

#ifdef RTS_SMOOTHER
#ifdef DOUBLE_CHECK_SMOOTHER
    smoother.addPassData(false, timeIndices.back(),
                         meanCurrent, meanHist,
                         covCurrent, covHist, covHC);
#endif
#endif

    logLik -= (y * y / S + log(S) + LOGTWOPI) / 2.0;

    meanHist += K * y;
    meanCurrent += HcHC.t() * (y / S);
    covHist -= K * H * covHist;
    covCurrent -= HcHC.t() * HcHC / S;
    covHC -= K * HcHC;

    smoother.addObservation(false,
                HcHC, K, H, y, S);

    /*
    HcHC = H * covHC;

    meanCurrent
    meanHist
    covCurrent
    covHist
    covHC

    KH=
    K * H
    HcHC.t() * H / S
    0
        0

    HtSiy=
    H * y / S
    0

    HtSiH=
    H.t() * H / S    
    0
    0
    0

    */


  }

  MyMatrix getMean(void) const
  {
    return(meanHist);
  }

  MyMatrix getCov(void) const
  {
    return(covHist);
  }

  double getLogLik(void) const
  {
    return(logLik);
  }

  MyMatrix& getDLogLik(void)
  {
    return(DlogLik);
  }

  MyMatrix& getDDLogLik(void)
  {
    return(DDlogLik);
  }

  
  void removeTime(const double timeStart, const bool allowEmpty)
  {
    uword removeSize = 0;
    vector<double>::iterator removeEnd = timeIndices.begin();
    //while( (*(removeEnd) < timeStart) && (removeEnd+2 != timeIndices.end()) && (removeEnd+1 != timeIndices.end()) && (removeEnd != timeIndices.end()) )
    while((removeEnd != timeIndices.end()) && (*(removeEnd) < timeStart))
      {
        removeEnd++;
        removeSize++;
      }

    if ((!allowEmpty) &&
        (removeSize >= meanHist.n_rows))
      {
        removeSize = (meanHist.n_rows > 0) ? (meanHist.n_rows - 1) : 0;
        removeEnd = timeIndices.begin() += removeSize;
      }

    if (removeSize == 0)
        return;

    if (removeEnd == timeIndices.end())
      {
        removeSize--;
        removeEnd--;
      }

    if (removeSize == 0)
        return;

    timeIndices.erase(timeIndices.begin(), removeEnd);

    removeSize *= seriesSize;

    if (meanHist.n_rows < removeSize)
      {
        throw std::runtime_error("Bad removal");
      }

    for (uword i = 0; i < Dparams.size(); i++)
      {
    trimMatrix(DmeanHist[i], removeSize, 0);
    trimMatrix(DcovHist[i], removeSize, removeSize);
    trimMatrix(DcovHC[i], removeSize, 0);
      }

    if (hasHess)
      {
    uword ij = 0;
    for (uword i = 0; i < Dparams.size(); i++)
      {
        for (uword j = i; j < Dparams.size(); j++)
          {
        trimMatrix(DDmeanHist[ij], removeSize, 0);
        trimMatrix(DDcovHist[ij], removeSize, removeSize);
        trimMatrix(DDcovHC[ij], removeSize, 0);
        ij++;
          }
      }
      }

#ifdef RTS_SMOOTHER
    smoother.addPassData(false, timeStart,
                         meanCurrent, meanHist,
                         covCurrent, covHist, covHC);
#endif

    trimMatrix(meanHist, removeSize, 0);
    trimMatrix(covHist, removeSize, removeSize);
    trimMatrix(covHC, removeSize, 0);

#ifdef RTS_SMOOTHER
    smoother.addPassData(true, timeStart,
                         meanCurrent, meanHist,
                         covCurrent, covHist, covHC);
#endif
    
    smoother.removeTime(removeSize);
  }

private:
  bool partialScores;

  double logLik;

  size_t seriesSize;
  vector<double> timeIndices;

  const OUparam& params;

  const vector<OUparam>& Dparams;

  const vector<OUparam>& DDparams;

  MyMatrix mu;
  MyMatrix P;
  MyMatrix Sigma;

  MyMatrix Ident;
  MyMatrix PpP;
  MyMatrix sPpP;
  MyMatrix sP;

  MatrixExp expP;

  MyMatrix covHist;
  MyMatrix covHC;
  MyMatrix covCurrent;

  MyMatrix meanHist;
  MyMatrix meanCurrent;

  MyMatrix DlogLik;

  MatrixExpFrechet DexpP;

  vector<MyMatrix> DPpP;
  vector<MyMatrix> DsPpP;
  vector<MyMatrix> DsP;
  vector<MyMatrix> Dmu;

  vector<MyMatrix> DmeanHist;
  vector<MyMatrix> DmeanCurrent;
  vector<MyMatrix> DcovHist;
  vector<MyMatrix> DcovHC;
  vector<MyMatrix> DcovCurrent;

  MyMatrix DDlogLik;

  MatrixExpDoubleFrechet DDexpP;

  vector<MyMatrix> DDPpP;
  vector<MyMatrix> DDsPpP;
  vector<MyMatrix> DDsP;
  vector<MyMatrix> DDmu;

  vector<MyMatrix> DDmeanHist;
  vector<MyMatrix> DDmeanCurrent;
  vector<MyMatrix> DDcovHist;
  vector<MyMatrix> DDcovHC;
  vector<MyMatrix> DDcovCurrent;
  
  double lastTime;
  MyMatrix expPnegtime;
  MyMatrix IexpPpPnegtime;

  MyMatrix DexpPnegtime;
  MyMatrix DIexpPpPnegtime;

  MyMatrix DDexpPnegtime;
  MyMatrix DDIexpPpPnegtime;

  modifiedBrysonFrazier smoother;

  bool hasHess;
};


class OUsimulation
{
public:
  OUsimulation(const MyMatrix iParams,
           const inputVector& iTimeIndices,
           // const inputIntegerVector& iProcessIndices,
           const size_t iProcessIndexSize,
           const size_t iHouseSize,
           // const inputIntegerVector& iSeriesIndices,
           // const inputIntegerVector& iHouseIndices,
           // const inputVector& iTimeStart,
           // const inputVector& iTimeEnd,
           // const MyMatrix& iExoData,
           // const inputVector& iVariances,
           // const inputVector& iVarCenters,
           const bool iHasStart,
           const bool iHasExo,
           const bool iVarByHouse,
           const inputVector& iVarCenters,
           const inputMatrix& iStartConstraint,
           const inputMatrix& iMeanConstraint,
           const inputMatrix& iExtraMeanConstraint,
           const inputMatrix& iReversionConstraint,
           const inputMatrix& iSimultaneousConstraint,
           const inputMatrix& iSigmaSqrtConstraint,
           const inputMatrix& iExtraVarAddConstraint,
           const inputMatrix& iExtraVarMultConstraint,
           const inputMatrix& iExtraVarPowConstraint,
           const inputMatrix& iExoConstraint)
    :timeIndices(iTimeIndices)
    // ,processIndices(iProcessIndices)
    // ,processIndexSize(iProcessIndexSize)
    // ,seriesIndices(iSeriesIndices)
    // ,houseIndices(iHouseIndices)
    // ,timeStart(iTimeStart)
    // ,timeEnd(iTimeEnd)
    // ,exoData(iExoData)
    // ,variances(iVariances)
    // ,varCenters(iVarCenters)
    ,paramSpace(iProcessIndexSize,
        iHouseSize,
        iHasStart,
        iHasExo,
        iVarByHouse,
        iVarCenters,
        iStartConstraint,
        iMeanConstraint,
        iExtraMeanConstraint,
        iReversionConstraint,
        iSimultaneousConstraint,
        iSigmaSqrtConstraint,
        iExtraVarAddConstraint,
        iExtraVarMultConstraint,
        iExtraVarPowConstraint,
        iExoConstraint)
    ,sampleCurrent(iMeanConstraint.INPUTNROW,timeIndices.size(),iProcessIndexSize,fill::zeros)
    ,sampleHist(iMeanConstraint.INPUTNROW,timeIndices.size()-1,iProcessIndexSize,fill::zeros)
  {
    bool fail;
    params = paramSpace.getParams(iParams, fail);
    vector<OUparam> derivParams;
    vector<OUparam> doubleDerivParams;
    for (uword k = 0; k < ((uword) iProcessIndexSize); k++)
      {
    OUprocess series(params.process(k),
             derivParams,
             doubleDerivParams,
//             timeIndices(0),
             false,
             false);
    series.reset(timeIndices(0));
    sampleCurrent.slice(k).col(0) = series.sampleCurrent();
    if (params.mean.n_cols > 1)
      sampleCurrent.slice(k).col(0) += params.mean.col(k);
    else
      sampleCurrent.slice(k).col(0) += params.mean;
    for (int i = 1; i < timeIndices.size(); i++)
      {
        series.addTime(timeIndices(i));
        sampleCurrent.slice(k).col(i) = series.sampleCurrent();
        if (timeIndices(i)-timeIndices(i-1) > 1e-10)
          sampleHist.slice(k).col(i-1) = series.sampleHist();
        else
          sampleHist.slice(k).col(i-1) = sampleCurrent.slice(k).col(i);
        if (params.mean.n_cols > 1)
          {
        sampleCurrent.slice(k).col(i) += params.mean.col(k);
        sampleHist.slice(k).col(i-1) += params.mean.col(k);
          }
        else
          {
        sampleCurrent.slice(k).col(i) += params.mean;
        sampleHist.slice(k).col(i-1) += params.mean;
          }
        series.removeTime(timeIndices(i), true);
      }
      }
        
  }

  inputVector getSampleData(const inputIntegerVector& processIndices,
             const inputIntegerVector& seriesIndices,
             const inputIntegerVector& houseIndices,
             const inputVector& timeStart,
             const inputVector& timeEnd,
             const MyMatrix& exoData,
             const inputVector& variances)
  {
    inputVector data(timeStart.size());

    for (int i = 0; i < timeStart.size(); i++)
      {
    double obs = 0;
    double obsLength = 0;
    for (int j = 0; j < timeIndices.size(); j++)
      {
        if (timeEnd(i) < timeIndices(j))
          {
        break;
          }
        else
          {
        if (timeStart(i) == timeIndices(j))
          {
            obs = sampleCurrent(seriesIndices(i),j,processIndices(i));
            break;
          }
        else if ((timeStart(i) < timeIndices(j)) && (j > 0))
          {
            obs += sampleHist(seriesIndices(i),j-1,processIndices(i))
              * (timeIndices(j) - timeIndices(j-1));
            obsLength += (timeIndices(j) - timeIndices(j-1));
          }
          }
      }
    if (obsLength > 0)
      obs /= obsLength;

    if (houseIndices(i) >= 0)
      obs += params.extraMean(houseIndices(i),
                  (processIndices(i) < ((int) params.extraMean.n_cols)) ? processIndices(i) : 0);
    
    if (params.hasExo)
      obs += as_scalar(exoData.row(i) * params.exo);

    double var = variances(i);
        if ((houseIndices(i) >= 0) & (var > 0))
      {
        double varPow = params.getVarCenter(houseIndices(i),seriesIndices(i)) *
          std::pow(var/params.getVarCenter(houseIndices(i),seriesIndices(i)),
               params.getExtraVarPow(houseIndices(i), seriesIndices(i)));
        var = params.getExtraVarMult(houseIndices(i), seriesIndices(i))
          * varPow
          + params.getExtraVarAdd(houseIndices(i), seriesIndices(i));
      }
    if (var>0)
      data(i) = R::rnorm(obs, std::sqrt(var));
    else
      data(i) = obs;
      }
    return data;
  }

  MyCube& getSampleCurrent()
  {
    return sampleCurrent;
  }
  
  MyCube& getSampleHist()
  {
    return sampleHist;
  }
  
  // MyMatrix getSampleCurrentObs()
  // {
  //   return sampleCurrentObs;
  // }
  
  // MyMatrix getSampleHistObs()
  // {
  //   return sampleHistObs;
  // }
  
private:
  const inputVector& timeIndices;
  // const inputIntegerVector& processIndices;
  // const size_t processIndexSize;
  // const inputIntegerVector& seriesIndices;
  // const inputIntegerVector& houseIndices;
  // const inputVector& timeStart;
  // const inputVector& timeEnd;
  // const MyMatrix& exoData;
  // const inputVector& variances;
  // const inputVector& varCenters;
  const OUparamspace paramSpace;
  MyCube sampleCurrent;
  MyCube sampleHist;
  // MyMatrix sampleCurrentObs;
  // MyMatrix sampleHistObs;
  OUparam params;
};

class OUsuper
{
public:
  OUsuper(const inputVector& iTimeIndices,
      const inputIntegerVector& iProcessIndices,
      const size_t iProcessIndexSize,
      const size_t iHouseSize,
      const inputIntegerVector& iSeriesIndices,
      const inputIntegerVector& iHouseIndices,
      const inputVector& iTimeStart,
      const inputVector& iTimeEnd,
      const MyMatrix& iExoData,
      const inputVector& iObservations,
      const inputVector& iVariances,
      const inputVector& iVarCenters,
      const bool iHasStart,
      const bool iHasExo,
      const bool iVarByHouse,
      const inputMatrix& iStartConstraint,
      const inputMatrix& iMeanConstraint,
      const inputMatrix& iExtraMeanConstraint,
      const inputMatrix& iReversionConstraint,
      const inputMatrix& iSimultaneousConstraint,
      const inputMatrix& iSigmaSqrtConstraint,
      const inputMatrix& iExtraVarAddConstraint,
      const inputMatrix& iExtraVarMultConstraint,
      const inputMatrix& iExtraVarPowConstraint,
      const inputMatrix& iExoConstraint,
	  const inputMatrix& iPriorMatrix,
          const double iEpsilon = 0,
          const size_t iGrainSize = 0,
          const vector<double>& iHistTimes = vector<double>())
    :mEpsilon(iEpsilon)
    ,processIndexSize(iProcessIndexSize)
//    ,timeIndices(iTimeIndices)
//    ,processIndices(iProcessIndices)
//    ,seriesIndices(iSeriesIndices)
//    ,houseIndices(iHouseIndices)
//    ,timeStart(iTimeStart)
//    ,timeCumStart(iTimeCumStart)
//    ,timeEnd(iTimeEnd)
//    ,exoData(iExoData)
//    ,observations(iObservations)
//    ,variances(iVariances)
    ,paramSpace(iProcessIndexSize,
        iHouseSize,
        iHasStart,
        iHasExo,
        iVarByHouse,
        iVarCenters,
        iStartConstraint,
        iMeanConstraint,
        iExtraMeanConstraint,
        iReversionConstraint,
        iSimultaneousConstraint,
        iSigmaSqrtConstraint,
        iExtraVarAddConstraint,
        iExtraVarMultConstraint,
        iExtraVarPowConstraint,
        iExoConstraint)
    ,priors(OUprior::getPriorFromMatrix(iPriorMatrix))
    ,grainSize(iGrainSize)
    ,points(processIndexSize)
    ,hasHist(iHistTimes.size() > 0)
  {
    for (size_t k = 0; k < processIndexSize; k++)
      {
        vector<seriesPoint>& processPoints = points[k];
        set<double> splitTimes;
        for (INPUTINT i = 0; i < iTimeStart.size(); i++)
          {
            if (iProcessIndices(i) == k)
              {
                seriesPoint newObs(k,
                                   iSeriesIndices(i),
                                   iHouseIndices(i),
                                   iTimeStart(i),
                                   iTimeEnd(i),
                                   iObservations(i),
                                   iVariances(i),
                                   iExoData.row(i));
                splitTimes.insert(newObs.timeStart);
                processPoints.push_back(newObs);
              }
          }
        for (size_t i = 0; i < iHistTimes.size(); ++i)
          {
            seriesPoint newObs(k, iHistTimes[i], true);
            processPoints.push_back(newObs);
          }
        for (vector<seriesPoint>::iterator ptx = processPoints.begin();
             ptx != processPoints.end();
             ++ptx)
          {
            splitTimes.erase(ptx->timeEnd);
          }
        for (set<double>::const_iterator ix = splitTimes.begin(); ix != splitTimes.end(); ++ix)
          {
            seriesPoint newObs(k, *ix, false);
            processPoints.push_back(newObs);
          }
        std::sort(processPoints.begin(), processPoints.end());
        for (vector<seriesPoint>::iterator ptx = processPoints.begin();
             ptx != processPoints.end();
             ++ptx)
          {
            ptx->cumTimeStart = std::min_element(ptx, processPoints.end(), startTimeLessThan)->timeStart;
          }
      }
  }

  OUsuper(const OUsuper& ou)
    :mEpsilon(ou.mEpsilon)
    ,processIndexSize(ou.processIndexSize)
//    ,timeIndices(ou.timeIndices)
//    ,processIndices(ou.processIndices)
//    ,seriesIndices(ou.seriesIndices)
//    ,houseIndices(ou.houseIndices)
//    ,timeStart(ou.timeStart)
//    ,timeCumStart(ou.timeCumStart)
//    ,timeEnd(ou.timeEnd)
//    ,exoData(ou.exoData)
//    ,observations(ou.observations)
//    ,variances(ou.variances)
    ,paramSpace(ou.paramSpace)
    ,grainSize(ou.grainSize)
    ,points(ou.points)
    ,hasHist(ou.hasHist)
  {
  }

  double getEpsilon(void)
  {
    return mEpsilon;
  }

  void printParams(const MyMatrix& params)
  {
    bool fail;
    OUparam newParams = paramSpace.getParams(params, fail);
    if (fail)
      {
    COUT << endl << "Failed to set parameters." << endl << endl;
    return;
      }
    COUT << endl;
    COUT << (newParams.hasStart ? "Start" : "No start") << endl;
    COUT << newParams.start.t() << endl;
    COUT << newParams.mean.t() << endl;
    COUT << newParams.extraMean.t() << endl;
    COUT << newParams.reversion << endl;
    COUT << newParams.sigma << endl;
    COUT << newParams.extraVarAdd.t() << endl;
    COUT << newParams.extraVarMult.t() << endl;
    COUT << newParams.extraVarPow.t() << endl;
    COUT << (newParams.hasExo ? "Exo" : "No exo") << endl;
    COUT << newParams.exo.t() << endl;
    COUT << endl;
  }
  
  void estimate(double& logLik, vector<MyMatrix>& score, vector<MyMatrix>& hess,
          const MyMatrix& params, const uword i = 0, const uword j = 0,
          const bool partialScores = false, const bool getHess = false) const
  {
    score.clear();
    hess.clear();

#ifndef MONOCAR_NO_DEBUG
    assert(j <= params.n_rows);
#endif

    logLik = 0.0;
    
    try
      {

    bool showOutput = false;
    if (partialScores && (debugLevel >= 1) && (debugLevel <= 1))
      showOutput = true;
    if (getHess && (debugLevel >= 0) && (debugLevel <= 1))
      showOutput = true;
    long totalIterations = 0;
    for (size_t ptnum = 0; ptnum < points.size(); ++ptnum)
      totalIterations += points[ptnum].size();
    if (totalIterations <= 0)
      showOutput = false;
    long usedIterations = 0;
    long usedScreen = 0;
    if (showOutput)
      {
        COUT << endl << endl;
        COUT << "Computing";
        if (partialScores && (debugLevel >= 1) && (debugLevel <= 1))
          {
        COUT << " partial scores";
        if (getHess)
          COUT << " and";
          }
        if (getHess && (debugLevel >= 0) && (debugLevel <= 1))
          COUT << " Hessian";
        COUT << "..." << endl;
        for (int kx = 0; kx < screenWidth; kx++)
          COUT << "_" << flush;
        COUT << endl;
      }

        bool fail = false;
        OUparam newParams = paramSpace.getParams(params, fail);

        if (fail)
          {
            logLik = -(numeric_limits<double>::infinity());
            return;
          }

        vector<OUparam> derivParams;
        for (uword iDP = i; iDP < j; iDP++)
            derivParams.push_back(paramSpace.getDParams(params, iDP));
        vector<OUparam> derivParamsSeries(derivParams.size());

        vector<OUparam> doubleDerivParams;
        if (getHess)
            for (uword iD1P = i; iD1P < j; iD1P++)
                for (uword iD2P = iD1P; iD2P < j; iD2P++)
                    doubleDerivParams.push_back(paramSpace.getDDParams(params, iD1P, iD2P));
        vector<OUparam> doubleDerivParamsSeries(doubleDerivParams.size());

        if (det(newParams.reversion) == 0.0)
          {
            logLik = -(numeric_limits<double>::infinity());
            return;
          }

    for (size_t k = 0; k < points.size(); k++)
      {
        const vector<seriesPoint>& processPoints(points[k]);

        OUparam newParamsSeries = newParams.process(k);
        for (size_t m = 0; m < derivParams.size(); ++m)
            derivParamsSeries[m] = derivParams[m].process(k);
        for (size_t m = 0; m < doubleDerivParams.size(); ++m)
            doubleDerivParamsSeries[m] = doubleDerivParams[m].process(k);

        OUprocess series(newParamsSeries,
                         derivParamsSeries,
                         doubleDerivParamsSeries,
                         hasHist,
                         partialScores);

        double currentTime = -(numeric_limits<double>::infinity());
        if (processPoints.size() > 0)
            currentTime = processPoints.front().cumTimeStart;
        series.reset(currentTime);

        for (vector<seriesPoint>::const_iterator ptx = processPoints.begin();
             ptx != processPoints.end();
             ++ptx)
          {
            if (showOutput)
              {
                usedIterations++;
                while (usedScreen * totalIterations < long(screenWidth) * usedIterations)
                  {
                    usedScreen++;
                    COUT << "=" << flush;
                  }
                CHECK_USER_INTERRUPT;
              }
            series.removeTime(ptx->cumTimeStart, false);
            if (ptx->timeEnd > currentTime)
              {
                series.addTime(ptx->timeEnd);
                currentTime = ptx->timeEnd;
              }
            series.removeTime(ptx->cumTimeStart, true);
            series.addPoint(*ptx);
          }
        
        logLik += series.getLogLik();

        if (hasHist)
          {
            vector<double> inputTime;
            MyCube inputMean, inputCov;
            series.smooth(inputTime, inputMean, inputCov);

            MyMatrix outputTime, outputMean, outputStdDev;

            outputTime.set_size(inputTime.size(), 1);
            outputMean.set_size(inputMean.n_slices, inputMean.n_rows);
            outputStdDev.set_size(inputCov.n_slices, inputCov.n_rows);
            for (size_t i = 0; i < inputTime.size(); ++i)
              {
                outputTime[i] = inputTime[i];
                outputMean.row(i) = inputMean.slice(i).t();
                outputStdDev.row(i) = arma::sqrt(inputCov.slice(i).diag().t());
              }

            score.push_back(outputMean);
            hess.push_back(outputStdDev);
          }
	if (j > i)
          {
            score.push_back(series.getDLogLik());
            if (getHess)
                hess.push_back(series.getDDLogLik());
          }
      }
    for (std::vector<OUprior>::const_iterator prior_ptr = priors.begin();
	 prior_ptr != priors.end();
	 ++prior_ptr)
      {
	logLik += prior_ptr->getLogProb(newParams);
	if (j > i)
	  {
	    MyMatrix prior_score(derivParams.size(), 1);
	    for (size_t i = 0; i < derivParams.size(); ++i)
	      {
		prior_score[i] = prior_ptr->getDLogProb(newParams, derivParams[i]);
	      }
	    score.push_back(prior_score);
	  }
      }
    if (showOutput)
      {
        usedIterations++;
        while (usedScreen < long(screenWidth))
          {
	    usedScreen++;
	    COUT << "=" << flush;
          }
        COUT << endl;
        CHECK_USER_INTERRUPT;
      }
    return;
      }
  catch (const std::runtime_error& error)
    {
      if (debugLevel > 1)
	COUT << "Runtime error in likelihood!" << endl;
      else if (debugLevel >= 0)
	printDot("!");
      logLik = -(numeric_limits<double>::infinity());
      return;
    }

  }

    void getHessian(double& logLik, MyMatrix& score, MyMatrix& hess,
                    const MyMatrix& params, const uword i = 0, const uword j = 0) const
    {
      vector<MyMatrix> scoreVector, hessVector;
      estimate(logLik, scoreVector, hessVector, params, i, j, false, true);
      score.zeros(j-i, 1);
      for (vector<MyMatrix>::const_iterator ix = scoreVector.begin(); ix != scoreVector.end(); ++ix)
          score += *ix;
      hess.zeros(j-i, j-i);
      for (vector<MyMatrix>::const_iterator ix = hessVector.begin(); ix != hessVector.end(); ++ix)
          hess += *ix;
    }

    void getScore(double& logLik, MyMatrix& score,
                  const MyMatrix& params,
                  const uword i = 0, const uword j = 0,
                  const bool partialScores = false) const
    {
      vector<MyMatrix> scoreVector, hessVector;
      estimate(logLik, scoreVector, hessVector, params, i, j, partialScores, false);
      if (partialScores)
        {
          score.zeros(j-i, 0);
          for (vector<MyMatrix>::const_iterator ix = scoreVector.begin(); ix != scoreVector.end(); ++ix)
              score = join_rows(score, *ix);
        }
      else
        {
          score.zeros(j-i, 1);
          for (vector<MyMatrix>::const_iterator ix = scoreVector.begin(); ix != scoreVector.end(); ++ix)
              score += *ix;
        }
    }
    
    double getLikelihood(const MyMatrix& params) const
    {
      double logLik;
      vector<MyMatrix> scoreVector, hessVector;
      estimate(logLik, scoreVector, hessVector, params);
      return logLik;
    }

    double getLikelihood(const MyMatrix& params, const double epsi, const uword i, const double epsj = 0.0, const uword j = 0) const
    {
      double logLik;
      vector<MyMatrix> scoreVector, hessVector;

      MyMatrix newParams(params);

      newParams(i, 0) += epsi;
      newParams(j, 0) += epsj;

      estimate(logLik, scoreVector, hessVector, newParams);
      return logLik;
    }
  
  void getHist(const MyMatrix& params,
          vector<vector<double> >& outputTime,
          vector<MyMatrix>& outputMean,
          vector<MyMatrix>& outputVar)
  {
    double logLik;
    estimate(logLik, outputMean, outputVar, params, 0, 0, false, false);

    outputTime.clear();
    for (size_t k = 0; k < points.size(); ++k)
      {
        vector<double> seriesTime;
        for (size_t i = 0; i < points[k].size(); ++i)
            if (points[k][i].is_record())
                seriesTime.push_back(points[k][i].timeEnd);
        outputTime.push_back(seriesTime);
      }

    return;
  }


  const size_t getGrainSize(void) const
  {
    return grainSize;
  }

private:
  const double mEpsilon;

  const size_t processIndexSize;
//  const inputVector timeIndices;
//  const inputIntegerVector processIndices;
//  const inputIntegerVector seriesIndices;
//  const inputIntegerVector houseIndices;
//  const inputVector timeStart;
//  const inputVector timeCumStart;
//  const inputVector timeEnd;
//
//  const MyMatrix exoData;
//
//  const inputVector observations;
//  const inputVector variances;

  const OUparamspace paramSpace;

  const std::vector<OUprior> priors;

  const size_t grainSize;

  vector<vector<seriesPoint> > points;
  const bool hasHist;
};


#ifdef USE_NLOPT
extern "C" {
  double numderiv(OUsuper* ou, MyMatrix& param, double eps, uword i)
  {
    double logLikPlus = -ou->getLikelihood(param, eps, i);
    double logLikMinus = -ou->getLikelihood(param, -eps, i);
    return (logLikPlus - logLikMinus) / (2 * eps);
  }

  double numderivside(OUsuper* ou, MyMatrix& param, double eps, uword i)
  {
    double logLikZero = -ou->getLikelihood(param);
    double logLikPlus = -ou->getLikelihood(param, eps, i);
    return (logLikPlus - logLikZero) / (eps);
  }

  double numsecderiv(OUsuper* ou, MyMatrix param, double eps, uword i, uword j)
  {
    double llpp = -ou->getLikelihood(param, eps, i, eps, j);
    double llpm = -ou->getLikelihood(param, eps, i, -eps, j);
    double llmm = -ou->getLikelihood(param, -eps, i, -eps, j);
    double llmp = -ou->getLikelihood(param, -eps, i, eps, j);
    return ((llpp - llpm)/(eps+eps) + (llmm - llmp)/(eps+eps))/(eps+eps);
  }

  double numsecderivdiag(OUsuper* ou, MyMatrix param, double eps, uword i)
  {
    double llo = -ou->getLikelihood(param);
    double llpp = -ou->getLikelihood(param, eps+eps, i);
    double llmm = -ou->getLikelihood(param, -(eps+eps), i);
    return ((llpp - llo)/(eps+eps) + (llmm - llo)/(eps+eps))/(eps+eps);
  }

  static double oufun(uword n, const double *p, void *ex)
  {
    CHECK_USER_INTERRUPT;

    MyMatrix param(n, 1);
    for (uword i = 0; i < n; i++)
      param(i, 0) = p[i];
    
    OUsuper* ou = (OUsuper*) ex;
    double logLik = ou->getLikelihood(param);
    if (debugLevel > 4)
      ou->printParams(param);
    if (debugLevel > 3)
      COUT << "Log-likelihood: " << logLik << endl;
    return(-logLik);
  }

  static void ougrad_analytic(uword n, const double *p, double *gr, void *ex)
  {
    MyMatrix param(n, 1);
    for (uword i = 0; i < n; i++)
      param(i, 0) = p[i];
    
    OUsuper* ou = (OUsuper*) ex;

    if (debugLevel > 2)
      ou->printParams(param);

    double logLik;
    MyMatrix score;

    ou->getScore(logLik, score, param, 0, n);    
#ifndef MONOCAR_NO_DEBUG
    assert(n == score.n_rows);
#endif
    for (uword i = 0; i < n; i++)
      gr[i] = -score(i, 0);
  }

  static void ougrad_numeric(uword n, const double *p, double *gr, void *ex)
  {
    MyMatrix param(n, 1);
    for (uword i = 0; i < n; i++)
      param(i, 0) = p[i];
    
    OUsuper* ou = (OUsuper*) ex;

    const double eps = ou->getEpsilon();

    if (debugLevel > 2)
      ou->printParams(param);
      
    if (debugLevel > 1)
      COUT << "Grad (n):" << flush;

    CHECK_USER_INTERRUPT;
    for (uword i = 0; i < n; i++)
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

    if (debugLevel > 1)
      COUT << endl;
  }

  static void ougrad(uword n, const double *p, double *gr, void *ex)
    {
      OUsuper* ou = (OUsuper*) ex;
      const double eps = ou->getEpsilon();

      if (eps <= 0.0)
    ougrad_analytic(n, p, gr, ex);
      else
    ougrad_numeric(n, p, gr, ex);

    if (debugLevel > 1)
      {
    COUT << "Gradient:" << flush;
    for (uword i = 0; i < n; i++)
      COUT << " " << gr[i] << flush;
    COUT << endl;
      }
    }
  
  static void ougradloud(uword n, const double *p, double *gr, void *ex)
    {
      if (debugLevel == 1)
    printDot(":");
      if (debugLevel == 0)
    printDot(".");
      ougrad(n, p, gr, ex);
    }

  static void ouhess_analytic(uword n, const double *p,
             double *hess,
             void *ex)
  {
    MyMatrix param(n, 1);
    for (uword i = 0; i < n; i++)
      param(i, 0) = p[i];
    
    OUsuper* ou = (OUsuper*) ex;

    if (debugLevel > 2)
      ou->printParams(param);

    if (debugLevel > 1) {
      COUT << "Computing Hessian" << endl;
      // for (uword i = 0; i < n; i++)
      //     for (uword j = i; j < n; j++)
      //       COUT << "_";
      // COUT << endl;
    }

    double logLik;
    MyMatrix score;
    MyMatrix hessianMatrix;

    ou->getHessian(logLik, score, hessianMatrix, param, 0, n, false, true);
#ifndef MONOCAR_NO_DEBUG
    assert(n == score.n_rows);
#endif
    // for (uword i = 0; i < n; i++)
    //   gr[i] = -score(i, 0);
    uword ij = 0;
    for (uword i = 0; i < n; i++)
      {
    for (uword j = 0; j < n; j++)
      {
        hess[ij] = -hessianMatrix(i, j);
        ij++;
      }
      }

    CHECK_USER_INTERRUPT;
    if (debugLevel > 1)
      printDot(".");
    if (debugLevel == 1 || debugLevel == 0)
      printDot("*");
  }

  static void ouhess_numeric(uword n, const double *p,
             double *hess,
             void *ex)
  {
    OUsuper* ou = (OUsuper*) ex;

    const double eps = ou->getEpsilon();

#ifndef MONOCAR_NO_DEBUG
      assert(eps != 0);
#endif

    if (debugLevel > 1) {
      COUT << "Computing Hessian" << endl;
      for (uword i = 0; i < n; i++)
        for (uword j = i; j < n; j++)
          COUT << "_";
      COUT << endl;
    }

    MyMatrix param(n, 1);
    for (uword i = 0; i < n; i++)
      param(i, 0) = p[i];


      for (uword i = 0; i < n; i++)
    {
      hess[i+i*n] = (4 * numsecderivdiag(ou, param, eps/2, i) - numsecderivdiag(ou, param, eps, i))/3;
      for (uword j = i+1; j < n; j++)
        {
          CHECK_USER_INTERRUPT;
          hess[i+j*n] = (4 * numsecderiv(ou, param, eps/2, i, j) - numsecderiv(ou, param, eps, i, j))/3;
          hess[j+i*n] = hess[i+j*n];
        }
      if (debugLevel > 1)
        printDot(".");
      if (debugLevel == 1 || debugLevel == 0)
        printDot("*");
    }
      
      if (debugLevel > 1)
    COUT << endl;
  }

  static void ouhess(uword n, const double *p,
             double *hess,
             void *ex)
    {
      OUsuper* ou = (OUsuper*) ex;
      const double eps = ou->getEpsilon();
      if (eps <= 0.0)
    ouhess_analytic(n, p, hess, ex);
      else
    ouhess_numeric(n, p, hess, ex);
    }

}

static double oufg(unsigned n, const double *xin, double *grad, void *ex)
{
  double loglik = oufun(n, xin, ex);
  if (grad) {
    ougrad(n, xin, grad, ex);
    if (debugLevel == 1)
      printDot(":");
    if (debugLevel == 0)
      printDot(".");
  }
  else {
    if (debugLevel == 1)
      printDot(".");
  }
  return loglik;
}
#endif
