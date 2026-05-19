/* Matrix exponential using Pade approximation */

MyMatrixCol MatrixExpC_pm(const uint iM)
{
  MyMatrixCol C_pm(iM+1);
  C_pm[iM]=1.0;
  for (uint j = iM; j > 0; j--)
    C_pm[j-1] = C_pm[j] * (double)(j) / (double)((iM - j + 1) * (iM + j));
  return C_pm;
}

class MatrixExp
{
public:
  MatrixExp(const MyMatrix& iH)
    :mH(iH)
    ,mHnorm(sum(abs(iH),1).max())
  {
#ifndef MOU_NO_DEBUG
    assert(mH.is_square());
    assert(mM >= 3);
#endif
  }

  MatrixExp(const uint iN)
    :mH(iN,iN,fill::zeros)
    ,mHnorm(0)
  {
  }

  MatrixExp(const MatrixExp& iMatrixExp)
    :mH(iMatrixExp.getH())
    ,mHnorm(iMatrixExp.getHnorm())
  {
  }

  void init(void)
  {
    mHnorm = sum(abs(mH),1).max();
  }

  MyMatrix operator()(double scale) const
  {
    double Hnorm = abs(scale * mHnorm);
    if (Hnorm == 0.0)
      {
#ifndef MOU_NO_DEBRG
        assert(mH.is_finite());
#endif
	return eye(mH.n_rows,mH.n_rows);
      }
    uint s = 0;
    if (Hnorm >= 0.25)
      {
    	s = (uint)(log2(Hnorm)+2.0);
    	scale /= exp2(s);
      }
    MyMatrix R = scale * mH;
    MyMatrix Rsq = R * R;
    MyMatrix Q = mC_pm[0] * Rsq;
    MyMatrix P = mC_pm[1] * Rsq;
    Q.diag() += mC_pm[2];
    P.diag() += mC_pm[3];
    uint j = 3;
    while (++j < mC_pm.n_elem)
      {
	Q *= Rsq;
	Q.diag() += mC_pm[j++];
	if (j < mC_pm.n_elem)
	  {
	    P *= Rsq;
	    P.diag() += mC_pm[j];
	  }
      }
    if (j < mC_pm.n_elem)
      {
	Q *= R;
	Q -= P;
	R = -2.0 * solve(Q, P);
	R.diag() -= 1.0;
      }
    else
      {
	P *= R;
        Q -= P;
	R = 2.0 * solve(Q, P);
	R.diag() += 1.0;
      }
    for (j = 0; j < s; j++)
      R *= R;
    return R;
  }

  MyMatrix& getH(void)
  {
    return mH;
  }
  const MyMatrix& getH(void) const
  {
    return mH;
  }
  const double getHnorm(void) const
  {
    return mHnorm;
  }
 private:
  MyMatrix mH;
  double mHnorm;
  static const MyMatrixCol mC_pm;
};

// Argument must be >=3
const MyMatrixCol MatrixExp::mC_pm = MatrixExpC_pm(6);
