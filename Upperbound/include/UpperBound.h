#include <stdlib.h>
#include <vector>
#include <itpp/itbase.h>

#define EXPINF 500.f

#define LOG_EXP_MINUS(x,y) (((x)>((y)+EXPINF))?(x):((y)+log1p(expm1((x)-(y))-1)))
#define LOG_EXP_PLUS(x,y) ((std::abs((x)-(y))>EXPINF)?(std::max(x,y)):(std::min(x,y)+log1p(expm1(std::abs((x)-(y)))+1)))


class UpperBound{

	const int N;
	const int K;
	const int T;
	const int P;
	const double SNR;
	double &m_ResultHolder;
	const double beta;
	double m_Konst;
	std::vector<itpp::mat> m_DervKonstMats;
	
	itpp::mat l_M11;
	itpp::mat l_M12;
	itpp::mat l_M11k;
	itpp::mat l_M11d;

	double SupX(const double &lambda , itpp::vec &retSupXVal);
	double FunctionEvaluate(const itpp::vec &dmat, const double &lambda);
	double IterateL(const double &lStart,const double & lEnd,const double &lIter,double &lMatch,itpp::vec &xMatch);
	double IterateX(const double &lambda,const itpp::vec &xMin,const itpp::vec &xMax,const double &xIter , itpp::vec &xMatch);

	void VecMin(itpp::vec &iVec,const double &iVal);
	
	
	static int m_LPrecision;
	static float m_LMinStart;
	static float m_LMaxStart;
	static float m_LDeltaStart;
	static float m_LDeltaDivide;
	
	static int m_XPrecision;
	static float m_XMinStart;
	static float m_XMaxStart;
	static float m_XDeltaStart;
	static float m_XDeltaDivide;


public:
	UpperBound(const int &n,const int &k,const int &t,const double &snr,
		double &resultHolder
				); // constructor
	void operator()(); // Thread running

};



