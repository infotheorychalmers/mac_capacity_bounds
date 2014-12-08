#include <stdlib.h>

#include <itpp/itbase.h>

#define EXPINF 500.f

#define LOG_EXP_MINUS(x,y) (((x)>((y)+EXPINF))?(x):((y)+log1p(expm1((x)-(y))-1)))
#define LOG_EXP_PLUS(x,y) ((std::abs((x)-(y))>EXPINF)?(std::max(x,y)):(std::min(x,y)+log1p(expm1(std::abs((x)-(y)))+1)))



class MutualInformation{
	
	const int N;
	const int U;
	const std::vector<int> &A;
	int K;
	const int T;
	const double SNR;
	const size_t MCX;
	const size_t MCY;
	const itpp::Array<itpp::vec> &m_EigV_Vals;
	const itpp::Array<itpp::vec> &m_EigM_Vals;
	const itpp::vec &m_LogDenomVals;
	double &m_ResultHolder;
	
	double ln_gamma_inc(const double &iX , const unsigned int &iN);
	
	
public:
	MutualInformation(const int &n,const int &u,const std::vector<int> &a,const int &t,const double &snr,const size_t &mcx,const size_t &mcy,
		double &resultHolder,
		const itpp::Array<itpp::vec> &eigv,
		const itpp::Array<itpp::vec> &eigm,
		const itpp::vec &logdenom
                      ); // constructor
	void operator()(); // Thread running
	
};



