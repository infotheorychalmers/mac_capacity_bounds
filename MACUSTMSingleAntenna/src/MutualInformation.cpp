
#include <iostream>

#include <itpp/itbase.h>
#include <gsl/gsl_cdf.h>

#include <MutualInformation.h>

MutualInformation::MutualInformation(const int &n,const int &k,const int &t,const double &snr,const size_t &mcx,const size_t &mcy,
	double &resultHolder,
	const itpp::Array<itpp::vec> &eigV,
	const itpp::vec &logdenom
	):
N(n),
K(k),
T(t),
SNR(snr),
MCX(mcx),
MCY(mcy),
m_ResultHolder(resultHolder),
m_EigV_Vals(eigV),
m_LogDenomVals(logdenom)
{

}

void MutualInformation::operator()(){

	double tempVal = 0.f;
	double beta = T*SNR/K;
	double beta_dash = 1/beta;
	const double & sqrt_beta = sqrt(beta);

	double hyVal = 0.f;
	double hy_xVal = 0.f;
	const int &P = std::max(K,N);
	const int &L = std::min(K,N);


	for( size_t mcyIter=0 ; mcyIter<MCY ; mcyIter++ )
	{

		itpp::cmat X(K,T);
		itpp::randn_c(K,T,X);

		for( int rIter = 0; rIter < K ; ++rIter )
		{
			double normVal = 0.f;
			for( int cIter = 0 ; cIter < T ; ++cIter )
			{
				normVal += std::norm(X(rIter,cIter));
			}
			normVal = sqrt(normVal);
			for( int cIter = 0 ; cIter < T ; ++cIter )
			{
				X(rIter,cIter) /= normVal;
			}
		}

		itpp::cmat S(N,K);
		itpp::randn_c(N,K,S);

		itpp::cmat W(N,T);
		itpp::randn_c(N,T,W);

		const itpp::cmat & Y = sqrt_beta*S*X+W;

		itpp::vec eigY;
		itpp::eig_sym(Y*Y.hermitian_transpose(),eigY);

		tempVal = itpp::sum(eigY)+(T-N)*itpp::sum(itpp::log(eigY));
		for(int i=N-1;i>=1;i--)
		{
			double & eigYi = eigY[i];
			for(int j=i-1;j>=0;j--)
			{
				tempVal += std::log(eigYi-eigY[j]);
			}
		}

		double logExpect = 0.0;

		for( size_t mcxIter=0 ; mcxIter<MCX ; mcxIter++ )
		{
			itpp::mat M(P,P);
			const itpp::vec & eigV = m_EigV_Vals(mcxIter); 
			for( int rIter = P-1 ; rIter >=0 ; --rIter )
			{
				double &y_i = eigY[rIter];
				for( int cIter = P-1 ; cIter >=0 ; --cIter )
				{
					const double & l_j = eigV[cIter];
					const double & b_j = l_j/(beta_dash+l_j);
					const double & m_ij = y_i*b_j;
					M(P-1-rIter,P-1-cIter) = m_ij+((m_ij<=0.0)?0.0:(log1p(-gsl_cdf_poisson_P(T-P-1,m_ij))));
				}
			}
			itpp::mat S(P,P);
			S.ones();
			double logDet = M(0,0);
			for( int step = 0; step < (P-1) ; step ++ )
			{
				for( int rIter = (step+1) ; rIter < P ; rIter++ )
				{
					for( int cIter = (step+1) ; cIter < P ; cIter++ )
					{
						const double &E2 = M(rIter,step)+M(step,cIter)-M(step,step);
						const double &S2 = S(rIter,step)*S(step,cIter)*S(step,step);
						double E1 = M(rIter,cIter);
						double S1 = S(rIter,cIter);
						const double &diffSign = (((E1-E2)>0)?1:(-1));
						S(rIter,cIter) = (((S1-S2+(S1+S2)*diffSign)>0)?1:(-1));
						const double &diff = std::abs(E1-E2);
						M(rIter,cIter) = (diff<EXPINF)?(std::min(E1,E2)+log1p(expm1(diff)-S1*S2)):(std::max(E1,E2));
					}
				}
				logDet += M(step+1,step+1);
			}

			const double &accVal = logDet+(T-N-1)*itpp::sum(itpp::log(beta_dash+eigV))-m_LogDenomVals(mcxIter);
			logExpect = ((mcxIter)?(itpp::log_add(accVal,logExpect)):(accVal));

		}
		hyVal += (( tempVal-logExpect + std::log(MCX) ) / MCY);
	}
	for( int iter = (T-L+1) ; iter <= T ; ++iter )
	{
		hyVal -= std::log(itpp::gamma(iter));
	}
	hyVal -= K*(K-1)*std::log(beta)/2;

	for( size_t mcxIter=0 ; mcxIter<MCX ; mcxIter++ )
	{
		const itpp::vec & eigV = m_EigV_Vals(mcxIter);
		hy_xVal += N*(itpp::sum(itpp::log(beta_dash+eigV))/ MCX);
	}
	hy_xVal += N*T;

	m_ResultHolder = (hyVal - hy_xVal) / ( T* std::log(2));
}


