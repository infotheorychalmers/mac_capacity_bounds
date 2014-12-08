
#include <iostream>

#include <itpp/itbase.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>

#include <MutualInformation.h>

MutualInformation::MutualInformation(const int &n,const int &k,const int &t,const double &snr,const size_t &mcx,const size_t &mcy,
	double &resultHolder
	):
N(n),
K(k),
T(t),
SNR(snr),
MCX(mcx),
MCY(mcy),
m_ResultHolder(resultHolder)
{

}

void MutualInformation::operator()(){

	double tempVal = 0.f;
	double beta = T*SNR/K;
	const double & sqrt_beta = sqrt(beta);
	const double & beta_dash = beta/(1+beta);

	double hyVal = 0.f;
	double hy_xVal = 0.f;
	const int &P = std::max(K,N);
	const int &L = std::min(K,N);


	for( size_t mcyIter=0 ; mcyIter<MCY ; mcyIter++ )
	{

		itpp::cmat X(K,T);
		itpp::randn_c(K,T,X);
		itpp::cmat U(K,K);
		itpp::cmat V(T,T);
		itpp::vec D(K);
		itpp::svd(X,U,D,V);
		for(int row=0;row<K;row++){
			for( int col=0;col<T;col++){
				X(row,col)=V(row,col);
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
		itpp::mat M(P,P);
		for( int rIter = P-1 ; rIter >=0 ; --rIter )
		{
			double &y_i = eigY[rIter];
			for( int cIter = P-1 ; cIter >=0 ; --cIter )
			{
				const double &localVal = y_i*beta_dash;
				M(P-1-rIter,P-1-cIter) = cIter*log(y_i)+localVal+((localVal<=0.0)?0.0:(log1p(-gsl_cdf_poisson_P(T-P-1-cIter,localVal))));
			}
		}
		itpp::mat SignMat(P,P);
		SignMat.ones();
		double logDet = M(0,0);
		for( int step = 0; step < (P-1) ; step ++ )
		{
			for( int rIter = (step+1) ; rIter < P ; rIter++ )
			{
				for( int cIter = (step+1) ; cIter < P ; cIter++ )
				{
					const double &E2 = M(rIter,step)+M(step,cIter)-M(step,step);
					const double &S2 = SignMat(rIter,step)*SignMat(step,cIter)*SignMat(step,step);
					double E1 = M(rIter,cIter);
					double S1 = SignMat(rIter,cIter);
					const double &diffSign = (((E1-E2)>0)?1:(-1));
					SignMat(rIter,cIter) = (((S1-S2+(S1+S2)*diffSign)>0)?1:(-1));
					const double &diff = std::abs(E1-E2);
					M(rIter,cIter) = (diff<EXPINF)?(std::min(E1,E2)+log1p(expm1(diff)-S1*S2)):(std::max(E1,E2));
				}
			}
			logDet += M(step+1,step+1);
		}
		tempVal -= logDet;


		hyVal += (tempVal/MCY);
	}
	hyVal = N*K*log(beta)+K*(T-K-N)*log(beta_dash)+hyVal;
	for( int iter=0;iter<P;iter++)
	{
		hyVal += (gsl_sf_lngamma(iter+1)-gsl_sf_lngamma(iter+1+T-P));
	}
	hy_xVal = N*T+N*K*log(1+beta);

	m_ResultHolder = (hyVal - hy_xVal) / ( T* std::log(2));
}


