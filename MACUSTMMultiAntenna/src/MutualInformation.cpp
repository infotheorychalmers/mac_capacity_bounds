
#include <iostream>

#include <itpp/itbase.h>
#include <gsl/gsl_cdf.h>

#include <MutualInformation.h>

MutualInformation::MutualInformation(const int &n,const int &u,const std::vector<int> &a,const int &t,const double &snr,const size_t &mcx,const size_t &mcy,
	double &resultHolder,
	const itpp::Array<itpp::vec> &eigv,
	const itpp::Array<itpp::vec> &eigm,
	const itpp::vec &logdenom
	):
N(n),
U(u),
A(a),
T(t),
SNR(snr),
MCX(mcx),
MCY(mcy),
m_ResultHolder(resultHolder),
m_EigV_Vals(eigv),
m_EigM_Vals(eigm),
m_LogDenomVals(logdenom)
{
	K = 0;
	for (int iter=0; iter<U; ++iter) {
		K+=A[iter];
	}
}


inline double MutualInformation::ln_gamma_inc(const double &iX , const unsigned int &iN)
{
	return iX+((iX<=0.0)?0.0:(log1p(-gsl_cdf_poisson_P(iN,iX))));
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
	
	itpp::Array<itpp::vec> BxVals(MCX);
	itpp::vec logDenomVals(MCX);
	
	for( size_t mcxIter=0 ; mcxIter<MCX ; mcxIter++ ){
		const itpp::vec &eigV = m_EigV_Vals(mcxIter);
		const itpp::vec &eigM = m_EigM_Vals(mcxIter);
		itpp::vec &Bx = BxVals(mcxIter);
		const itpp::vec & ik_eigV = beta_dash + eigV;
		Bx = itpp::elem_div( eigV ,  ik_eigV);
		double & logdenom = logDenomVals[mcxIter];
		logdenom = m_LogDenomVals[mcxIter];
		const int & E = Bx.size();
		for(int i=E-1;i>=1;i--)
		{
			double & Bxi = Bx[i];
			const double & eigMi = eigM[i];
			for(int j=i-1;j>=0;j--)
			{
				logdenom += eigMi*eigM[j]*std::log(Bxi-Bx[j]);
			}
			logdenom += (T-K-N)*eigMi*std::log(Bxi);
		}
		logdenom += (T-K-N)*eigM[0]*std::log(Bx[0]);
		hy_xVal += N*(itpp::sum(itpp::elem_mult(itpp::log(ik_eigV),eigM))/ MCX);
	}
	hy_xVal += N*T;
	
	for( size_t mcyIter=0 ; mcyIter<MCY ; mcyIter++ )
	{
		
		itpp::cmat V(K,T);
		int tIter = 0;
		for (int uIter = 0 ; uIter < U ; ++uIter) {
			int aCount = A[uIter];
			itpp::cmat uX(aCount,T);
			itpp::randn_c(aCount,T,uX);
			itpp::cmat uU(aCount,aCount);
			itpp::cmat uV(T,T);
			itpp::vec uD(aCount);
			itpp::svd(uX,uU,uD,uV);
			for(int row=0;row<aCount;row++){
				for( int col=0;col<T;col++){
					V(tIter,col)=uV(row,col);
				}
				++tIter;
			}
		}
		
		itpp::cmat S(N,K);
		itpp::randn_c(N,K,S);
		
		itpp::cmat W(N,T);
		itpp::randn_c(N,T,W);
		
		const itpp::cmat & Y = sqrt_beta*S*V+W;
		
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
		
		double logExpect = 0.f;
		double logDet = 0.f;
		itpp::mat M(P,P);
		itpp::mat SignMat(P,P);
		
		for( size_t mcxIter=0 ; mcxIter<MCX ; mcxIter++ )
		{
			
			const itpp::vec &Bx = BxVals(mcxIter);
			const itpp::vec &eigM = m_EigM_Vals(mcxIter);
			const int & E = Bx.size();
			for( int rIter = P-1 ; rIter >=0 ; --rIter )
			{
				const double &y_i = eigY[rIter];
				const double &ln_y_i = std::log(y_i);
				int cIter = P-1;
				for (int eIter = E-1; eIter>=0; --eIter) {
					const int & eigMult = eigM[eIter];
					for (int col = 0 ; col < eigMult; ++col) {
						M(P-1-rIter,P-1-cIter) = ln_gamma_inc(y_i*Bx[eIter],T-P-1-col)+col*ln_y_i;
						SignMat(P-1-rIter,P-1-cIter) = 1.f;
						--cIter;
					}
				}
			}
			
			logDet = M(0,0);
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
			
			const double &accVal = logDet-logDenomVals(mcxIter);
			logExpect = ((mcxIter)?(itpp::log_add(accVal,logExpect)):(accVal));
			
		}
		hyVal += (( tempVal-logExpect + std::log(MCX) ) / MCY);
	}
	
	for( int iter = (T-L+1) ; iter <= T ; ++iter )
	{
		hyVal -= std::log(itpp::gamma(iter));
	}
	
	m_ResultHolder = (hyVal - hy_xVal) / ( T* std::log(2));
}


