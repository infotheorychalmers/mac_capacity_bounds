
#include <iostream>

#include <itpp/itbase.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>

#include <UpperBound.h>

UpperBound::UpperBound(const int &n,const int &k,const int &t,const double &snr,
	double &resultHolder
	):
N(n),
K(k),
T(t),
P(std::max(n,k)),
SNR(snr),
m_ResultHolder(resultHolder),
beta(snr*t/k),
m_DervKonstMats(),
l_M11(P,P),
l_M12(P,T-P),
l_M11k(P,P),
l_M11d(P,P)
{
	
	
	m_Konst = -P*T+P*P*log(beta)+P*T/beta;
	for( int iter=0;iter<P;++iter)
	{
		m_Konst += (gsl_sf_lngamma(iter+1)-gsl_sf_lngamma(iter+1+T-P));
	}
	
	itpp::mat M22(T-P,T-P);
	for(int row=0;row<(T-P-1);++row)
	{
		for(int col=0;col<(T-P);++col)
		{
			M22(row,col) = ((col>row)?0:((1-2*((row-col)%2))/gsl_sf_beta(row-col+1,T-P-row-1)));
		}
	}
	for(int col=0;col<T-P;++col){
		M22(T-P-1,col) = 1-2*((T-P-col-1)%2);
	}
	
	itpp::mat M21(T-P,P);
	for(int row=0;row<(T-P-1);++row){
		for (int col=0; col<P ; ++col) {
			M21(row,col) = 1/gsl_sf_beta(P-col,T-P-row-1);
		}
	}
	for (int col=0; col<P ; ++col) {
		M21(T-P-1,col) = 1;
	}
	
	itpp::mat M21d(T-P,P);
	for(int row=0;row<(T-P-1);++row){
		for (int col=0; col<P ; ++col) {
			M21d(row,col) = gsl_sf_psi_int(T-row-col-1)/gsl_sf_beta(P-col,T-P-row-1);
		}
	}
	for (int col=0; col<P ; ++col) {
		M21d(T-P-1,col) = gsl_sf_psi_int(P-col);
	}
	
	itpp::mat baseKonstMat(T-P,P);
	itpp::backslash(M22,M21,baseKonstMat);
	
	itpp::mat derivKonstMatrix(T-P,P);
	itpp::backslash(M22,M21d,derivKonstMatrix);
	
	for (int k=0; k<P ; ++k) {
		itpp::mat curMat(baseKonstMat);
		for (int row=0; row<(T-P) ; ++row) {
			curMat(row,k) = derivKonstMatrix(row,k);
		}
		m_DervKonstMats.push_back(curMat);
	}
	
	
}

int UpperBound::m_LPrecision = 8;
float UpperBound::m_LMinStart = 0.f;
float UpperBound::m_LMaxStart = 1000.f;
float UpperBound::m_LDeltaStart = 100.f;
float UpperBound::m_LDeltaDivide = 10.f;

int UpperBound::m_XPrecision = 8;
float UpperBound::m_XMinStart = 1.0f;
float UpperBound::m_XMaxStart = 15.f;
float UpperBound::m_XDeltaStart = 0.5f;
float UpperBound::m_XDeltaDivide = 10.f;

void UpperBound::operator()(){
	itpp::vec xMatch(P);
	double delta = m_LDeltaStart;
	double lMin = m_LMinStart;
	double lMax = m_LMaxStart;
	double lMatch = 0.f;
	double infVal = IterateL(lMin,lMax,delta,lMatch,xMatch);
	double prevInfVal = infVal;
	for(int iter=0;iter<m_LPrecision;++iter){
		const double nDelta = delta/m_LDeltaDivide;
		lMin = std::max<double>(lMatch-delta,m_LMinStart);
		lMax = lMatch + delta;
		infVal = IterateL(lMin,lMax,nDelta,lMatch,xMatch);
		prevInfVal = infVal;
		delta = nDelta;
	}
	
	std::cout<<" infVal = "<<infVal<<", l = "<<lMatch<<", x = "<<xMatch<<"\n";
	
	m_ResultHolder = infVal;
	
	
}

inline double UpperBound::SupX(const double &lambda , itpp::vec &xMatch)
{
	const double &maxVal = 31.f*T/16.f;
	const double &minVal = T/16.f;
	double delta = T/16.f;
	itpp::vec xMin(P);
	for(int iter=0;iter<P;++iter){
		xMin(iter) = minVal;
	}
	itpp::vec xMax(P);
	for(int iter=0;iter<P;++iter){
		xMax(iter) = maxVal;
	}
	
	double supVal = IterateX(lambda,xMin,xMax,delta,xMatch);
	double prevSupVal = supVal;
	for(int iter=0;iter<m_XPrecision;++iter){
		const double & nDelta = delta/m_XDeltaDivide;
		xMin = xMatch-delta; VecMin(xMin,minVal);
		supVal = IterateX(lambda,xMin,xMatch+delta,nDelta,xMatch);
		prevSupVal = supVal;
		delta = nDelta;
	}
	return supVal;
}

inline double UpperBound::IterateL(const double &lStart,const double & lEnd,const double &lIter,double &lMatch,itpp::vec &xMatch)
{
	const double & firstVal = SupX(lStart,xMatch);
	lMatch = lStart;
	double infVal = firstVal;
	for(double lambda =  lStart+lIter ; lambda < lEnd ; lambda += lIter)
	{
		const double & retVal = SupX(lambda,xMatch);
		if( retVal<infVal ){
			infVal = retVal;
			lMatch = lambda;
		}
	}
    //	std::cout<<" lStart = "<<lStart<<" lEnd = "<<lEnd<<" lMatch = "<<lMatch;
	return infVal;
	
}


inline double UpperBound::IterateX(const double &lambda,const itpp::vec  &xMin,const itpp::vec &xMax,const double &xIter , itpp::vec &xMatch)
{
	itpp::vec u_Max = xMax;
	for(int iter=1;iter<P;++iter){
		u_Max(iter) = std::min<double>(u_Max(iter-1)-xIter,u_Max(iter));
	}
	itpp::vec u_Min = xMin;
	for(int iter=P-2;iter>=0;--iter){
		u_Min(iter) = std::max<double>(u_Min(iter+1)+xIter,u_Min(iter));
	}
	itpp::vec dmat = u_Min;
	double supVal = FunctionEvaluate(SNR*dmat/P,lambda);
	xMatch = dmat;
	bool proceed;
	do{
		proceed = false;
		for(int iter=0;iter<P;++iter)
		{
			if(!((dmat(iter)+xIter)>u_Max(iter))){
				dmat(iter) += xIter;
				for(int uIter = iter-1 ; uIter>=0;--uIter){
					dmat(uIter) = std::max<double>(dmat(uIter+1)+xIter,u_Min(uIter));
				}
				proceed = true;
                //std::cout<<" sum = "<<(T*P-itpp::sum(dmat))<<"\n";
				const double & rVal = FunctionEvaluate(SNR*dmat/P,lambda);
				if(rVal>supVal)
				{
					supVal = rVal;
					xMatch = dmat;
				}
				break;
			}
		}
	}while(proceed);
	return (m_Konst+supVal)/(T*log(2));
	
}

inline void UpperBound::VecMin(itpp::vec &iVec,const double &iVal)
{
	for(int iter=iVec.size()-1;iter>=0;--iter){
		iVec(iter) = std::max<double>(iVec(iter),iVal);
	}
}

inline double UpperBound::FunctionEvaluate(const itpp::vec &d_sq, const double &lambda)
{
	const itpp::vec &I_d_sq = 1+d_sq;
	const double &trace_d_sq = itpp::sum(d_sq);
	const double &det_d_sq = itpp::prod(d_sq);
	const double &det_I_d_sq = itpp::prod(I_d_sq);
	
	double v_d_sq = 1.f;
	
	for(int iter1=0;iter1<(P-1);++iter1)
	{
		const double &di = d_sq(iter1);
		for(int iter2=iter1+1;iter2<P;++iter2)
		{
			v_d_sq *= (di-d_sq(iter2));
		}
	}
	
	for (int row = 0 ; row < P ; ++row) {
		const double &neg_bi_inv = I_d_sq(row);
		for(int col=0;col<(T-P);++col){
			l_M12(row,col) = (1-2*((T-P-col-1)%2))/pow(neg_bi_inv,T-P-col-1);
		}
	}
	
	for (int row=0; row<P ; ++row) {
		const double &neg_ai_inv = I_d_sq(row);
		const double &log_neg_ai_inv = log(neg_ai_inv);
		for (int col=0; col<P  ; ++col) {
			const double & baseVal = pow(neg_ai_inv,P-col);
			l_M11(row,col) = baseVal;
			l_M11d(row,col) = baseVal*(log_neg_ai_inv+gsl_sf_psi_int(P-col));
		}
	}
	
	double detSum = 0.f;
	l_M11k = l_M11;
	for (int k=0; k<P ; ++k) {
		for (int row = 0; row<P ; ++row) {
			l_M11k(row,k) = l_M11d(row,k);
		}
		detSum += itpp::det(l_M11k-l_M12*m_DervKonstMats[k]);
		for (int row = 0; row<P ; ++row) {
			l_M11k(row,k) = l_M11(row,k);
		}
	}
	
	double retVal = P*trace_d_sq/beta;
	retVal -= P*log(det_I_d_sq);
	retVal += lambda*(T*SNR-trace_d_sq);
	retVal += (T-P)*pow(det_I_d_sq,T-P-1)*detSum/(v_d_sq*pow(det_d_sq,T-P));
	return retVal;
	
}


