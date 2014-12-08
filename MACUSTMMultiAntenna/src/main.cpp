#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>

#include <boost/thread.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/program_options.hpp>
#include <itpp/itcomm.h>
#include <gsl/gsl_sf_gamma.h>

#include <MutualInformation.h>

int main( int argc , char ** argv )
{
	
	int N = 4;
	int U = 4;
	int T = 8;
	size_t MCX = 1000;
	size_t MCY = 1000;
	double startSNR = 0;
	double endSNR = 40;
	double deltaSNR = 1;
	std::vector<int> A;
	A.push_back(1);
	
	boost::program_options::options_description desc("Allowed options");
	desc.add_options()
	("U",boost::program_options::value<int>()->implicit_value(U)," Number of users ")
	("A",boost::program_options::value<std::vector<int> >()->multitoken()->implicit_value(A,"1")," Antenna distribution space separated values")
//	("N",boost::program_options::value<int>()->implicit_value(N)," Number of receive anternna ")
	("T",boost::program_options::value<int>()->implicit_value(T)," Coherence time ")
	("MCX",boost::program_options::value<int>()->implicit_value(MCX)," Monte carlo ")
	("startSNR",boost::program_options::value<double>()->implicit_value(startSNR)," SNR start value ")
	("endSNR",boost::program_options::value<double>()->implicit_value(endSNR)," SNR end vaule ")
	("deltaSNR",boost::program_options::value<double>()->implicit_value(deltaSNR)," SNR delta value ")
	("help"," Print help ")
	;
	
	boost::program_options::variables_map varMap;
	boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), varMap);
	boost::program_options::notify(varMap);
	
	if (varMap.count("help")) {
		std::cout << desc << "\n";
		return 1;
	}
	
	
	if( varMap.count("U") )
	{
		U = varMap["U"].as<int>();
	}
	
//	if( varMap.count("N") )
//	{
//		N = varMap["N"].as<int>();
//	}

	
	if( varMap.count("T") )
	{
		T = varMap["T"].as<int>();
	}
	
	if( varMap.count("MCX") )
	{
		MCX = varMap["MCX"].as<int>();
		MCY = MCX;
	}
	
	if( varMap.count("startSNR") )
	{
		startSNR = varMap["startSNR"].as<double>();
	}
	
	if( varMap.count("endSNR") )
	{
		endSNR = varMap["endSNR"].as<double>();
	}
	
	if( varMap.count("deltaSNR") )
	{
		deltaSNR = varMap["deltaSNR"].as<double>();
	}
	
	if( varMap.count("A") )
	{
		A = varMap["A"].as<std::vector<int> >();
	}
	
	int userCount = A.size();
	
	if( userCount!=U)
	{
		A.resize(U);
		for(int iter=userCount;iter<U;++iter)
		{
			A[iter] = A[userCount-1];
		}
	}
	
	int K = 0;
	for (int iter=0; iter<U; ++iter) {
		K += A[iter];
	}
	
	N = K; // symmetrical
	
	int SNRCount = (endSNR-startSNR)/deltaSNR + 1;
	
	std::cout<<" U = "<<U;
	std::cout<<", A = [ ";
		for (int iter=0; iter<U; ++iter) {
			std::cout<<A[iter]<<" ";
		}
		std::cout<<"], K = "<<K<<", N = "<<N<<", T = "<<T<<", MCX = "<<MCX<<", SNR = "<<startSNR<<":"<<deltaSNR<<":"<<endSNR<<"\n";
std::cout<<"\n";

itpp::Array<itpp::vec> eigV(MCX);
itpp::Array<itpp::vec> eigM(MCX);
itpp::vec logDenom(MCX);

double epsilon = 5e-15;

for ( int mcIter = 0 ; mcIter < MCX ; mcIter++ )
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
	
	itpp::vec & eigVRef = eigV(mcIter);
	itpp::vec & eigMRef = eigM(mcIter);
	itpp::vec eigVals(K);
	itpp::eig_sym(V*V.hermitian_transpose(),eigVals);
	
	double & logDenomRef  = logDenom[mcIter];
	int ecount = 0;
	eigVRef.ins(ecount,eigVals[0]);
	eigMRef.ins(ecount,1);
	logDenomRef = 0;
	for (int iter=1; iter<K; ++iter) {
		if( (eigVals[iter]-eigVRef[ecount])< epsilon )
		{
			++eigMRef[ecount];
			logDenomRef += gsl_sf_lngamma(eigMRef[ecount]);
		}else{
			++ecount;
			eigVRef.ins(ecount,eigVals[iter]);
			eigMRef.ins(ecount,1);
		}
		
	}
	logDenomRef += N*itpp::sum(itpp::elem_mult(itpp::log(eigVRef),eigMRef));
	
}

boost::thread_group threadGroupObj;


itpp::vec resultHolder(SNRCount);
for( int iter = 0;iter<SNRCount;iter++)
{
	
	MutualInformation mutInfoObj(N,U,A,T,itpp::inv_dB(startSNR+iter*deltaSNR),MCX,MCY,resultHolder(iter),eigV,eigM,logDenom);
	boost::thread *curThread = new boost::thread(mutInfoObj);
	threadGroupObj.add_thread(curThread);
	
}

threadGroupObj.join_all();


	// writing to file
boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
std::stringstream ss;
ss<<"InfoValsUniform_U"<<U<<"_A_";
for (int iter=0; iter<U; ++iter) {
	ss<<A[iter]<<"_";
}
ss<<"_N"<<N<<"_T"<<T<<"_"<<boost::posix_time::to_iso_string(now)<<".it";
itpp::it_file outFile;
outFile.open(ss.str().c_str());
outFile<<itpp::Name("InfoValsUniform")<<resultHolder;
std::cout<<resultHolder<<" \n Done running everything \n";
return 0;

}
