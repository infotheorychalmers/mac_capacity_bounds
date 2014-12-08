#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>

#include <boost/thread.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <itpp/itcomm.h>

#include <MutualInformation.h>
#include <boost/program_options.hpp>

int main( int argc , char ** argv )
{

	int N = 4;
	int K = 4;
	int T = 8;
	size_t MCX = 1000;
	size_t MCY = 1000;
	double startSNR = 0;
	double endSNR = 40;
	double deltaSNR = 1;

	boost::program_options::options_description desc("Allowed options");
	desc.add_options()
	("K",boost::program_options::value<int>()->implicit_value(K)," Number of users ")
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

	if( varMap.count("K") )
	{
		K = varMap["K"].as<int>();
	}

//if( varMap.count("N") )
//{
//	N = varMap["N"].as<int>();
//}

N = K; // symmetric

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

int SNRCount = (endSNR-startSNR)/deltaSNR + 1;

std::cout<<" K = "<<K<<", N = "<<N<<", T = "<<T<<", MCX = "<<MCX<<", SNR = "<<startSNR<<":"<<deltaSNR<<":"<<endSNR<<"\n";

itpp::Array<itpp::vec> eigV(MCX);
itpp::vec logDenom(MCX);

for ( int mcIter = 0 ; mcIter < MCX ; mcIter++ )
{
	itpp::cmat V(K,T);
	itpp::randn_c(K,T,V);
	itpp::vec & eigVRef = eigV(mcIter);
	itpp::eig_sym(V*V.hermitian_transpose(),eigVRef);

	double &tempVal = logDenom[mcIter];

	for(int i=K-1;i>=1;i--)
	{
		double &eigVRefi = eigVRef[i];
		for(int j=i-1;j>=0;j--)
		{
			tempVal += std::log(eigVRefi-eigVRef[j]);
		}
	}
	tempVal += (T-K)*itpp::sum(itpp::log(eigVRef));

}

boost::thread_group threadGroupObj;


itpp::vec resultHolder(SNRCount);
for( int iter = 0;iter<SNRCount;iter++)
{

	MutualInformation  mutInfoObj(N,K,T,itpp::inv_dB(startSNR+iter*deltaSNR),MCX,MCY,resultHolder(iter),eigV,logDenom);
	boost::thread *curThread = new boost::thread(mutInfoObj);
	threadGroupObj.add_thread(curThread);

}

threadGroupObj.join_all();


// writing to file
boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
std::stringstream ss;
ss<<"InfoValsGaussian_K"<<K<<"_N"<<N<<"_T"<<T<<"_"<<boost::posix_time::to_iso_string(now)<<".it";
itpp::it_file outFile;
outFile.open(ss.str().c_str());
outFile<<itpp::Name("InfoValsGaussian")<<resultHolder;
std::cout<<resultHolder<<" \n Done running everything \n";
return 0;

}
