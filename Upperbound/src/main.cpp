#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>

#include <boost/thread.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <itpp/itcomm.h>

#include <UpperBound.h>
#include <boost/program_options.hpp>

int main( int argc , char ** argv )
{

	int N = 2;
	int K = 2;
	int T = 4;
	double endSNR = 0;
	double startSNR = 0;
	double deltaSNR = 1;

	boost::program_options::options_description desc("Allowed options");
	desc.add_options()
	("K",boost::program_options::value<int>()->implicit_value(K)," Number of users ")
//	("N",boost::program_options::value<int>()->implicit_value(N)," Number of receive anternna ")
	("T",boost::program_options::value<int>()->implicit_value(T)," Coherence time ")
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
N = K; // symetrical 

if( varMap.count("T") )
{
	T = varMap["T"].as<int>();
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

std::cout<<" K = "<<K<<", N = "<<N<<", T = "<<T<<", SNR = "<<startSNR<<":"<<deltaSNR<<":"<<endSNR<<"\n";

boost::thread_group threadGroupObj;

itpp::vec resultHolder(SNRCount);
for( int iter = 0;iter<SNRCount;iter++)
{

	UpperBound  mutInfoObj(N,K,T,itpp::inv_dB(startSNR+iter*deltaSNR),resultHolder(iter));
	boost::thread *curThread = new boost::thread(mutInfoObj);
	threadGroupObj.add_thread(curThread);

}

threadGroupObj.join_all();

boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
std::stringstream ss;
ss<<"InfoValsUpperBound_K"<<K<<"_N"<<N<<"_T"<<T<<"_"<<boost::posix_time::to_iso_string(now)<<".it";
itpp::it_file outFile;
outFile.open(ss.str().c_str());
outFile<<itpp::Name("InfoValsUpperBound")<<resultHolder;
std::cout<<resultHolder<<" \n Done running everything \n";

return 0;

}



