#include <iostream>
#include <vector>
#include <utility>
#include <boost/random.hpp>
#include <fstream>
#include <chrono>
#include <list>
#include <algorithm>
#include <functional>

typedef double double_t;
typedef std::function<double_t(double_t, double_t)> pFunc;
typedef std::pair<double_t, double_t> coordinate_t;
typedef std::vector<coordinate_t> trajectory_t;

trajectory_t getTrajectory(
        double_t tBegin,
        double_t tEnd,
        pFunc a_func,
        pFunc b_func,
        double_t x0,
        size_t N
        )
{
    static bool firstCall = true;
    unsigned seed = 0;
    if( firstCall )
    {
        seed = std::chrono::system_clock::now().time_since_epoch().count();
        firstCall = false;
    }

    static boost::mt19937 generator( seed );
    static boost::normal_distribution<double_t> norm( 0., 1. );
    static boost::variate_generator< boost::mt19937, boost::normal_distribution<double> >
            normalSampler( generator, norm );

    trajectory_t answer( N, coordinate_t(tBegin, x0) );

    double_t tStep = (tEnd - tBegin) / (N - 1);

    for( size_t i = 1; i < N; ++i )
    {
        double curT = tStep * i + tBegin;

        double epsilon = normalSampler();

        answer[i].first = curT;
        answer[i].second =
                answer[i - 1].second +
                a_func( curT, answer[i - 1].second ) * tStep +
                epsilon * b_func( curT, answer[i - 1].second ) * sqrt( tStep );
    }

    return answer;
}

void writeTrajextoryToStream( trajectory_t &tr, std::ostream &str )
{
    for( size_t i = 0; i < tr.size(); ++i )
        str <<
               tr[i].first <<
               ',' <<
               tr[i].second << std::endl;
}

void getMeanTrajectory( const std::vector<trajectory_t> &trajectories, trajectory_t &mean )
{
    for( size_t i = 0; i < mean.size(); ++i )
    {
        mean[i].second = 0;
        for( size_t j = 0; j < trajectories.size(); ++j )
        {
            mean[i].first = trajectories[j][i].first;
            mean[i].second += trajectories[j][i].second / trajectories.size();
        }
    }
}

int main( int argc, char** argv )
{
    //if( argc != 5 )
    //    throw std::logic_error( "bad params" );

    //const size_t NUM_OF_TRAJECTORIES = atoi( argv[1] );//20;
    //const double_t n_param = strtod( argv[2], nullptr );
    //const double_t x0 = strtod( argv[3], nullptr );
    //const std::string meanName = argv[4];

    //if( NUM_OF_TRAJECTORIES % 1000 )
    //    throw std::logic_error( "number of trajectories must has 1000 as devide" );

    //std::cout << "x0 = " << x0 << "\t n = " << n_param << std::endl;
    const size_t NUM_OF_TRAJECTORIES = 1000;
	double_t n_param_begin = 1.;
	double_t n_param_end = 7.;
	size_t numExperiments = 100;
	double_t n_param = n_param_begin;	
	const double_t x0 = 1.5;	
	
	auto a_func = [&]( double_t t, double_t X ) -> double_t
    {
        return (n_param == 1) ? -X : -2. * X / n_param;
    };

    auto b_func = [&]( double_t t, double_t X ) -> double_t
    {
        return ( n_param > 1 ) ? X * sqrt( 2. / (n_param * (n_param - 1) ) ) : 1;
    };
	
	trajectory_t difference( numExperiments, std::make_pair( n_param_begin, 0 ) );
	
	for( size_t curN = 0; curN < numExperiments; ++curN )
	{
		n_param = n_param_begin + (n_param_end - n_param_begin) / numExperiments * curN;
		
		std::vector<trajectory_t> trajectories(1000);
		std::vector<trajectory_t> means( NUM_OF_TRAJECTORIES / 1000, trajectory_t( 10000 ) );
		
		for( size_t i = 1; i <= NUM_OF_TRAJECTORIES; ++i )
		{
			//std::fstream fStream( "res_5" + std::to_string( i ) + ".csv", std::ios::out );
			trajectory_t tr = getTrajectory(
						0,
						2,
						a_func,
						b_func,
						pow( x0, n_param ),
						10000
						);
			trajectories[(i-1) % 1000] = tr;
			if( i > 0 && !(i % 1000) )
				getMeanTrajectory( trajectories, means[ i / 1000 - 1] );

			//writeTrajextoryToStream( tr, fStream );
		}

		trajectory_t mean( 10000 );
		getMeanTrajectory( means, mean );

		//std::fstream meanStr( "mean" + meanName + ".csv", std::ios::out );
		
		difference[curN].first = n_param;
		
		trajectory_t::value_type &point = *std::max_element(
			mean.begin(),
			mean.end(),
            [&]( trajectory_t::value_type &pointFst, trajectory_t::value_type &pointSnd )
			{
                return fabs( pow( pointFst.second, n_param ) - pow(x0, n_param) * exp( pointFst.first ) ) -
                    fabs( pow( pointSnd.second, n_param ) - pow(x0, n_param) * exp( pointSnd.first ) );
			} );
		
		difference[curN].second = pow( point.second, n_param ) - pow( x0, n_param ) * exp( - point.first );

        std::cerr << "curN = " << curN << std::endl;
	}
	
	std::fstream diffStr( "difference.csv", std::ios::out );
	
    writeTrajextoryToStream( difference, diffStr );

    return EXIT_SUCCESS;
}
