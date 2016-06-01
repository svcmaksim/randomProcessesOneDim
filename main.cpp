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

int main( int argc, char** argv )
{
    if( argc != 4 )
        throw std::logic_error( "bad params" );

    const size_t NUM_OF_TRAJECTORIES = atoi( argv[1] );//20;
    const double_t n_param = strtod( argv[2], nullptr );
    const double_t x0 = strtod( argv[3], nullptr );

    std::cout << "x0 = " << x0 << "\t n = " << n_param << std::endl;

    auto a_func = [&]( double_t t, double_t X ) -> double_t
    {
        return -3./2. * X / n_param;
    };

    auto b_func = [&]( double_t t, double_t X ) -> double_t
    {
        return X  / sqrt( n_param /* * (n_param - 1) */);
    };

    //auto a_func = [&]( double_t t, double_t X ) -> double_t
    //{
    //    return (-t /*- 2.0*/) * X / (1 + pow((/*2.0 -*/ t), 2)) - (X >= 0) ? 1. : -1.;
    //};

    //auto b_func = [&]( double_t t, double_t X ) -> double_t
    //{
    //    return sqrt(fabs(X));
    //};*/

    std::vector<trajectory_t> trajectories;
    for( size_t i = 0; i < NUM_OF_TRAJECTORIES; ++i )
    {
        std::fstream fStream( "res_5" + std::to_string( i ) + ".csv", std::ios::out );
        trajectory_t tr = getTrajectory(
                    0,
                    2,
                    a_func,
                    b_func,
                    pow( x0, n_param ),
                    1000
                    );
        trajectories.push_back( tr );
        writeTrajextoryToStream( tr, fStream );
    }

    trajectory_t mean( 1000 );

    for( size_t i = 0; i < 1000; ++i )
    {
        mean[i].first = 0.002 * i;
        mean[i].second = 0;
        for( size_t j = 0; j < trajectories.size(); ++j )
            mean[i].second += trajectories[j][i].second / trajectories.size();
    }

    std::fstream meanStr( "mean5.csv", std::ios::out );

    writeTrajextoryToStream( mean, meanStr );

    return EXIT_SUCCESS;
}
