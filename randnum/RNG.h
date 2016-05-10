/*
 * =====================================================================================
 *
 *       Filename:  RNG.h
 *
 *    Description:  Random Number Generator class
 *
 *        Created:  05/09/2016 12:00:05 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dilawar Singh (), dilawars@ncbs.res.in
 *   Organization:  NCBS Bangalore
 *
 * =====================================================================================
 */


#ifndef  RNG_INC
#define  RNG_INC

#ifdef  USE_BOOST
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/random_device.hpp>
#else      /* -----  not USE_BOOST  ----- */

#ifdef  ENABLE_CPP11
#include <random>
#elif USE_GSL      /* -----  not ENABLE_CPP11 and using GSL  ----- */
#include <ctime>
#include <gsl/gsl_rng.h>
#endif     /* -----  not ENABLE_CPP11  ----- */

#endif     /* -----  not USE_BOOST  ----- */

#include <limits>

namespace moose {

/* 
 * =====================================================================================
 *        Class:  RNG
 *  Description:  Random number generator class.
 * =====================================================================================
 */

template < typename T >
class RNG
{
    public:
        // ====================  LIFECYCLE     =======================================
        RNG ()                           /* constructor      */
        {
            // Setup a random seed if possible.
#ifdef  ENABLE_CPP11 
            std::random_device rd;
            setSeed( rd() );
#elif USE_BOOST
            boost::random::random_device rd;
            setSeed( rd() );
#else      /* -----  not ENABLE_CPP11  ----- */

            gsl_r_ = gsl_rng_alloc( gsl_rng_default );
            gsl_rng_set( gsl_r_, time(NULL) );

#endif     /* -----  not ENABLE_CPP11  ----- */

        }

        RNG ( const RNG &other ); /* copy constructor */

        ~RNG ()                                     /* destructor       */
        {

#if defined(USE_BOOST) || defined(ENABLE_CPP11) 
#else
            gsl_rng_free( gsl_r_ );
#endif

        }

        /* ====================  ACCESSORS     ======================================= */
        T getSeed( void )
        {
            return seed_;
        }

        /* ====================  MUTATORS      ======================================= */
        void setSeed( const T seed )
        {
#if defined(USE_BOOST) || defined(ENABLE_CPP11)
            seed_ = seed;
            rng_.seed( seed_ );
#else 
            gsl_rng_set(gsl_r_, seed );
#endif
        }

        /**
         * @brief Generate a uniformly distributed random number between a and b.
         *
         * @param a Lower limit (inclusive)
         * @param b Upper limit (exclusive).
         */
        T uniform( const T a, const T b)
        {
            size_t maxInt = std::numeric_limits<int>::max();

#if defined(USE_BOOST) || defined(ENABLE_CPP11)
            return ( (b - a ) * dist_( rng_ ) / maxInt ) + a;
#else
            return ( (b -a ) * gsl_rng_get( gsl_r_ ) / gsl_rng_max( gsl_r_ ) + a );
#endif
        }

        /**
         * @brief Return a uniformly distributed random number between 0 and 1
         * (exclusive).
         *
         * @return randum number.
         */
        T uniform( void )
        {
#if defined(USE_BOOST) || defined(ENABLE_CPP11)
            return dist_( rng_ ) / std::numeric_limits<int>::max();
#else
            return gsl_rng_uniform( gsl_r_ );
#endif
        }


    private:
        /* ====================  DATA MEMBERS  ======================================= */
        T res_;
        T seed_;

#if USE_BOOST
        boost::random::mt19937 rng_;
        boost::random::uniform_int_distribution<> dist_;
#elif ENABLE_CPP11
        std::mt19937 rng_;
        std::uniform_int_distribution<> dist_;
#else      /* -----  not ENABLE_CPP11  ----- */
        gsl_rng* gsl_r_;
#endif     /* -----  not ENABLE_CPP11  ----- */

}; /* -----  end of template class RNG  ----- */


}                                               /* namespace moose ends  */

#endif   /* ----- #ifndef RNG_INC  ----- */
