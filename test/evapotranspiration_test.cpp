#include "test_pch.h"
#include "core/utctime_utilities.h"
#include <armadillo>
#include <vector>
#include <chrono>
#include <boost/math/constants/constants.hpp>
#include "core/radiation.h"
#include <cmath>
#include <random>
#include <tuple>

namespace shyft::core{
    namespace evapotranspiration{
        /**\brief Evapotranspiration models
         * Penman-Monteith, ref.: ASCE-EWRI The ASCE Standardized Reference Evapotranspiration Equation Environmental and
         * Water Resources Institute (EWRI) of the American Society of Civil Engineers Task Com- mittee on Standardization of
         * Reference Evapotranspiration Calculation, ASCE, Washington, DC, Environmental and Water Resources Institute (EWRI) of
         * the American Society of Civil Engineers Task Com- mittee on Standardization of Reference Evapotranspiration Calculation,
         * ASCE, Washington, DC, 2005*/
        namespace penman_monteith {
            using namespace std;
            /**ASCE-EWRI hourly based model
             * it is stated taht this model can be used with the smaller periods as well, but we need to account for fractional hour
             * TODO: discuss how to manipulate between models and timesteps*/
             struct parameter{};
             struct response {
                 double et_ref = 0.0; // reference evapotranspiration [mm/h]
             };

             template<class P, class R>
             struct calculator {
                 P param;

                 explicit calculator(const P &p) : param(p) {}

                 // manipulating with time-steps
//                 utctimespan delta = time_axis.delta();
//                std::cout<< to_seconds(delta)  <<"  ===== "<<std::endl;
                 double ret(){
                     return 0.0;
                 //    return (inverse_lambda*Delta*(Rn-G)+gamma*Cn/(T+273)*u2*(es-ea))/(Delta+gamma*(1+Cd*u2)); // standardized reference crop evapotranspiration for short (ETos) or tall (ETrs) surfaces, [mm/h]
                 }

                 void reference_evapotranspiration_asce(R& response){

                     response.et_ref = 0.0;
                     return;
                 };

             private:
                 // ============================================ //
                 // constants for asce-ewri standard equation
                 /**Values for Cn and Cd in Eq. 1 ref.: ASCE-EWRI p.21*/
                 const double Cni[2][3] = {{900.0,37.0,37.0}, {1600.0,66.0,66.0}};// Cn short reference: daily -- 900; hourly day -- 37 hourly night -- 37
                 // CN tall reference: daily 1600, hourly day 66; houry night
                 const double Cdi[2][3] = {{0.34, 0.24, 0.96}, {0.38, 0.25, 1.7}}; // Cd short reference: daily -- 0.34; hourly 0.24; hourly night 0.96
                 // Cd tall reference: dai;y 0.38, houry day 0.25, hiurly night 1.7
                 /**ASCE Penman-Monteith Terms Standardized for application of the RET equation*/
                 const double zveg[2] = {0.12,0.5}; // vegetation height short and tall reference
                 const double zh[2] = {2.0,2.0}; //[m], height of air temperature and humidity measurements 1.5-2.5 m
                 const double zw[2] = {2.0,2.0};//[m], height corresponding to wind speed
                 const double zo[2] = {0.08,0.08}; //[m], zero plane displacement height
                 const double lambda[2] = {2.45,2.45}; //[MJ/kg], latent heat of vaorization
                 const double rs_dailly[2] = {70.0, 45.0}; //[s/m], surface resistance, daily, short and tall reference
                 const double rs_h_day[2] = {50.0,30.0}; //[s/m], surface resistance, hourly during daytime, short and tall
                 const double rs_h_night[2] = {200.0,200.0};//[s/m], surface resistance, hourly during nighttime, short and tall
                 // predicting daytime Rnet > 0; nighttime Rnet <= 0.0
                 // ============================================ //
                 // constants for asce penman-monteith method
                 const double rho_w = 1000;//[kg/m^3]
                 const double ktime = 3600; // unit conversion for ET in mm/h
                 double h = 0.1; // mean height of vegetation --> parameter
                 double d(double h) { return 0.67 * h;} // zero displacement height, [m]
                 double zom(double h) {return 0.123*h;} // roughness length governing momentum transfer, [m]
                 double zoh(double h) {return 0.0123*h;}// roughness height for transfer of heat and vapor, [m]
                 double kappa = 0.41; // von Karman's constant

                 double resistnace_aerodynamic(double h, double zw, double uz){
                     return log((zw - d(h))/zom(h))*log((zoh(h)-d(h))/zoh(h))/kappa/kappa/uz;
                 }

                 double lai_a(double h){
                     double lai = 24.0*h; // for clipped grass
                     lai = 5.5 + 1.5*log(h);// for alfalfa /// TODO find condition to switch
                     return 0.5*lai;
                 }
                 double resistance_surface(double h){
                     double rl = 0.1; // effective stomatal resistnace of a well-illuminated leaf [s/m] TODO: move to parameter or find equation
                     return rl/lai_a(h);
                 }

                 double vaporization_latent_heat(double temperature){
                     return 2.501 - (2.361*0.001)*temperature;
                 }

                 /**\brief atmospheric air density
                  * \param pressure [kPa]
                  * \param temperature [degC]
                  * \param actual vapor pressure, [kPa]*/
                 double density_air(double pressure, double temperature, double ea){
                     double Tk = 273.16+temperature;
                     double Tkv = Tk/(1-0.378*ea/pressure);
                     return 3.486*pressure/Tkv;
                 }

             };


        }
    }

}
namespace shyft::test {



}


TEST_SUITE("evapotranspiration") {
    using shyft::core::radiation::parameter;
    using shyft::core::radiation::response;
    using shyft::core::radiation::calculator;
//    using shyft::core::radiation::surface_normal;
    using shyft::core::calendar;
    using shyft::core::utctime;
    // test basics: creation, etc



    TEST_CASE("hourly"){

            std::cout<<"hello, world"<<std::endl;

    }

}
