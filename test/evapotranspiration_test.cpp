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
#include "core/hydro_functions.h"

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
            using namespace shyft::core::hydro_functions;
            /**ASCE-EWRI hourly based model
             * it is stated taht this model can be used with the smaller periods as well, but we need to account for fractional hour
             * TODO: discuss how to manipulate between models and timesteps*/
             struct parameter{
                 double lai = 0.1; //leaf area index
                 double hveg = 0.1; // vegetation height, [m]
                 double rl = 0.1; // effective stomatal resistance, [s/m]
                 parameter(double hveg=0.1, double lai = 0.1,double rl=0.1) : hveg(hveg), lai(lai), rl(rl) {}
             };
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

                 void reference_evapotranspiration_asce(R& response, double temperature, double rhumidity, double elevation=0.1, double height_ws=2.0, double windspeed=0.0){

                     double net_radiation = 1.0;
                     double pressure = atm_pressure(elevation, 293.15);// recommended value for T0 during growing season is 20 gradC, see eq. B.8
                     double lambda_rho = vaporization_latent_heat(temperature)*rho_w;
                     double delta = svp_slope(temperature);
                     double avp = actual_vp(temperature, rhumidity);
                     double nominator = delta * (net_radiation + soil_heat_flux(net_radiation)) +
                             ktime*density_air(pressure, temperature, avp)*cp/resistance_aerodynamic(height_ws,windspeed)*(svp(temperature)-avp);
                     double denominator = delta + gamma_pm(pressure, temperature)*(1+resistance_surface(param.hveg)/resistance_aerodynamic(height_ws,windspeed));
                     response.et_ref = nominator/denominator/lambda_rho;
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
                 const double rho_w = 1;//[Mg/m^3]
                 const double cp = 1.013*0.001; // specific heat of moist air, [MJ/kg/gradC]
                 const double ktime = 3600; // unit conversion for ET in mm/h
                 //double h = 0.1; // mean height of vegetation --> parameter
                 double d =  0.67 * param.hveg; // zero displacement height, [m]
                 double zom = 0.123*param.hveg; // roughness length governing momentum transfer, [m]
                 double zoh =  0.0123*param.hveg;// roughness height for transfer of heat and vapor, [m]
                 double kappa = 0.41; // von Karman's constant

                 /**\brief aerodynamic resistance, [s/m], eq. B.2
                  * \param hws -- wind speed measurements height, [m]
                  * \param ws -- wind speed at height hws, [m/s]*/
                 double resistance_aerodynamic(double hws, double ws){
                     return log((hws - d)/zom)*log((zoh-d)/zoh)/kappa/kappa/ws;
                 }

                 /**\brief active leaf-area index (LAI), eq.B.4
                  * \param hveg -- vegetation height, [m]*/
                 double lai_a(double hveg){
                     double lai = 24.0*hveg; // for clipped grass
                     lai = 5.5 + 1.5*log(hveg);// for alfalfa /// TODO find condition to switch
                     return 0.5*lai;
                 }
                 /**\brief bulk surface resistance, [s/m], eq. B.3
                  * \param zveg -- vegetation height,[m]*/
                 double resistance_surface(double hveg=0.1){
                     return param.rl/lai_a(hveg);
                 }

                 /// TODO: move all to the hydro_functions???

                 /**\brief Latent heat of vaporization, eq. B.7; ASCE-EWRI*/
                 double vaporization_latent_heat(double temperature){
                     return 2.501 - (2.361*0.001)*temperature;
                 }
//                 /**\brief atmospheric pressure, ASCE-EWRI, eq.B.8
//                  * */
//                 double pressure_air(double elevation){
//                     const double p0 = 101.325;// pressure at sea level, [kPa]
//                     const double T0 = 293.15; // temperature at sea level for growing season, [K]
//                     const double alpha1 = 0.0065; // lapse rate, [K/m]
//                     const double z0 = 0; // sea level, [m]
//                     const double Rgas = 287; // universal gas constant, [J/kg/K]
//                     const double g = 9.807; // gravitational acceleration, [m/s2]
//                     return p0*pow((T0-alpha1*(elevation-z0))/T0,g/alpha1/Rgas);
//                 }

                 /**\brief atmospheric air density, ASCE-EWRI, eq. B.10
                  * \param pressure [kPa]
                  * \param temperature [degC]
                  * \param actual vapor pressure, [kPa]*/
                 double density_air(double pressure, double temperature, double ea){
                     double Tk = 273.16+temperature;
                     double Tkv = Tk/(1-0.378*ea/pressure); // B.11
                     return 3.486*pressure/Tkv; // B.10
                 }

                 /**\brief psychrometric constant, [kPa/gradC]
                  * \param atmospheric pressure, [kPa]
                  * \param air temperature, [degC]*/
                 double gamma_pm(double pressure, double temperature){

                     double lambda = vaporization_latent_heat(temperature);
                     const double epsilon = 0.622; // ratio of molecular weight of water vapor/dry air
                     return cp*pressure/epsilon/lambda;
                 }

                 /**\brief soil heat flux density (G) for hourly periods, eq.B.13*/
                 double soil_heat_flux(double net_radiation){
                     double kg = 0.4;
                     if (net_radiation <= 0.0)
                         kg = 2.0;
                     return kg*exp(-0.5*lai_a(param.zveg))*net_radiation;
                 }

                 /**\brief wind speed adjustment for measurement height, eq.B.14*/
                 double ws_adjustment(double height_measure, double ws_measure){
                     return ws_measure*log((2-d)/zom)/log((height_measure-d)/zom);
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
