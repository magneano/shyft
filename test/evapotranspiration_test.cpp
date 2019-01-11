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
                 double lai = 2.0; //leaf area index, defaulted value taken from MORECS-metoff 1981 model for grass
                 double hveg = 0.15; // vegetation height, [m] (grass from MORECS)
                 double rl = 50.0; // effective stomatal resistance, [s/m] (calculated for grass from MORECS data)
                 parameter(double hveg=0.12, double lai = 2.0,double rl=50.0) : hveg(hveg), lai(lai), rl(rl) {}
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

                /**\brief computes ASCE PENMAN-MONTEITH reference evapotranspiation, [mm/h] based on eq. B.1 from ASCE-EWRI
                 * \param response -- updating response structure
                 * \param net_radiation, [MJ/m2/h] -- net radiation
                 * \param temperature [degC] -- air temperature
                 * \param rhumidity [-] -- relative humidity
                 * \param elevation [m]
                 * \param height_ws [m] -- wind speed measurement height
                 * \param windspeed [m/s at height_ws]
                 * \param height_t [m] -- height of humidity,temperature measuerements*/
                 void reference_evapotranspiration_asce(R& response, double net_radiation, double temperature, double rhumidity, double elevation=0.1, double height_ws=2.0, double windspeed=0.1, double height_t=2.0){

                     double pressure = atm_pressure(elevation, 293.15);// recommended value for T0 during growing season is 20 gradC, see eq. B.8
                     double lambda_rho = vaporization_latent_heat(temperature)*rho_w;
//                     std::cout<<"lambda_rho: "<<lambda_rho<<std::endl;
                     double delta = svp_slope(temperature);
//                     std::cout<<"delta: "<<delta<<std::endl;
                     double avp = actual_vp(temperature, rhumidity);
//                     std::cout<<"avp: "<<avp<<std::endl;
                     double G = soil_heat_flux(net_radiation);
//                     std::cout<<"G: "<<G<<std::endl;
                     double rho = density_air(pressure, temperature, avp);
//                     std::cout<<"rho: "<<rho<<std::endl;
                     double ra = resistance_aerodynamic(height_ws,windspeed);
//                     std::cout<<"ra: "<<ra<<std::endl;
                     double sat_vp = svp(temperature);
//                     std::cout<<"sat_vp: "<<sat_vp<<std::endl;
                     double nominator = delta * (net_radiation + soil_heat_flux(net_radiation)) +
                             ktime*density_air(pressure, temperature, avp)*cp/resistance_aerodynamic(height_ws,windspeed, height_t)*(svp(temperature)-avp);
//                     std::cout<<"nominator: "<<nominator<<std::endl;
                     double denominator = delta + gamma_pm(pressure, temperature)*(1+resistance_surface(param.hveg)/resistance_aerodynamic(height_ws,windspeed));
//                     std::cout<<"denominator: "<<denominator<<std::endl;
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
                  * \param ws -- wind speed at height hws, [m/s]
                  * \param ht [m] -- height of humidity and or temperature measurements*/
                 double resistance_aerodynamic(double hws, double ws, double ht=2.0){
                     return log((hws - d)/zom)*log((ht-d)/zoh)/kappa/kappa/min(ws,0.01);
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
                 double resistance_surface(double hveg=0.1){ /// TODO discuss the incorporation of MORECS model here, is it worth to account for seasons?
                     return param.rl/lai_a(hveg);
                 }

                 /// TODO: move all to the hydro_functions???

                 /**\brief Latent heat of vaporization, eq. B.7; ASCE-EWRI
                  * \param temperature [degC]*/
                 double vaporization_latent_heat(double temperature = 20.0){
                     return 2.501 - (2.361*0.001)*temperature;//with default value gives 2.45;
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
                 double gamma_pm(double pressure, double temperature=20.0){

                     double lambda = vaporization_latent_heat(temperature);
                     const double epsilon = 0.622; // ratio of molecular weight of water vapor/dry air
                     return cp*pressure/epsilon/lambda;
                 }

                 /**\brief soil heat flux density (G) for hourly periods, eq.B.13
                  * [net_radiation] should be either MJ/m2/h or Wt/m2*/
                 double soil_heat_flux(double net_radiation){
                     double kg = 0.4;
                     if (net_radiation <= 0.0)
                         kg = 2.0;
                     return kg*exp(-0.5*lai_a(param.hveg))*net_radiation;
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

    class trapezoidal_average {
    private:
        double area = 0.0;
        double f_a = 0.0;; // Left hand side of next integration subinterval
        double t_start = 0.0; // Start of integration period
        double t_a = 0.0; // Left hand side time of next integration subinterval
    public:
        explicit trapezoidal_average() {}

        /** \brief initialize must be called to reset states before being used during ode integration.
         */
        void initialize(double f0, double t_start) {
            this->f_a = f0;
            this->t_start = t_start;
            t_a = t_start;
            area = 0.0;
        }

        /** \brief Add contribution to average using a simple trapezoidal rule
         *
         * See: http://en.wikipedia.org/wiki/Numerical_integration
         */
        void add(double f, double t) {
            area += 0.5*(f_a + f)*(t - t_a);
            f_a = f;
            t_a = t;
        }

        double result() const { return area/(t_a - t_start); }
    };

}


TEST_SUITE("evapotranspiration") {
    using namespace shyft::core;
//    using shyft::core::radiation::surface_normal;
    using shyft::core::calendar;
    using shyft::core::utctime;
    using namespace shyft::test;
    // test basics: creation, etc



    TEST_CASE("hourly"){

        //===========================//
        // getting radiation  //
        radiation::parameter rad_p;
        radiation::response rad_r;
        rad_p.albedo = 0.2;
        rad_p.turbidity = 1.0;
        radiation::calculator<radiation::parameter,radiation::response> rad(rad_p);
        calendar utc_cal;
        // Greeley, Colorado weather station
        double lat = 40.41;
        double elevation = 1462.4;
        double ht = 1.68;
        double hws = 3.0;
        utctime t;
        // checking for horizontal surface Eugene, OR, p.64, fig.1b
        arma::vec surface_normal({0.0,0.0,1.0});
        double slope = 0.0;
        double aspect = 0.0;
        utctime ta;

        ta = utc_cal.time(2000, 06, 2, 00, 00, 0, 0);
        //rad.psw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0);




        //evapotranspiraiton PM
        evapotranspiration::penman_monteith::parameter pm_p;
        evapotranspiration::penman_monteith::response pm_r;
        evapotranspiration::penman_monteith::calculator<evapotranspiration::penman_monteith::parameter,evapotranspiration::penman_monteith::response> pm_calculator(pm_p);
        trapezoidal_average av_et_ref;
        av_et_ref.initialize(pm_r.et_ref, 0.0);
        // ref.: ASCE=EWRI Appendix C: Example Calculation of ET
        double temperature[23] = {16.5, 15.4, 15.5, 13.5, 13.2, 16.2, 20.0, 22.9, 26.4, 28.2, 29.8, 30.9, 31.8, 32.5, 32.9, 32.4, 30.2, 30.6, 28.3, 25.9, 23.9, 20.1, 19.9};
        double vap_pressure[23] = {1.26, 1.34, 1.31, 1.26, 1.24, 1.31, 1.36, 1.39, 1.25, 1.17, 1.03, 1.02, 0.98, 0.87, 0.86, 0.93, 1.14, 1.27, 1.27, 1.17, 1.20, 1.10, 1.05};
        double windspeed[23] = {0.5, 1.0, 0.68, 0.69, 0.29, 1.24, 1.28, 0.88, 0.72, 1.52, 1.97, 2.07, 2.76, 2.990, 3.10, 2.77, 3.41, 2.78, 2.95, 3.27, 2.86, 2.7, 2.0};
        double svp[23];
        double rhumidity[23];
        for (int i=0;i<23;i++){
            svp[i] = hydro_functions::svp(temperature[i]);
            //std::cout<<"svp: "<<svp[i]<<std::endl; // matched
            rhumidity[i] = svp[i]*100/vap_pressure[i];
        }

        rad.net_radiation(rad_r, lat, ta, slope, aspect, temperature[0], rhumidity[0], elevation);
        for (int h = 1; h < 24; ++h) {
            t = utc_cal.time(2000, 06, 2, h, 00, 0, 0); // June
            rad.net_radiation(rad_r, lat, t, slope,aspect, temperature[h-1], rhumidity[h-1], elevation);
            std::cout<<"lw: "<<rad_r.lw_radiation<<std::endl;
            std::cout<<"sw: "<<rad_r.sw_radiation<<std::endl;
            std::cout<<"net: "<<rad_r.net_radiation<<std::endl;
            std::cout<<"-----------------"<<std::endl;
            pm_calculator.reference_evapotranspiration_asce(pm_r,rad_r.net_radiation*0.0036,temperature[h-1], rhumidity[h-1],elevation, hws, windspeed[h-1],ht);
            std::cout<<"et_ref: "<<pm_r.et_ref<<std::endl;
            std::cout<<"-----------------"<<std::endl;
            av_et_ref.add(pm_r.et_ref, h);

        }
        std::cout << "et_ref_soh_av: " << av_et_ref.result() << std::endl;
        FAST_CHECK_EQ(av_et_ref.result(), doctest::Approx(0.4).epsilon(0.005));

    }

}
