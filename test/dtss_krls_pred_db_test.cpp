#include "test_pch.h"

#include "core/dtss_krls.h"

#include <boost/filesystem.hpp>
namespace  fs=boost::filesystem;


TEST_SUITE("dtss_krls_pred_db") {

TEST_CASE("dtss_krls_io_routines") {
    using krls_io = shyft::dtss::krls_pred_db_io;

    shyft::core::calendar utc;

    fs::path datapath = fs::temp_directory_path() / fs::unique_path();
    fs::create_directories(datapath);
    const std::string ffp = (datapath / fs::path{ "krls_datafile.b" }).string();
    std::unique_ptr<std::FILE, decltype(&std::fclose)> fh{ std::fopen(ffp.c_str(), "w+b"), &std::fclose };

    const std::string source_url = "shyft://container/path/to/ts";
    const auto period = shyft::core::utcperiod{ utc.time(2017, 1, 1), utc.time(2018, 1, 1) };
    const auto krls_dt = utc.DAY;
    const auto point_fx = shyft::time_series::ts_point_fx::POINT_INSTANT_VALUE;
    const std::size_t krls_dict_size = 100000;
    const double krls_tolerance = 0.0001;
    const double rbfk_gamma = 0.001;

    // create the default kernel
    const auto predictor = krls_io::create_rbf_file(fh.get(), source_url, period, krls_dt, point_fx, krls_dict_size, krls_tolerance, rbfk_gamma);

    // attempt to use io functions to read the file back, all with default skip

    FAST_CHECK_UNARY(krls_io::can_read_file(fh.get()));

    auto read_krls_header = krls_io::read_header(fh.get());
    FAST_CHECK_EQ(read_krls_header.scaling, shyft::core::to_seconds64(krls_dt) );
    FAST_CHECK_EQ(read_krls_header.tolerance, krls_tolerance);
    FAST_CHECK_EQ(read_krls_header.point_fx, point_fx);
    FAST_CHECK_EQ(read_krls_header.t_start, shyft::core::to_seconds64(period.start));
    FAST_CHECK_EQ(read_krls_header.t_end, shyft::core::to_seconds64(period.end));

    auto read_source_url = krls_io::read_source_url(fh.get());
    FAST_CHECK_EQ(read_source_url, source_url);

    auto read_kernel_type = krls_io::read_predictor_kernel_type(fh.get());
    FAST_CHECK_EQ(read_kernel_type, shyft::dtss::krls_kernel_type_identifiers::radial_basis_kernel);

    auto read_rbf_header = krls_io::read_predictor_rbf_header(fh.get());
    FAST_CHECK_EQ(read_rbf_header.gamma, rbfk_gamma);

    auto read_predictor = krls_io::read_predictor_rbf_predictor(fh.get());
    FAST_CHECK_EQ(read_predictor.get_predictor_dt(), predictor.get_predictor_dt());
    FAST_CHECK_EQ(read_predictor.get_tolerance(), predictor.get_tolerance());
    FAST_CHECK_EQ(read_predictor.get_dictionary_size(), predictor.get_dictionary_size());
    FAST_CHECK_EQ(read_predictor.get_max_dictionary_size(), predictor.get_max_dictionary_size());
    FAST_CHECK_EQ(read_predictor.get_rbf_gamma(), predictor.get_rbf_gamma());
    FAST_CHECK_EQ(read_predictor.get_predicted_ts_point_policy(), predictor.get_predicted_ts_point_policy());
}

TEST_CASE("dtss_krls_db_register_update_read_rbf_series") {
    using calendar = shyft::core::calendar;
    using krls_db_t = shyft::dtss::krls_pred_db;
    using ts_vector_t = shyft::time_series::dd::ats_vector;
    using utcperiod = shyft::core::utcperiod;

    calendar utc;

    fs::path krls_server_path = fs::temp_directory_path() / fs::unique_path();

    const std::string krls_name = "test-series";
    const std::string source_url = "test-series-url";
    const auto earliest_time = utc.time(2017, 1, 1);
    const auto middle_time = utc.time(2017, 3, 1);
    const auto latest_time = utc.time(2017, 5, 1);
    auto period = utcperiod{ earliest_time, middle_time };
    const auto krls_dt = 6*calendar::HOUR;
    const auto point_fx = shyft::time_series::ts_point_fx::POINT_INSTANT_VALUE;
    const std::size_t krls_dict_size = 1000000;
    const double krls_tolerance = 0.001;
    const double rbfk_gamma = 0.001;

    // server callback
    int data_dt = 2;
    int callback_call_count = 0;
    double data_offset = 0;
    std::function<ts_vector_t(const std::string &, utcperiod, bool, bool)> server_read_cb =
        [&callback_call_count, &period, &data_offset, point_fx, data_dt](const std::string &, utcperiod, bool, bool) -> ts_vector_t
    {
        callback_call_count += 1;

        // construct time axis
        std::size_t n = (period.end - period.start)/(data_dt*calendar::HOUR);
        shyft::time_axis::fixed_dt ta{ period.start, data_dt*calendar::HOUR, n };
        // generate values
        std::vector<double> values{};
        values.reserve(n);
        for ( std::size_t i = 0; i < n; ++i ) {
            values.emplace_back(data_offset + callback_call_count*i);
        }
        // setup return
        ts_vector_t result{};
        result.reserve(1);
        result.emplace_back(ta, values, point_fx);
        return result;
    };

    // setup the krls container
    krls_db_t krls_db{ krls_server_path.string(), server_read_cb };

    // register a new series
    krls_db.register_rbf_series(krls_name, source_url, period, krls_dt, point_fx, krls_dict_size, krls_tolerance, rbfk_gamma);
    // --------------------
    FAST_CHECK_EQ( callback_call_count, 1 );
    FAST_CHECK_UNARY( fs::exists(krls_server_path / krls_name) );

    // predict using the written data
    std::size_t n = (period.end - period.start)/calendar::HOUR;
    shyft::time_axis::generic_dt predict_ta{ period.start, calendar::HOUR, n };
    auto result_ts = krls_db.predict_time_series(krls_name, predict_ta);
    // --------------------
    FAST_CHECK_EQ( callback_call_count, 1 );
    FAST_REQUIRE_EQ( n, result_ts.size() );
    for ( std::size_t i = 10; i < n; ++i ) {  // skip 10 samples because of initial error
        FAST_CHECK_EQ( result_ts.value(i), doctest::Approx(callback_call_count*i/double(data_dt)).epsilon(0.1).scale(1) );  // epsilon is percent difference
    }
    data_offset = result_ts.value(n-1);

    // update the data with one more year
    period = utcperiod{ middle_time, latest_time };
    krls_db.update_rbf_series(krls_name, period);
    FAST_CHECK_EQ( callback_call_count, 2 );
    // --------------------
    n = (period.end - period.start)/calendar::HOUR;
    predict_ta = shyft::time_axis::generic_dt{ period.start, calendar::HOUR, n };
    result_ts = krls_db.predict_time_series(krls_name, predict_ta);
    // --------------------
    FAST_REQUIRE_EQ( n, result_ts.size() );
    for ( std::size_t i = 10; i < n; ++i ) {  // skip 10 samples because of initial error
        FAST_CHECK_EQ( result_ts.value(i), doctest::Approx(data_offset + callback_call_count*i/double(data_dt)).epsilon(0.1).scale(1) );  // epsilon is percent difference
    }
}

TEST_CASE("dtss_krls_db_register_update_read_rbf_series_using_save_read") {
    using calendar = shyft::core::calendar;
    using gta_t = shyft::time_axis::generic_dt;
    using gts_t = shyft::time_series::point_ts<gta_t>;
    using krls_db_t = shyft::dtss::krls_pred_db;
    using ts_vector_t = shyft::time_series::dd::ats_vector;
    using utcperiod = shyft::core::utcperiod;


    calendar utc;

    fs::path krls_server_path = fs::temp_directory_path() / fs::unique_path();

    const std::string krls_name = "test-series";
    const std::string source_url = "test-series-url";
    const auto earliest_time = utc.time(2017, 1, 1);
    const auto middle_time = utc.time(2017, 3, 1);
    const auto latest_time = utc.time(2017, 5, 1);
    auto period = utcperiod{ earliest_time, middle_time };
    const auto krls_dt = 6*calendar::HOUR;
    const auto point_fx = shyft::time_series::ts_point_fx::POINT_INSTANT_VALUE;
    const std::size_t krls_dict_size = 1000000;
    const double krls_tolerance = 0.001;
    const double rbfk_gamma = 0.001;

    // server callback
    int data_dt = 2;
    int callback_call_count = 0;
    double data_offset = 0;
    std::function<ts_vector_t(const std::string &, utcperiod, bool, bool)> server_read_cb =
        [&callback_call_count, &period, &data_offset, point_fx, data_dt](const std::string &, utcperiod, bool, bool) -> ts_vector_t
    {
        callback_call_count += 1;

        // construct time axis
        std::size_t n = (period.end - period.start)/(data_dt*calendar::HOUR);
        gta_t ta{ period.start, data_dt*calendar::HOUR, n };
        // generate values
        std::vector<double> values{};
        values.reserve(n);
        for ( std::size_t i = 0; i < n; ++i ) {
            values.emplace_back(data_offset + callback_call_count*i);
        }
        // setup return
        ts_vector_t result{};
        result.reserve(1);
        result.emplace_back(ta, values, point_fx);
        return result;
    };

    // setup the krls container
    krls_db_t krls_db{ krls_server_path.string(), server_read_cb };

    // register a new series
    std::size_t n = (period.end - period.start)/calendar::HOUR;
    krls_db.save(krls_name,
        gts_t{ gta_t{ period.start, calendar::HOUR, n }, 0., point_fx }, false,
        krls_db_t::queries_t{
            {"source_url", source_url},
        {"dt_scaling", std::to_string(shyft::core::to_seconds64(krls_dt))},
        {"point_fx", "instant"},
        {"krls_dict_size", std::to_string(krls_dict_size)},
        {"tolerance", std::to_string(krls_tolerance)},
        {"gamma", std::to_string(rbfk_gamma)}
        });
    // --------------------
    FAST_CHECK_EQ( callback_call_count, 1 );
    FAST_CHECK_UNARY( fs::exists(krls_server_path / krls_name) );

    // predict using the written data
    auto result_ts = krls_db.read(krls_name, period, krls_db_t::queries_t{
        {"dt", std::to_string(shyft::core::to_seconds64(calendar::HOUR))},
        });
    // --------------------
    FAST_CHECK_EQ( callback_call_count, 1 );
    FAST_REQUIRE_EQ( n, result_ts.size() );
    for ( std::size_t i = 10; i < n; ++i ) {  // skip 10 samples because of initial error
        FAST_CHECK_EQ( result_ts.value(i), doctest::Approx(callback_call_count*i/double(data_dt)).epsilon(0.1).scale(1) );  // epsilon is percent difference
    }
    data_offset = result_ts.value(n-1);

    // update the data with one more year
    period = utcperiod{ middle_time, latest_time };
    n = (period.end - period.start)/calendar::HOUR;
    krls_db.save(krls_name, gts_t{ gta_t{ period.start, calendar::HOUR, n }, 0., point_fx },
        false, krls_db_t::queries_t{});
    FAST_CHECK_EQ( callback_call_count, 2 );
    // --------------------
    n = (period.end - period.start)/calendar::HOUR;
    result_ts = krls_db.read(krls_name, period, krls_db_t::queries_t{
        {"dt", std::to_string(shyft::core::to_seconds64(calendar::HOUR))},
        });
    // --------------------
    FAST_REQUIRE_EQ( n, result_ts.size() );
    for ( std::size_t i = 10; i < n; ++i ) {  // skip 10 samples because of initial error
        FAST_CHECK_EQ( result_ts.value(i), doctest::Approx(data_offset + callback_call_count*i/double(data_dt)).epsilon(0.1).scale(1) );  // epsilon is percent difference
    }
}

TEST_CASE("dtss_krls_pred_db_move") {
    using calendar = shyft::core::calendar;
    using gta_t = shyft::time_axis::generic_dt;
    using gts_t = shyft::time_series::point_ts<gta_t>;
    using krls_db_t = shyft::dtss::krls_pred_db;
    using ts_vector_t = shyft::time_series::dd::ats_vector;
    using utcperiod = shyft::core::utcperiod;
    
    
    calendar utc;
    
    fs::path krls_server_path = fs::temp_directory_path() / fs::unique_path();
    
    const std::string krls_name = "test-series";
    const std::string move_to_name = "destination-test-series";
    const std::string source_url = "test-series-url";
    const auto earliest_time = utc.time(2017, 1, 1);
    const auto middle_time = utc.time(2017, 3, 1);
    auto period = utcperiod{ earliest_time, middle_time };
    const auto krls_dt = 6*calendar::HOUR;
    const auto point_fx = shyft::time_series::ts_point_fx::POINT_INSTANT_VALUE;
    const std::size_t krls_dict_size = 1000000;
    const double krls_tolerance = 0.001;
    const double rbfk_gamma = 0.001;
    
    // server callback
    int data_dt = 2;
    int callback_call_count = 0;
    double data_offset = 0;
    std::function<ts_vector_t(const std::string &, utcperiod, bool, bool)> server_read_cb =
    [&callback_call_count, &period, &data_offset, point_fx, data_dt](const std::string &, utcperiod, bool, bool) -> ts_vector_t
    {
        callback_call_count += 1;
    
        // construct time axis
        std::size_t n = (period.end - period.start)/(data_dt*calendar::HOUR);
        gta_t ta{ period.start, data_dt*calendar::HOUR, n };
        // generate values
        std::vector<double> values{};
        values.reserve(n);
        for ( std::size_t i = 0; i < n; ++i ) {
            values.emplace_back(data_offset + callback_call_count*i);
        }
        // setup return
        ts_vector_t result{};
        result.reserve(1);
        result.emplace_back(ta, values, point_fx);
        return result;
    };
    
    // setup the krls container
    krls_db_t krls_db{ krls_server_path.string(), server_read_cb };
    
    // register a new series
    std::size_t n = (period.end - period.start)/calendar::HOUR;
    krls_db.save(krls_name,
        gts_t{ gta_t{ period.start, calendar::HOUR, n }, 0., point_fx }, false,
        krls_db_t::queries_t{
            {"source_url", source_url},
            {"dt_scaling", std::to_string(shyft::core::to_seconds64(krls_dt))},
            {"point_fx", "instant"},
            {"krls_dict_size", std::to_string(krls_dict_size)},
            {"tolerance", std::to_string(krls_tolerance)},
            {"gamma", std::to_string(rbfk_gamma)}
        });
    // --------------------
    FAST_CHECK_EQ( callback_call_count, 1 );
    FAST_CHECK_UNARY( fs::exists(krls_server_path / krls_name) );
    
    krls_db.save(krls_name,
        gts_t{ gta_t{ period.start, calendar::HOUR, n }, 0., point_fx }, false,
        krls_db_t::queries_t{
            {"destination", move_to_name}
        }
    );
    // --------------------
    FAST_CHECK_EQ( callback_call_count, 1 );  // did not call again
    FAST_CHECK_UNARY( fs::exists(krls_server_path / move_to_name) );  // new name exist
    FAST_CHECK_UNARY_FALSE( fs::exists(krls_server_path / krls_name) );  // old name does not exist
}

}
