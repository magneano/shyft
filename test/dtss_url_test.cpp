#include "test_pch.h"
#include "core/dtss_url.h"
//#include <vector>
//#include <memory>
//#include <string>
#include <array>

using namespace std;

TEST_SUITE("dtss_url") {
    
TEST_CASE("shyft_url") {

    using shyft::dtss::shyft_url;

    FAST_CHECK_EQ( shyft_url("abc","123"), string("shyft://abc/123") );

    std::map<std::string, std::string> m1{ std::make_pair("foo", "bar"), std::make_pair("baz", "") };
    FAST_CHECK_EQ( shyft_url("abc","123", m1), string("shyft://abc/123?baz=&foo=bar") );

    // queries are encoded
    std::map<std::string, std::string> m2{ std::make_pair("$$$", "=_="), std::make_pair("baz", "   ") };
    FAST_CHECK_EQ( shyft_url("abc","123", m2), string("shyft://abc/123?%24%24%24=%3D_%3D&baz=+++") );
}
TEST_CASE("extract_shyft_url_container") {

    using shyft::dtss::extract_shyft_url_container;

    std::string extracted_1 = extract_shyft_url_container("shyft://abc/something/else");
    std::string extracted_2 = extract_shyft_url_container("shyft://abc/something/else?query=string&here=foo");
    FAST_CHECK_EQ( extracted_1, string("abc") );
    FAST_CHECK_EQ( extracted_2, string("abc") );

    std::string extracted_3 = extract_shyft_url_container("grugge");
    std::string extracted_4 = extract_shyft_url_container("grugge?query");
    FAST_CHECK_EQ( extracted_3, string{} );
    FAST_CHECK_EQ( extracted_4, string{} );

}
TEST_CASE("extract_shyft_url_path") {

    using shyft::dtss::extract_shyft_url_path;

    std::string extracted_1 = extract_shyft_url_path("shyft://abc/something/else");
    std::string extracted_2 = extract_shyft_url_path("shyft://abc/something/else?query=string&here=foo");
    FAST_CHECK_EQ( extracted_1, string("something/else") );
    FAST_CHECK_EQ( extracted_2, string("something/else") );

    std::string extracted_3 = extract_shyft_url_path("grugge://abc/something/else");
    std::string extracted_4 = extract_shyft_url_path("grugge?query");
    FAST_CHECK_EQ( extracted_3, string{} );
    FAST_CHECK_EQ( extracted_4, string{} );

}
TEST_CASE("extract_shyft_url_query") {

    using shyft::dtss::extract_shyft_url_query_parameters;

    auto m1 = extract_shyft_url_query_parameters("shyft://abc/something/else?query=string&here=");
    FAST_CHECK_EQ( m1.size(), 2 );
    FAST_REQUIRE_EQ( m1.count("query"), 1 );
    FAST_CHECK_EQ( m1["query"], std::string{"string"} );
    FAST_REQUIRE_EQ( m1.count("here"), 1 );
    FAST_CHECK_EQ( m1["here"], std::string{""} );

    auto m2 = extract_shyft_url_query_parameters("shyft://abc/something/else?query=string&here=foo");
    FAST_CHECK_EQ( m2.size(), 2 );
    FAST_REQUIRE_EQ( m2.count("query"), 1 );
    FAST_CHECK_EQ( m2["query"], std::string{"string"} );
    FAST_REQUIRE_EQ( m2.count("here"), 1 );
    FAST_CHECK_EQ( m2["here"], std::string{"foo"} );

    auto m3 = extract_shyft_url_query_parameters("grugge");
    FAST_CHECK_EQ( m3.size(), 0 );

    auto m4 = extract_shyft_url_query_parameters("grugge?query");
    FAST_CHECK_EQ( m4.size(), 0 );

    // queries are decoded
    auto m5 = extract_shyft_url_query_parameters("shyft://abc/123?%24%24%24=%3D_%3D&baz=+++");
    FAST_CHECK_EQ( m5.size(), 2 );
    FAST_REQUIRE_EQ( m5.count("$$$"), 1 );
    FAST_CHECK_EQ( m5["$$$"], std::string{"=_="} );
    FAST_REQUIRE_EQ( m5.count("baz"), 1 );
    FAST_CHECK_EQ( m5["baz"], std::string{"   "} );
}
TEST_CASE("remove_shyft_url_queries") {

    using shyft::dtss::remove_shyft_url_queries;

    auto url1 = remove_shyft_url_queries("shyft://abc/something/else?query=string&here=");
    FAST_CHECK_EQ( url1, "shyft://abc/something/else" );

    auto url2 = remove_shyft_url_queries("grugge");
    FAST_CHECK_EQ( url2, "" );

    auto url3 = remove_shyft_url_queries("grugge?query");
    FAST_CHECK_EQ( url3, "" );
}
TEST_CASE("filter_shyft_url_parsed_queries") {

    using shyft::dtss::filter_shyft_url_parsed_queries;

    std::array<std::string, 2> to_remove{{ "key02", "key06" }};
    std::map<std::string, std::string> queries{{
        {"key01", "value01" }, {"key02", "value02" },
        {"key03", "value03" }, {"key04", "value04" },
        {"key05", "value05" }, {"key06", "value06" }
    }};

    filter_shyft_url_parsed_queries(queries, to_remove);

    FAST_CHECK_EQ( queries.size(), 4 );
    FAST_CHECK_EQ( queries.find("key02"), queries.cend() );
    FAST_CHECK_EQ( queries.find("key06"), queries.cend() );
}
TEST_CASE("urlencode_urldecode") {

    using shyft::dtss::urlencode;
    using shyft::dtss::urldecode;

    // map of the reserved characters from RFC3986 (Jan. 2005)
    const std::map<std::string, std::string> reserved{
        { "!", "%21" }, { "#", "%23" }, { "$", "%24" },
        { "&", "%26" }, { "'", "%27" }, { "(", "%28" },
        { ")", "%29" }, { "*", "%2A" }, { "+", "%2B" },
        { ",", "%2C" }, { "/", "%2F" }, { ":", "%3A" },
        { ";", "%3B" }, { "=", "%3D" }, { "?", "%3F" },
        { "@", "%40" }, { "[", "%5B" }, { "]", "%5D" },
    };

    // all unreserved characters should not be encoded specially
    const std::string inp_01{ "abcdefghijklmnopqrstuvwxyz0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ-_.~" };
    const auto enc_01 = urlencode(inp_01);
    const auto dec_01 = urldecode(enc_01);
    FAST_CHECK_EQ( enc_01, inp_01 );
    FAST_CHECK_EQ( dec_01, inp_01 );

    // norwegian characters
    const std::string inp_02{ "æøåÆØÅ" };
    const auto enc_02 = urlencode(inp_02);
    const auto dec_02 = urldecode(enc_02);
    FAST_CHECK_EQ( dec_02, inp_02 );

    // explicitly stated reserved characters are encoded correctly
    for ( const auto & [ inp_03, encoded ] : reserved ) {
        const auto enc_03 = urlencode(inp_03);
        const auto dec_03 = urldecode(enc_03);
        FAST_CHECK_EQ( enc_03, encoded );
        FAST_CHECK_EQ( dec_03, inp_03 );
    }

    // + encoded SP character
    const std::string inp_04{ " " };
    const auto enc_04 = urlencode(inp_04);  // space_plus=true by default
    const auto dec_04 = urldecode(enc_04);  // space_plus=true by default
    FAST_CHECK_EQ( enc_04, "+" );
    FAST_CHECK_EQ( dec_04, inp_04 );

    // percent encoding of SP character
    const std::string inp_05{ " " };
    const auto enc_05 = urlencode(inp_05, false);
    const auto dec_05 = urldecode(enc_05, false);
    FAST_CHECK_EQ( enc_05, "%20" );
    FAST_CHECK_EQ( dec_05, inp_05 );

    // decode with mixed cased letters
    const std::string inp_06{ "%5b%28%2B%3d%2A%29%5d" };
    const std::string exp_06{ "[(+=*)]" };
    const auto dec_06 = urldecode(inp_06);
    FAST_CHECK_EQ( dec_06, exp_06 );
    
    // possible to decode encoded unreserved characters
    const std::string inp_07{ "%61%42%63%44%31%32%33" };
    const std::string exp_07{ "aBcD123" };
    const auto dec_07 = urldecode(inp_07);
    FAST_CHECK_EQ( dec_07, exp_07 );

    // percent itself is percent-encoded
    const std::string inp_08{ "%" };
    const std::string exp_08{ "%25" };
    const auto enc_08 = urlencode(inp_08);
    FAST_CHECK_EQ( enc_08, exp_08 );

    // error to decode unencoded reserved characters
    for ( const auto & [ inp_09, encoded ] : reserved ) {
        try {
            urldecode(inp_09, false);  // disable space as plus to test plus errors
        } catch ( const std::runtime_error & e ) {
            FAST_CHECK_EQ( std::string{e.what()}.find_first_of(std::string{"urldecode: "} + "invalid character '" + std::string{inp_09} + "'"), 0 );
            continue;
        }
        FAIL( std::string{ "No error when decoding reserved character '" } + std::string{inp_09} + "'" );
    }
    // verify it throws if invalid expansion at the end
    try {
        urldecode(string("abc%1"));
        FAST_CHECK_UNARY(false);
    } catch(const runtime_error& e) {
        FAST_CHECK_UNARY(true);
    }
    //verify it works if  %20 is the verylast thing in the string
    FAST_CHECK_EQ(urldecode(string("abc%20")),string("abc "));

}

}
