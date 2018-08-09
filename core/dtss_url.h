/** This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
See file COPYING for more details **/
#pragma once

#include <algorithm>
#include <iterator>
#include <map>
#include <limits>
#include <sstream>
#include <string>
#include <string_view>
#include <ctype.h>

namespace shyft {
namespace dtss {

// fwd
std::string urlencode(const std::string & text, bool space_plus=true);
std::string urldecode(const std::string & encoded, bool space_plus=true);


template < std::size_t N >
inline constexpr std::size_t char_str_length(const char(&a)[N]) noexcept {
    return N-1;  // don't count \0 char
}


// TODO: inline when vs implements P0386R2: Inline variables
constexpr char shyft_prefix[] = "shyft://";  ///< marks all internal handled urls


/** Construct a shyft-url from a container and a ts-name. */
inline std::string shyft_url(const std::string & container, const std::string & ts_name) {
    return std::string{ shyft_prefix } + container + "/" + ts_name;
}
/** Construct a shyft-url from a container, a ts-name, and a collection of query flags.
 *
 * Query keys and values are always urlencoded. The separating `?`, `&`, and `=` are not encoded.
 */
inline std::string shyft_url(const std::string & container, const std::string & ts_name, const std::map<std::string, std::string> & queries) {
    std::ostringstream str_s{ };
    if ( ! queries.empty() ) {
        str_s << "?";
        for ( auto it = queries.cbegin(); it != queries.cend(); ) {
            auto p = *it++;  // record current value, increment iterator
            str_s << urlencode(p.first) << '=' << urlencode(p.second) << (it != queries.cend() ? "&" : "");
        }
    }
    return shyft_url(container, ts_name) + str_s.str();
}


/** match & extract fast the following 'shyft://<container>/'
 * \param url like pattern above
 * \return <container> or empty string if no match
 */
inline std::string extract_shyft_url_container(const std::string & url) {
    if ( (url.size() < char_str_length(shyft_prefix) + 2) || !std::equal(std::begin(shyft_prefix), std::prev(std::end(shyft_prefix)), begin(url)) )
        return std::string{};
    // path after container?
    auto ce = url.find_first_of('/', char_str_length(shyft_prefix));  // container end
    if ( ce == std::string::npos )
        return std::string{};
    return url.substr(char_str_length(shyft_prefix), ce - char_str_length(shyft_prefix));
}


/** \brief Extract the path part from a shyft url.
 * \param  url  String url to extract the path from.
 * \return The string path if the url is a valid shyft url, else an empty string.
 */
inline std::string extract_shyft_url_path(const std::string & url) {
    if ( url.size() < (char_str_length(shyft_prefix) + 2) || ! std::equal(std::begin(shyft_prefix), std::prev(std::end(shyft_prefix)), begin(url)) )
        return std::string{};
    auto ce = url.find_first_of('/', char_str_length(shyft_prefix));  // container end
    if ( ce == url.npos )
        return std::string{};
    auto qs = url.find('?', ce);
    if ( qs == url.npos ) {
        return url.substr(ce + 1);
    } else {
        return url.substr(ce + 1, qs - ce - 1);
    }
}


/** Extract any query parameters from a url.
 *
 * The implementation ignores data until the first `?` character, afterwhich it parses the rest of
 * the url as a query string.
 *
 * The query string is assumed to be on the format `?key1=value1&key2=value2&key3=&key4=value4`.
 * This will be parsed into a map with four keys: `key1` through `key4`, where `key3` have a
 * empty string value, while the rest have respectivly values `value1`, `value2`, and `value4`.
 *
 * Both query keys and values are assumed to be urlencoded, thus urlencode is called on every
 * key and every value.
 *
 * If the url does not have any query parameters an empty mapping wil be returned.
 *
 * \param url String url to parse.
 * \return A map-type with the key-values from the url query string.
 */
inline std::map<std::string, std::string> extract_shyft_url_query_parameters(const std::string & url) {
    
    using map_t = std::map<std::string, std::string>;
    
    // locate query string if present
    auto it = std::find(url.cbegin(), url.cend(), '?');
    if ( it == url.cend() )
        return map_t{};
    std::advance(it, 1);  // skip to first character of the query string
    if ( it == url.cend() )
        return map_t{};

    // parse key=value pairs
    map_t queries{};
    std::string key_str{};
    std::string val_str{};
    while ( it != url.cend() ) {
        // locate key
        auto end = std::find(it, url.cend(), '=');
        if ( end == url.cend() ) {
            // no = after key, and at the end of the string -> do not add to map, then return
            break;
        }
        key_str.assign(urldecode({it, end}));  // found key!

        // locate value
        it = end + 1;  // skip past =
        end = std::find(it, url.cend(), '&');
        if ( std::distance(it, end) == 0 ) {
            // no value after = -> use empty string
            val_str.assign("");
        } else {
            // value after =
            val_str.assign(urldecode({it, end}));
        }

        // store and advance
        queries[key_str] = val_str;
        it = end;
        if ( it != url.cend() ) {
            std::advance(it, 1);
        }
    }
    return queries;
}


/** \brief Returns the shyft url without any query parameters.
 * \param  url  String url to remove queries from.
 * \return The string url without queries if the url is a valid shyft url, else an empty string.
 */
inline std::string remove_shyft_url_queries(std::string_view url) {
    if ( url.size() < (char_str_length(shyft_prefix) + 2) || ! std::equal(std::begin(shyft_prefix), std::prev(std::end(shyft_prefix)), begin(url)) )
        return std::string{};
    return std::string{ url.substr(0, url.find('?')) };
}


/** \brief Remove queries from the map of parsed queries.
 * \tparam  StrContainer  A container type supporting for-each style iterating,
 *      where the elements yielded are implicitly convertible to std::string.
 * \param  queries  Map of parsed queries.
 * \param  remove  The sequence of query keys to remove.
 */
template < typename StrContainer >
inline void filter_shyft_url_parsed_queries(std::map<std::string, std::string> & queries, const StrContainer & remove) {
    for ( const auto & str : remove ) {
        queries.erase(std::string{str});
    }
}

/** return true if supplied character should be unescaped according to RFC3986 */
inline bool is_unescaped_char(const char c) noexcept { return ('a' <= c && c <= 'z') || ('A' <= c && c <= 'Z') || ('0' <= c && c <= '9') || c == '-' || c == '_' || c == '.' || c == '~' ; }

/** \brief Percent-encode a string for use in URLs.
 *
 * All characters designated as reserved in RFC3986 (Jan. 2005), sec. 2.2 are always percent
 * encoded, while all character explitily stated as unreserved (RFC3986, Jan. 2005, sec. 2.3)
 * are never percent encoded. Space characters are encoded as `+` characters, while all other
 * characters are precent-encoded.
 *
 * The implementation only handles 8-bit character values. The behavior for multibyte characters
 * is unspecified.
 *
 * The reverse operation of this is `urldecode`.
 *
 * \param  text  Text string to encode.
 * \param  percent_plus  When true the `SP` character (ASCII 0x20) is encoded as `+` instead of
 *                       its percent encoding `%20`. Defaults to true.
 * \returns  Percent-encoded representation of the input.
 *
 * \sa  urldecode
 */

inline std::string urlencode(const std::string & text, bool space_plus /* = true */) {
    std::string encoded{};
    encoded.reserve(3*text.size());
    static const char *hex_c="0123456789ABCDEF";
    for ( const auto & c : text ) {
        if ( is_unescaped_char(c)) {
            encoded.push_back(c);
        } else if ( space_plus && c == ' ' ) {
            encoded.push_back('+');
        } else {
            unsigned char s = static_cast<unsigned char>(c);
            encoded.push_back('%');
            encoded.push_back(hex_c[(s&0xf0)>>4]);
            encoded.push_back(hex_c[s&0x0f]);
        }
    }
    return encoded;
}


/** \brief Decode a percent-encoded string to its original representation.
 *
 * All characters designated as unreserved in RFC3986 (Jan. 2005), sec. 2.3 are always passed
 * through unmodified, except where they are encountered while parsing a percent-encoded value.
 *
 * The implementation only handles 8-bit character values. The behavior for multibyte characters
 * is unspecified.
 * Additionally it is undefined if characters outside the range `0-9A-Fa-f` are encountered as
 * one of the two characters immidiatly succeeding a percent-sign.
 *
 * This is the reverse operation of `urlencode`.
 *
 * \param  encoded  Text string to decode.
 * \param  percent_plus  When true `+` characters are decoded to `SP` characters (ASCII 0x20).
 *                       When this is true and a `+` is encountered and exception is thrown.
 *                       The default value is true.
 *
 * \returns  Original representation of the encoded input.
 *
 * \throw  std::runtime_error  Thrown if unencoded characters outside the unreserved range is
 *                             encountered, includes `+` when the `percent_plus` argument is
 *                             false. The exception message contains the character and its
 *                             location in the string.
 *
 * \sa  urlencode
 */
inline std::string urldecode(const std::string & encoded, bool space_plus /* = true */) {
    std::string text{};
    text.reserve(encoded.size());
    for ( auto it = encoded.cbegin(); it != encoded.cend(); ++it ) {
        const auto& c = *it;
        if ( is_unescaped_char(c)) {
            text.push_back(c);
        } else if ( space_plus && c == '+' ) {
            text.push_back(' ');
        } else if ( c == '%' && ( encoded.cend()-it) >2) {
            char hex[3];hex[2]=0;++it;//skip percent
            hex[0]=*it++;hex[1]=*it;// leave iterator at the last digit so that end-loop increment does the job
            if(isxdigit(hex[0]) && isxdigit(hex[1])) {
                auto n = std::strtol( hex,nullptr,16);
                text.push_back(static_cast<char>(n));
            } else {
                throw std::runtime_error(std::string{"urldecode: "} + "invalid hex character '" + c + "' approx at position " + std::to_string(it-encoded.cbegin()));
            }
        } else {
            throw std::runtime_error(std::string{"urldecode: "} + "invalid character '" + c + "' at position " + std::to_string(it-encoded.cbegin()));
        }
    }

    return text;
}

}
}
