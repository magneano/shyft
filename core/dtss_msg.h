/** This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
See file COPYING for more details **/
#pragma once
#include <string>
#include <cstdint>
#include <exception>

namespace shyft {
namespace dtss {

/** \brief dtss message-types
 *
 * The message types used for the wire-communication of dtss.
 *
 */
enum class message_type : uint8_t {
	SERVER_EXCEPTION,
	EVALUATE_TS_VECTOR,
	EVALUATE_TS_VECTOR_PERCENTILES,
    FIND_TS,
    GET_TS_INFO,
	STORE_TS,
	CACHE_FLUSH,
	CACHE_STATS,
	EVALUATE_EXPRESSION,
	EVALUATE_EXPRESSION_PERCENTILES,
    MERGE_STORE_TS,
    REMOVE_TS,
	EVALUATE_TS_VECTOR_CLIP,
    EVALUATE_EXPRESSION_CLIP
	// EVALUATE_TS_VECTOR_HISTOGRAM //-- tsv,period,ta,bin_min,bin_max -> ts_vector[n_bins]
};

// ========================================

namespace msg {

/** stream utility functions for reading basic message-types/parts
 *
 */

template <class T>
message_type read_type(T& in) {
	int32_t mtype;
	in.read((char*)&mtype, sizeof(mtype));
	if (!in)
		throw dlib::socket_error(string("failed to read message type"));
	return (message_type)mtype;
}

template <class T>
void write_type(message_type mt, T& out) {
	int32_t mtype = (int32_t)mt;
	out.write((const char *)&mtype, sizeof(mtype));
	if (!out.good())
		throw dlib::socket_error(string("failed writing message type"));
}

template <class T>
void write_string(const std::string& s, T& out) {
	int32_t sz = s.size();
	out.write((const char*)&sz, sizeof(sz));
	out.write(s.data(), sz);
	if (!out.good())
		throw dlib::socket_error(string("failed writing string"));
}

template <class T>
std::string read_string(T& in) {
	std::int32_t sz;
	in.read((char*)&sz, sizeof(sz));
	if(!in)
		throw dlib::socket_error(string("failed reading size of string"));
	std::string msg(sz, '\0');
	in.read((char*)msg.data(), sz);
	if(!in)
		throw dlib::socket_error(string("failed reading string"));
	return msg;
}

template <class T>
void write_exception(const std::exception& e, T& out) {
	int32_t sz = strlen(e.what());
	out.write((const char*)&sz, sizeof(sz));
	out.write(e.what(), sz);
}

template <class T>
void send_exception(const std::exception& e, T& out) {
	write_type(message_type::SERVER_EXCEPTION, out);
	int32_t sz = strlen(e.what());
	out.write((const char*)&sz, sizeof(sz));
	out.write(e.what(), sz);
}

template <class T>
std::runtime_error read_exception(T& in) {
	int32_t sz;
	in.read((char*)&sz, sizeof(sz));
	if (!in)
		throw dlib::socket_error(string("failed reading exception size"));
	std::string msg(sz, '\0');
	in.read((char*)msg.data(), sz);
	if (!in)
		throw dlib::socket_error(string("failed reading exception data"));
	return std::runtime_error(msg);
}

}  // msg
}
}
