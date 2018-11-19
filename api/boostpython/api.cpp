/** This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
See file COPYING for more details **/
#include "boostpython_pch.h"
#include <fstream>

char const* version() {
   return "v4.x";
}
namespace expose {
    extern void api_geo_point();
    extern void api_geo_cell_data();
    extern void calendar_and_time();
    extern void vectors();
    extern void api_time_axis();
    extern void timeseries();
    extern void target_specification();
    extern void region_environment() ;
    extern void priestley_taylor();
    extern void actual_evapotranspiration();
    extern void gamma_snow();
    extern void universal_snow();
    extern void kirchner();
    extern void precipitation_correction();
    extern void hbv_snow();
    extern void hbv_physical_snow();
    extern void cell_environment();
    extern void interpolation();
    extern void skaugen_snow();
    extern void kalman();
	extern void hbv_soil();
	extern void hbv_tank();
	extern void hbv_actual_evapotranspiration();
	extern void glacier_melt();
	extern void routing();
	extern void dtss();
    extern void dtss_finalize();
    extern void api_cell_state_id();
    extern void radiation();


    static std::vector<char> byte_vector_from_file(std::string path) {
        using namespace std;
        ostringstream buf;
        ifstream input;
        input.open(path.c_str(),ios::in|ios::binary);
        if (input.is_open()) {
            buf << input.rdbuf();
            auto s = buf.str();
            return std::vector<char>(begin(s), end(s));
        } else
            throw runtime_error(string("failed to open file for read:") + path);
    }
    static std::vector<char> byte_vector_from_hex_str(const std::string& s) {
        using namespace std;
        vector<char> r;
        if(s.size()==0)
            return r;
        if(s.size() %2 )
            throw runtime_error(string("hex_str should be even-sized"));
        r.reserve(s.size()/2);
        for(size_t i=0;i<s.size();i+=2) {
            uint32_t b=0;
            if(sscanf(s.c_str()+i,"%02x",&b)==1) {
                r.push_back(b&0x00ffu);
            } else {
                throw runtime_error(string("illegal hex at ")+to_string(i)+ string(" offset "));
            }
        }
        return r;
    }
    static std::string byte_vector_to_hex_str(const std::vector<char>&b) {
        using namespace std;
        string r;
        r.reserve(b.size()*2);
        for(const unsigned char c:b){
            char s[10];
            sprintf(s,"%02x",c); // c
            r.push_back(s[0]);
            r.push_back(s[1]);
        }
        return r;
    }

    static void byte_vector_to_file(std::string path, const std::vector<char>&bytes) {
        using namespace std;
        ofstream out;
        out.open(path, ios::out | ios::binary | ios::trunc);
        if (out.is_open()) {
            out.write(bytes.data(), bytes.size());
            out.flush();
            out.close();
        } else {
            throw runtime_error(string("failed to open file for write:") + path);
        }
    }

    void api() {
        calendar_and_time();
        vectors();
        api_time_axis();
        api_geo_point();
        api_geo_cell_data();
        timeseries();
        target_specification();
        region_environment();
        precipitation_correction();
        priestley_taylor();
        actual_evapotranspiration();
        gamma_snow();
        skaugen_snow();
        hbv_snow();
        hbv_physical_snow();
        kirchner();
        cell_environment();
        interpolation();
        kalman();
		hbv_soil();
		hbv_tank();
		hbv_actual_evapotranspiration();
		glacier_melt();
		routing();
        api_cell_state_id();
        radiation();
        using namespace boost::python;
        def("byte_vector_from_file", byte_vector_from_file, (arg("path")), "reads specified file and returns its contents as a ByteVector");
        def("byte_vector_to_file", byte_vector_to_file, (arg("path"), arg("byte_vector")), "write the supplied ByteVector to file as specified by path");
        def("byte_vector_from_hex_str", byte_vector_from_hex_str, (arg("hex_str")), "converts earlier byte-vector hex-string to byte-vector");
        def("byte_vector_to_hex_str", byte_vector_to_hex_str, (arg("byte_vector")), "return hex-string of byte-vector");
        dtss();
    }
}
void finalize_api() {
    //extern void expose::dtss_finalize();
    expose::dtss_finalize();
}
#ifdef _WIN32
// there are some extra fuzz needed to get win working
// it goes here:

#define _WIN32_WINNT 0x0501 // XP and above
#include <sdkddkver.h>
#include <stdio.h>
#include <windows.h>
#pragma comment(lib,"ntdll.lib")
EXTERN_C NTSTATUS NTAPI NtSetInformationProcess(HANDLE, ULONG, PVOID, ULONG);
#if 0
typedef NTSTATUS(NTAPI *NtQueryInformationProcessFn)(HANDLE process, ULONG infoClass, void* data,ULONG dataSize, ULONG* outSize);
static NtQueryInformationProcessFn NtQueryInformationProcess;

typedef NTSTATUS(NTAPI *NtSetInformationProcessFn)(HANDLE process, ULONG infoClass, void* data,ULONG dataSize);
static NtSetInformationProcessFn NtSetInformationProcess;
#endif
// these values determined by poking around in the debugger - use at your own risk!
const DWORD ProcessInformationMemoryPriority = 0x27;
const DWORD ProcessInformationIoPriority = 0x21;
const DWORD DefaultMemoryPriority = 5;
const DWORD LowMemoryPriority = 3;
const DWORD DefaultIoPriority = 2;
const DWORD LowIoPriority = 1;
static void win_throw_on_error(DWORD r,std::string msg) {
    if (r!=0) throw std::runtime_error(msg + std::string(", error-code:")+std::to_string(GetLastError()));
}
static void win_set_priority(int p_class) {
    DWORD cpu, io, mem;
    if (p_class==0) {
        cpu = NORMAL_PRIORITY_CLASS;
        mem = DefaultMemoryPriority;
        io = DefaultIoPriority;
    } else if (p_class==-1) {
        cpu = BELOW_NORMAL_PRIORITY_CLASS;
        mem = LowMemoryPriority;
        io = LowIoPriority;
    } else {
        throw std::runtime_error("Only 0=normal and -1=low are supported for priority_class");
    }
    // locate the functions for querying/setting memory and IO priority
    // HMODULE ntdll = LoadLibrary("ntdll.dll");
    // NtSetInformationProcess = (NtSetInformationProcessFn)GetProcAddress(ntdll, "NtSetInformationProcess");

    auto me =GetCurrentProcess();// we don't need to release this accoring to https://msdn.microsoft.com/en-us/library/windows/desktop/ms683179(v=vs.85).aspx
    if (!me) throw std::runtime_error("Failed to get current process handle");
    win_throw_on_error(!SetPriorityClass(me, cpu), "Failed setting cpu-priority");
    win_throw_on_error(NtSetInformationProcess(me, ProcessInformationMemoryPriority,&mem, sizeof(mem)), "Failed setting mem-priority");
    win_throw_on_error(NtSetInformationProcess(me, ProcessInformationIoPriority, &io, sizeof(io)), "Failed setting mem-priority");
}

std::string win_short_path(const std::string& long_path) {
    long length = GetShortPathName(long_path.c_str(), NULL, 0);
    std::string r(length+1, '\0');
    length = GetShortPathName(long_path.c_str(),(char*) r.data(), length);
    r.resize(length);
    return r;
}
#else
//-- on linux we are not even close to having these problems
std::string win_short_path(const std::string& long_path) {
    return long_path;
}
void win_set_priority(int) {}
#endif

BOOST_PYTHON_MODULE(_api) {
    namespace py = boost::python;
    py::scope().attr("__doc__") = "Shyft python api providing basic types";
    py::def("version", version);
    py::docstring_options doc_options(true, true, false);// all except c++ signatures
    expose::api();
    // We register the function with the atexit module which
    // will be called _before_ the Python C-API is unloaded.
    // needed for proper clean-up on windows platform
    // otherwise python hangs on dlib::shared_ptr_thread_safe destruction
    py::def("win_short_path", win_short_path, py::arg("path"),
            doc_intro("WinApi function GetShortPath exposed to python")
            doc_intro("https://msdn.microsoft.com/en-us/library/windows/desktop/aa364989(v=vs.85).aspx")
            doc_intro("Note that it only works for file-paths that exists, returns null string for not-existing files")
            doc_parameters()
            doc_parameter("path","str","a long path form")
            doc_returns("short_path","str","windows 8.3 path string if on windows for *existing* files, otherwise same as input path")
            );
    py::def("win_set_priority", win_set_priority, py::arg("p_class"),
            doc_intro("Win32 Api function to set normal (=0) or low(=-1) priority")
            doc_intro("This is *very* specific to windows, and especially task-scheduler/bg.tasks")
            doc_intro("get by default a complete garbled priority leaving the process you run")
            doc_intro("close to useless when it comes to cpu,io and memory performance")
            doc_intro("The UI does not help fixing this, and even the xml edit of files does not solve memory/io issues")
            doc_parameters()
            doc_parameter("p_class","int","priority class, 0=normal, -1=below normal(slightly above useless")
    );
    py::def("_finalize", &finalize_api);
    py::object atexit = py::object(py::handle<>(PyImport_ImportModule("atexit")));
    py::object finalize = py::scope().attr("_finalize");
    atexit.attr("register")(finalize);
}

