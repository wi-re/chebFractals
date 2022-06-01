#pragma once
#include <boost/exception/diagnostic_information.hpp> 
#ifdef _WIN32
#include <signal.h>
#include <windows.h>
#include <eh.h>
#include <Psapi.h>
#include <sstream>
#include <string>
#include <vector>
class InfoFromSE
{
public:
	typedef unsigned int exception_code_t;
    static const char* opDescription(const ULONG opcode);
	static const char* seDescription(const exception_code_t& code);
	static std::string information(struct _EXCEPTION_POINTERS* ep, bool has_exception_code = false, exception_code_t code = 0);
};

namespace dbg {
    struct StackFrame {
        DWORD64 address;
        std::string name;
        std::string module;
        unsigned int line;
        std::string file;
    };

    inline void trace(const char* msg, ...);
    inline std::string basename(const std::string& file);
    inline std::vector<StackFrame> stack_trace();
    inline void handle_assert(const char* func, const char* cond);
    inline void fail(const char* func, const char* msg);
}
void translator(InfoFromSE::exception_code_t code, struct _EXCEPTION_POINTERS* ep);
LONG WINAPI ExpFilter(EXCEPTION_POINTERS* pExp, DWORD dwExpCode);

#define CATCH_DEFAULT                                                                                                    \
  catch (...) {                                                                                                        \
    std::cerr << "Caught exception while running simulation" << std::endl;                                             \
    std::cerr << boost::current_exception_diagnostic_information() << std::endl;                                       \
    auto trace = dbg::stack_trace();                                                                                   \
    for (auto t : trace) {                                                                                             \
      std::cerr << t.address << ": " << t.file << " - " << t.line << " => " << t.name << " @ " << t.module             \
                << std::endl;                                                                                          \
    }                                                                                                                  \
  }

#define CATCH_FN(x)                                                                                                    \
  catch (...) {                                                                                                        \
    std::cerr << "Caught exception while running simulation" << std::endl;                                             \
    std::cerr << boost::current_exception_diagnostic_information() << std::endl;                                       \
    auto trace = dbg::stack_trace();                                                                                   \
    for (auto t : trace) {                                                                                             \
      std::cerr << t.address << ": " << t.file << " - " << t.line << " => " << t.name << " @ " << t.module             \
                << std::endl;                                                                                          \
    }                                                                                                                  \
    x;                                                                                                                 \
  }
#endif
