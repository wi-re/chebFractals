#include <tools/exceptionHandling.h>
#include <iostream>



#ifdef _WIN32
    const char* InfoFromSE::opDescription(const ULONG opcode)
    {
        switch (opcode) {
        case 0: return "read";
        case 1: return "write";
        case 8: return "user-mode data execution prevention (DEP) violation";
        default: return "unknown";
        }
    }

    const char* InfoFromSE::seDescription(const exception_code_t& code)
    {
        switch (code) {
        case EXCEPTION_ACCESS_VIOLATION:         return "EXCEPTION_ACCESS_VIOLATION";
        case EXCEPTION_ARRAY_BOUNDS_EXCEEDED:    return "EXCEPTION_ARRAY_BOUNDS_EXCEEDED";
        case EXCEPTION_BREAKPOINT:               return "EXCEPTION_BREAKPOINT";
        case EXCEPTION_DATATYPE_MISALIGNMENT:    return "EXCEPTION_DATATYPE_MISALIGNMENT";
        case EXCEPTION_FLT_DENORMAL_OPERAND:     return "EXCEPTION_FLT_DENORMAL_OPERAND";
        case EXCEPTION_FLT_DIVIDE_BY_ZERO:       return "EXCEPTION_FLT_DIVIDE_BY_ZERO";
        case EXCEPTION_FLT_INEXACT_RESULT:       return "EXCEPTION_FLT_INEXACT_RESULT";
        case EXCEPTION_FLT_INVALID_OPERATION:    return "EXCEPTION_FLT_INVALID_OPERATION";
        case EXCEPTION_FLT_OVERFLOW:             return "EXCEPTION_FLT_OVERFLOW";
        case EXCEPTION_FLT_STACK_CHECK:          return "EXCEPTION_FLT_STACK_CHECK";
        case EXCEPTION_FLT_UNDERFLOW:            return "EXCEPTION_FLT_UNDERFLOW";
        case EXCEPTION_ILLEGAL_INSTRUCTION:      return "EXCEPTION_ILLEGAL_INSTRUCTION";
        case EXCEPTION_IN_PAGE_ERROR:            return "EXCEPTION_IN_PAGE_ERROR";
        case EXCEPTION_INT_DIVIDE_BY_ZERO:       return "EXCEPTION_INT_DIVIDE_BY_ZERO";
        case EXCEPTION_INT_OVERFLOW:             return "EXCEPTION_INT_OVERFLOW";
        case EXCEPTION_INVALID_DISPOSITION:      return "EXCEPTION_INVALID_DISPOSITION";
        case EXCEPTION_NONCONTINUABLE_EXCEPTION: return "EXCEPTION_NONCONTINUABLE_EXCEPTION";
        case EXCEPTION_PRIV_INSTRUCTION:         return "EXCEPTION_PRIV_INSTRUCTION";
        case EXCEPTION_SINGLE_STEP:              return "EXCEPTION_SINGLE_STEP";
        case EXCEPTION_STACK_OVERFLOW:           return "EXCEPTION_STACK_OVERFLOW";
        default: return "UNKNOWN EXCEPTION";
        }
    }

    std::string InfoFromSE::information(struct _EXCEPTION_POINTERS* ep, bool has_exception_code, exception_code_t code)
    {
        HMODULE hm;
        ::GetModuleHandleEx(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS, static_cast<LPCTSTR>(ep->ExceptionRecord->ExceptionAddress), &hm);
        MODULEINFO mi;
        ::GetModuleInformation(::GetCurrentProcess(), hm, &mi, sizeof(mi));
        char fn[MAX_PATH];
        ::GetModuleFileNameExA(::GetCurrentProcess(), hm, fn, MAX_PATH);

        std::ostringstream oss;
        oss << "SE " << (has_exception_code ? seDescription(code) : "") << " at address 0x" << std::hex << ep->ExceptionRecord->ExceptionAddress << std::dec
            << " inside " << fn << " loaded at base address 0x" << std::hex << mi.lpBaseOfDll << "\n";

        if (has_exception_code && (
            code == EXCEPTION_ACCESS_VIOLATION ||
            code == EXCEPTION_IN_PAGE_ERROR)) {
            oss << "Invalid operation: " << opDescription(ep->ExceptionRecord->ExceptionInformation[0]) << " at address 0x" << std::hex << ep->ExceptionRecord->ExceptionInformation[1] << std::dec << "\n";
        }

        if (has_exception_code && code == EXCEPTION_IN_PAGE_ERROR) {
            oss << "Underlying NTSTATUS code that resulted in the exception " << ep->ExceptionRecord->ExceptionInformation[2] << "\n";
        }

        return oss.str();
    }
#include <dbghelp.h>
#include <intrin.h>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <windows.h>

#pragma comment(lib, "dbghelp.lib")

#define DBG_TRACE(MSG, ...) ::dbg::trace(MSG, __VA_ARGS__)

#define DBG_SOFT_ASSERT(COND)                                                                                          \
  if ((COND) == false) {                                                                                               \
    DBG_TRACE(__FUNCTION__ ": Assertion '" #COND "' failed!\n");                                                       \
  }

#define DBG_ASSERT(COND)                                                                                               \
  if ((COND) == false) {                                                                                               \
    DBG_TRACE(__FUNCTION__ ": Assertion '" #COND "' failed!\n");                                                       \
    ::dbg::handle_assert(__FUNCTION__, #COND);                                                                         \
  }

#define DBG_FAIL(MSG)                                                                                                  \
  DBG_TRACE(__FUNCTION__ MSG "\n");                                                                                    \
  ::dbg::fail(__FUNCTION__, MSG);

namespace dbg {
    inline void trace(const char* msg, ...) {
        char buff[1024];

        va_list args;
        va_start(args, msg);
        vsnprintf(buff, 1024, msg, args);

        OutputDebugStringA(buff);

        va_end(args);
    }

    inline std::string basename(const std::string& file) {
        unsigned int i = file.find_last_of("\\/");
        if (i == std::string::npos) {
            return file;
        }
        else {
            return file.substr(i + 1);
        }
    }


    std::vector<StackFrame> stack_trace() {
#if _WIN64
        DWORD machine = IMAGE_FILE_MACHINE_AMD64;
#else
        DWORD machine = IMAGE_FILE_MACHINE_I386;
#endif
        HANDLE process = GetCurrentProcess();
        HANDLE thread = GetCurrentThread();

        if (SymInitialize(process, NULL, TRUE) == FALSE) {
            DBG_TRACE(__FUNCTION__ ": Failed to call SymInitialize.");
            return std::vector<StackFrame>();
        }

        SymSetOptions(SYMOPT_LOAD_LINES);

        CONTEXT context = {};
        context.ContextFlags = CONTEXT_FULL;
        RtlCaptureContext(&context);

#if _WIN64
        STACKFRAME frame = {};
        frame.AddrPC.Offset = context.Rip;
        frame.AddrPC.Mode = AddrModeFlat;
        frame.AddrFrame.Offset = context.Rbp;
        frame.AddrFrame.Mode = AddrModeFlat;
        frame.AddrStack.Offset = context.Rsp;
        frame.AddrStack.Mode = AddrModeFlat;
#else
        STACKFRAME frame = {};
        frame.AddrPC.Offset = context.Eip;
        frame.AddrPC.Mode = AddrModeFlat;
        frame.AddrFrame.Offset = context.Ebp;
        frame.AddrFrame.Mode = AddrModeFlat;
        frame.AddrStack.Offset = context.Esp;
        frame.AddrStack.Mode = AddrModeFlat;
#endif

        bool first = true;

        std::vector<StackFrame> frames;
        while (StackWalk(machine, process, thread, &frame, &context, NULL, SymFunctionTableAccess, SymGetModuleBase, NULL)) {
            StackFrame f = {};
            f.address = frame.AddrPC.Offset;

#if _WIN64
            DWORD64 moduleBase = 0;
#else
            DWORD moduleBase = 0;
#endif

            moduleBase = SymGetModuleBase(process, frame.AddrPC.Offset);

            char moduelBuff[MAX_PATH];
            if (moduleBase && GetModuleFileNameA((HINSTANCE)moduleBase, moduelBuff, MAX_PATH)) {
                f.module = basename(moduelBuff);
            }
            else {
                f.module = "Unknown Module";
            }
#if _WIN64
            DWORD64 offset = 0;
#else
            DWORD offset = 0;
#endif
            char symbolBuffer[sizeof(IMAGEHLP_SYMBOL) + 255];
            PIMAGEHLP_SYMBOL symbol = (PIMAGEHLP_SYMBOL)symbolBuffer;
            symbol->SizeOfStruct = (sizeof IMAGEHLP_SYMBOL) + 255;
            symbol->MaxNameLength = 254;

            if (SymGetSymFromAddr(process, frame.AddrPC.Offset, &offset, symbol)) {
                f.name = symbol->Name;
            }
            else {
                DWORD error = GetLastError();
                DBG_TRACE(__FUNCTION__ ": Failed to resolve address 0x%X: %u\n", frame.AddrPC.Offset, error);
                f.name = "Unknown Function";
            }

            IMAGEHLP_LINE line;
            line.SizeOfStruct = sizeof(IMAGEHLP_LINE);

            DWORD offset_ln = 0;
            if (SymGetLineFromAddr(process, frame.AddrPC.Offset, &offset_ln, &line)) {
                f.file = line.FileName;
                f.line = line.LineNumber;
            }
            else {
                DWORD error = GetLastError();
                DBG_TRACE(__FUNCTION__ ": Failed to resolve line for 0x%X: %u\n", frame.AddrPC.Offset, error);
                f.line = 0;
            }

            if (!first) {
                frames.push_back(f);
            }
            first = false;
        }

        SymCleanup(process);

        return frames;
    }

    inline void handle_assert(const char* func, const char* cond) {
        std::stringstream buff;
        buff << func << ": Assertion '" << cond << "' failed! \n";
        buff << "\n";

        std::vector<StackFrame> stack = stack_trace();
        buff << "Callstack: \n";
        for (unsigned int i = 0; i < stack.size(); i++) {
            buff << "0x" << std::hex << stack[i].address << ": " << stack[i].name << "(" << std::dec << stack[i].line << ") in "
                << stack[i].module << "\n";
        }

        MessageBoxA(NULL, buff.str().c_str(), "Assert Failed", MB_OK | MB_ICONSTOP);
        abort();
    }

    inline void fail(const char* func, const char* msg) {
        std::stringstream buff;
        buff << func << ":  General Software Fault: '" << msg << "'! \n";
        buff << "\n";

        std::vector<StackFrame> stack = stack_trace();
        buff << "Callstack: \n";
        for (unsigned int i = 0; i < stack.size(); i++) {
            buff << "0x" << std::hex << stack[i].address << ": " << stack[i].name << "(" << stack[i].line << ") in "
                << stack[i].module << "\n";
        }

        MessageBoxA(NULL, buff.str().c_str(), "General Software Fault", MB_OK | MB_ICONSTOP);
        abort();
    }
}


void translator(InfoFromSE::exception_code_t code, struct _EXCEPTION_POINTERS* ep)
{
    auto trace = dbg::stack_trace();
    for (auto t : trace) {
        std::cerr << t.address << ": " << t.file << " - " << t.line << " => " << t.name << " @ " << t.module << std::endl;
    }

    throw std::exception(InfoFromSE::information(ep, true, code).c_str());
}

#endif

// The exception filter function:
// LONG WINAPI ExpFilter(EXCEPTION_POINTERS* pExp, DWORD dwExpCode) {
//     auto trace = dbg::stack_trace();
//     for (auto t : trace) {
//         std::cerr << t.address << ": " << t.file << " - " << t.line << " => " << t.name << " @ " << t.module << std::endl;
//     }
//     //StackWalker sw;
//     //sw.ShowCallstack(GetCurrentThread(), pExp->ContextRecord);
//     return EXCEPTION_EXECUTE_HANDLER;
// }