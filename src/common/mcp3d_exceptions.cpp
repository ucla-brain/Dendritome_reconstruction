//
// Created by muyezhu on 3/26/18.
//

#include "mcp3d_exceptions.hpp"

using namespace std;

void mcp3d::PrintNested(const exception &e, ostream &out, int level)
{
    out << e.what() << endl;
    try {
        rethrow_if_nested(e);
    } catch(const exception& e) {
        PrintNested(e, out, level + 1);
    } catch(...) { cout << "catch all caught rethrowing level " << level << endl;}
}

void mcp3d::PrintNestedException(const exception_ptr &eptr, const char *file,
                                 int line, const char *function, ostream &out)
{
    try
    {
        if (eptr)
            rethrow_exception(eptr);
    }
    catch (const exception& e)
    {
        out << "nested exception at " << file << " line " << line
            << ", function " << function << endl << endl;
        mcp3d::PrintNested(e, out);
    }
}

void mcp3d::ReThrow(const exception_ptr& eptr, const char *file,
                    int line, const char *function)
{
    try
    {
        if (eptr)
            rethrow_exception(eptr);
    }
    catch (const mcp3d::MCPRuntimeError& e)
    {
        std::throw_with_nested(MCPRuntimeError("", file, line, function));
    }
    catch (const mcp3d::MCPOSError& e)
    {
        std::throw_with_nested(MCPOSError("", file, line, function));
    }
    catch (const mcp3d::MCPAssertionError& e)
    {
        std::throw_with_nested(MCPAssertionError("", file, line, function));
    }
    catch (const mcp3d::MCPInvalidArgument& e)
    {
        std::throw_with_nested(MCPInvalidArgument("", file, line, function));
    }
    catch (const mcp3d::MCPDomainError& e)
    {
        std::throw_with_nested(MCPDomainError("", file, line, function));
    }
    catch (const mcp3d::MCPOutOfRangeError& e)
    {
        std::throw_with_nested(MCPOutOfRangeError("", file, line, function));
    }
    catch (const mcp3d::MCPBadAlloc& e)
    {
        std::throw_with_nested(MCPBadAlloc("", file, line, function));
    }
    catch (const system_error& e)
    {
        std::string message(std::string("at ") + file + std::to_string(line) + " " + function);
        std::throw_with_nested(runtime_error(string("system_error").append(message)));
    }
    catch (const runtime_error& e)
    {
        std::string message(std::string("at ") + file + std::to_string(line) + " " + function);
        std::throw_with_nested(runtime_error(string("runtime_error ").append(message)));
    }
    catch (const invalid_argument& e)
    {
        std::string message(std::string("at ") + file + std::to_string(line) + " " + function);
        std::throw_with_nested(invalid_argument(string("invalid_argument ").append(message)));
    }
    catch (const bad_alloc& e)
    {
        std::string message(std::string("at ") + file + std::to_string(line) + " " + function);
        std::throw_with_nested(runtime_error(string("bad_alloc ").append(message)));
    }
    catch (...)
    {
        std::string message(std::string("at ") + file + std::to_string(line) + " " + function);
        std::throw_with_nested(runtime_error(string("unknown error ").append(message)));
    }
}

void mcp3d::MultiThreadExceptions::CaptureException()
{
    unique_lock<mutex> guard(lock_);
    e_ptr_ = current_exception();
}