/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


#include "shell_compiler.hpp"
#include "casadi/core/casadi_misc.hpp"
#include "casadi/core/casadi_meta.hpp"
#include "casadi/core/casadi_logger.hpp"
#include <fstream>

// Set default object file suffix
#ifndef OBJECT_FILE_SUFFIX
#define OBJECT_FILE_SUFFIX CasadiMeta::object_file_suffix()
#endif // OBJECT_FILE_SUFFIX

#include <cstdlib>
#include <cstring>

#ifdef _WIN32
#include <windows.h>
#include <map>
#ifdef CASADI_WITH_THREAD
#include <mutex>
#endif // CASADI_WITH_THREAD
#else // _WIN32
#include <spawn.h>
#include <sys/wait.h>
#include <cerrno>
#include <wordexp.h>
#ifdef __APPLE__
#include <crt_externs.h>
#else // __APPLE__
#include <dlfcn.h>
extern char **environ;
#endif // __APPLE__
#endif // _WIN32

namespace casadi {

namespace {
#ifdef _WIN32
    // Environment for the compiler and linker
    struct WinEnv {
      std::wstring block;
      std::wstring path;
    };

    // Run a command line without a shell, optionally capturing stdout
    int win_run(const std::string& cmd, const WinEnv* env, std::string* out=nullptr) {
      SECURITY_ATTRIBUTES sa = {sizeof(SECURITY_ATTRIBUTES), nullptr, TRUE};
      HANDLE rd = nullptr, wr = nullptr;
      STARTUPINFOA si;
      ZeroMemory(&si, sizeof(si));
      si.cb = sizeof(si);
      if (out) {
        if (!CreatePipe(&rd, &wr, &sa, 0)) return -1;
        SetHandleInformation(rd, HANDLE_FLAG_INHERIT, 0);
        si.dwFlags = STARTF_USESTDHANDLES;
        si.hStdInput = GetStdHandle(STD_INPUT_HANDLE);
        si.hStdOutput = wr;
        si.hStdError = GetStdHandle(STD_ERROR_HANDLE);
      }
      // No console window when the parent has none
      DWORD flags = GetConsoleWindow() ? 0 : CREATE_NO_WINDOW;
      if (env) flags |= CREATE_UNICODE_ENVIRONMENT;
      std::vector<char> buf(cmd.begin(), cmd.end());
      buf.push_back('\0');
      PROCESS_INFORMATION pi;
      BOOL ok = CreateProcessA(nullptr, buf.data(), nullptr, nullptr, TRUE, flags,
        env ? const_cast<wchar_t*>(env->block.data()) : nullptr, nullptr, &si, &pi);
      if (wr) CloseHandle(wr);
      if (!ok) {
        if (rd) CloseHandle(rd);
        return -1;
      }
      if (out) {
        char chunk[4096];
        DWORD n;
        while (ReadFile(rd, chunk, sizeof(chunk), &n, nullptr) && n>0) out->append(chunk, n);
        CloseHandle(rd);
      }
      WaitForSingleObject(pi.hProcess, INFINITE);
      DWORD ret = 1;
      GetExitCodeProcess(pi.hProcess, &ret);
      CloseHandle(pi.hProcess);
      CloseHandle(pi.hThread);
      return static_cast<int>(ret);
    }
#endif // _WIN32

    // Expand variables and split into arguments, like the shell did
    std::vector<std::string> expand(const std::string& s) {
#ifdef _WIN32
      // %VAR% as cmd.exe
      std::string e = s;
      DWORD n = ExpandEnvironmentStringsA(s.c_str(), nullptr, 0);
      if (n) {
        e.resize(n);
        n = ExpandEnvironmentStringsA(s.c_str(), &e[0], n);
        e.resize(n ? strlen(e.c_str()) : 0);
        if (!n) e = s;
      }
      // Split as the MSVC runtime
      std::vector<std::string> r;
      std::string arg;
      bool in_arg = false, in_quotes = false;
      for (size_t i=0; i<e.size(); ++i) {
        char c = e[i];
        if (c=='\\') {
          size_t n_backslash = 0;
          while (i<e.size() && e[i]=='\\') {
            n_backslash++;
            i++;
          }
          if (i<e.size() && e[i]=='"') {
            arg.append(n_backslash/2, '\\');
            if (n_backslash%2) {
              arg += '"';
            } else {
              i--;
            }
          } else {
            arg.append(n_backslash, '\\');
            i--;
          }
          in_arg = true;
        } else if (c=='"') {
          if (in_quotes && i+1<e.size() && e[i+1]=='"') {
            arg += e[++i];
          } else {
            in_quotes = !in_quotes;
          }
          in_arg = true;
        } else if ((c==' ' || c=='\t') && !in_quotes) {
          if (in_arg) r.push_back(arg);
          arg.clear();
          in_arg = false;
        } else {
          arg += c;
          in_arg = true;
        }
      }
      if (in_arg) r.push_back(arg);
      return r;
#else // _WIN32
      wordexp_t we;
      int ret = wordexp(s.c_str(), &we, WRDE_NOCMD);
      casadi_assert(ret!=WRDE_BADCHAR, "Unsupported character in '" + s + "': "
        "newline, |, &, ;, <, >, (, ), { and } are not interpreted");
      casadi_assert(ret!=WRDE_CMDSUB, "Command substitution not supported in '" + s + "'");
      casadi_assert(ret==0, "Cannot expand '" + s + "'");
      std::vector<std::string> r(we.we_wordv, we.we_wordv+we.we_wordc);
      wordfree(&we);
      return r;
#endif // _WIN32
    }

    // Quote an argument
    std::string quote(const std::string& s) {
#ifdef _WIN32
      // MSVC runtime rules
      if (!s.empty() && s.find_first_of(" \t\n\"")==std::string::npos) return s;
      std::string r = "\"";
      size_t n_backslash = 0;
      for (char c : s) {
        if (c=='\\') {
          n_backslash++;
        } else {
          if (c=='"') r.append(n_backslash+1, '\\');
          n_backslash = 0;
        }
        r += c;
      }
      r.append(n_backslash, '\\');
      return r + "\"";
#else // _WIN32
      // For display
      if (!s.empty() && s.find_first_not_of("abcdefghijklmnopqrstuvwxyz"
          "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789@%_-+=:,./")==std::string::npos) return s;
      std::string r = "'";
      for (char c : s) {
        if (c=='\'') {
          r += "'\\''";
        } else {
          r += c;
        }
      }
      return r + "'";
#endif // _WIN32
    }

    // Append a flag followed by its value; a trailing space makes them separate arguments
    void add_flag(std::vector<std::string>& args, const std::string& flag,
        const std::string& value) {
      std::vector<std::string> f = expand(flag);
      bool glued = !f.empty() && flag.back()!=' ';
      args.insert(args.end(), f.begin(), f.end());
      if (glued) {
        args.back() += value;
      } else {
        args.push_back(value);
      }
    }

    // Join arguments into a command line
    std::string cmdline(const std::vector<std::string>& args) {
      std::string r;
      for (const std::string& a : args) r += (r.empty() ? "" : " ") + quote(a);
      return r;
    }

    // Run a command
#ifdef _WIN32
    int run(const std::vector<std::string>& args, const WinEnv* env) {
      return win_run(cmdline(args), env);
#else // _WIN32
    int run(const std::vector<std::string>& args) {
      std::vector<char*> argv;
      for (const std::string& a : args) argv.push_back(const_cast<char*>(a.c_str()));
      argv.push_back(nullptr);
#ifdef __APPLE__
      char** env = *_NSGetEnviron();
#else // __APPLE__
      // Under RTLD_DEEPBIND, environ binds to a stale copy: take the executable's
      char** env = environ;
      void* h = dlopen(nullptr, RTLD_LAZY);
      if (h) {
        char*** p = reinterpret_cast<char***>(dlsym(h, "environ"));
        if (p) env = *p;
        dlclose(h);
      }
#endif // __APPLE__
      char* no_env[] = {nullptr};
      if (!env) env = no_env;
      pid_t pid;
      if (posix_spawnp(&pid, argv[0], nullptr, nullptr, argv.data(), env)) return -1;
      int status;
      while (waitpid(pid, &status, 0)==-1) {
        if (errno!=EINTR) return -1;
      }
      if (WIFEXITED(status)) return WEXITSTATUS(status);
      return 1;
#endif // _WIN32
    }

#ifdef _WIN32
    // Environment after running a batch file, cached
    const WinEnv* win_env(const std::string& bat, const std::string& args) {
      static std::map<std::string, WinEnv> cache;
#ifdef CASADI_WITH_THREAD
      static std::mutex mtx;
      std::lock_guard<std::mutex> lock(mtx);
#endif // CASADI_WITH_THREAD
      std::string key = bat + " " + args;
      auto it = cache.find(key);
      if (it!=cache.end()) return &it->second;

      const char* comspec = getenv("ComSpec");
      std::string cmd = quote(comspec ? comspec : "cmd.exe") + " /u /d /s /c \""
        + quote(bat) + " " + args + " >nul && set\"";
      std::string out;
      casadi_assert(win_run(cmd, nullptr, &out)==0,
        "Environment setup failed. Tried '" + cmd + "'");

      // Output is UTF-16, one NAME=value per line
      std::wstring wout(reinterpret_cast<const wchar_t*>(out.data()), out.size()/sizeof(wchar_t));
      WinEnv env;
      size_t pos = 0;
      while (pos<wout.size()) {
        size_t end = wout.find(L'\n', pos);
        if (end==std::wstring::npos) end = wout.size();
        std::wstring line = wout.substr(pos, end-pos);
        pos = end+1;
        if (!line.empty() && line.back()==L'\r') line.pop_back();
        size_t eq = line.find(L'=');
        if (eq==std::wstring::npos || eq==0) continue;
        if (_wcsicmp(line.substr(0, eq).c_str(), L"PATH")==0) env.path = line.substr(eq+1);
        env.block += line;
        env.block += L'\0';
      }
      casadi_assert(!env.block.empty(), "Environment setup failed. Tried '" + cmd + "'");
      env.block += L'\0';
      return &(cache[key] = env);
    }

    // Locate vcvarsall.bat of the latest Visual Studio with a C++ toolset
    std::string win_vcvarsall() {
      const char* pf = getenv("ProgramFiles(x86)");
      if (!pf) pf = getenv("ProgramFiles");
      if (!pf) return "";
      std::string vswhere = pf;
      vswhere += "\\Microsoft Visual Studio\\Installer\\vswhere.exe";
#if defined(_M_ARM64) || defined(__aarch64__)
      std::string component = "Microsoft.VisualStudio.Component.VC.Tools.ARM64";
#else
      std::string component = "Microsoft.VisualStudio.Component.VC.Tools.x86.x64";
#endif
      std::string out;
      if (win_run(quote(vswhere) + " -latest -products * -requires " + component
          + " -property installationPath", nullptr, &out)) return "";
      size_t end = out.find_first_of("\r\n");
      if (end!=std::string::npos) out.resize(end);
      if (out.empty()) return "";
      return out + "\\VC\\Auxiliary\\Build\\vcvarsall.bat";
    }

    // Find a program on the PATH of an environment
    bool win_find(const std::string& prog, const WinEnv* env, std::string& full) {
      wchar_t wprog[MAX_PATH], wfull[MAX_PATH];
      char afull[MAX_PATH];
      if (!MultiByteToWideChar(CP_ACP, 0, prog.c_str(), -1, wprog, MAX_PATH)) return false;
      DWORD n = SearchPathW(env ? env->path.c_str() : nullptr, wprog, L".exe",
        MAX_PATH, wfull, nullptr);
      if (n==0 || n>=MAX_PATH) return false;
      if (!WideCharToMultiByte(CP_ACP, 0, wfull, -1, afull, MAX_PATH, nullptr, nullptr)) {
        return false;
      }
      full = afull;
      return true;
    }

    // Use the full path of the program
    void win_resolve(std::vector<std::string>& args, const WinEnv* env) {
      std::string full;
      if (!args.empty() && win_find(args[0], env, full)) args[0] = full;
    }
#endif // _WIN32
}  // namespace

  extern "C"
  int CASADI_IMPORTER_SHELL_EXPORT
  casadi_register_importer_shell(ImporterInternal::Plugin* plugin) {
    plugin->creator = ShellCompiler::creator;
    plugin->name = "shell";
    plugin->doc = ShellCompiler::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &ShellCompiler::options_;
    return 0;
  }

  extern "C"
  void CASADI_IMPORTER_SHELL_EXPORT casadi_load_importer_shell() {
    ImporterInternal::registerPlugin(casadi_register_importer_shell);
  }

#ifdef _WIN32
  std::string ShellCompiler::vcvarsall_;
  bool ShellCompiler::vcvarsall_searched_ = false;
#ifdef CASADI_WITH_THREAD
  std::mutex ShellCompiler::mutex_vcvarsall_;
#endif // CASADI_WITH_THREAD
#endif // _WIN32

#ifdef _WIN32
  std::string ShellCompiler::vcvarsall() {
#ifdef CASADI_WITH_THREAD
    std::lock_guard<std::mutex> lock(mutex_vcvarsall_);
#endif // CASADI_WITH_THREAD
    if (!vcvarsall_searched_) {
      vcvarsall_ = win_vcvarsall();
      vcvarsall_searched_ = true;
    }
    return vcvarsall_;
  }
#endif // _WIN32

  ShellCompiler::ShellCompiler(const std::string& name) :
    ImporterInternal(name) {
      handle_ = nullptr;
  }

  ShellCompiler::~ShellCompiler() {
    if (handle_) close_shared_library(handle_);

    if (cleanup_) {
      if (remove(bin_name_.c_str())) casadi_warning("Failed to remove " + bin_name_);
      if (remove(obj_name_.c_str())) casadi_warning("Failed to remove " + obj_name_);
      for (const std::string& s : extra_suffixes_) {
        std::string name = base_name_+s;
        remove(name.c_str());
      }
    }
  }

  const Options ShellCompiler::options_
  = {{&ImporterInternal::options_},
     {{"compiler",
       {OT_STRING,
        "Compiler command"}},
      {"linker",
       {OT_STRING,
        "Linker command"}},
      {"directory",
       {OT_STRING,
        "Directory to put temporary objects in. Must end with a file separator."}},
      {"compiler_setup",
       {OT_STRING,
        "Compiler setup command. Intended to be fixed."
        " The 'flag' option is the preferred way to set"
        " custom flags."}},
      {"linker_setup",
       {OT_STRING,
        "Linker setup command. Intended to be fixed."
        " The 'flag' option is the preferred way to set"
        " custom flags."}},
      {"compiler_flags",
       {OT_STRINGVECTOR,
        "Alias for 'compiler_flags'"}},
      {"flags",
        {OT_STRINGVECTOR,
        "Compile flags for the JIT compiler. Default: None"}},
      {"linker_flags",
       {OT_STRINGVECTOR,
        "Linker flags for the JIT compiler. Default: None"}},
      {"cleanup",
       {OT_BOOL,
        "Cleanup temporary files when unloading. Default: true"}},
      {"compiler_output_flag",
       {OT_STRING,
       "Compiler flag to denote object output. Default: '-o '"}},
      {"compiler_include_flag",
       {OT_STRING,
       "Compiler flag to add an include directory. Default: '-I' ('/I' on MSVC)"}},
      {"include_dirs",
       {OT_STRINGVECTOR,
       "List of include directories, each passed to the compiler prefixed with "
       "'compiler_include_flag'. OS-agnostic alternative to passing '-I...' via 'flags'. "
       "Default: None"}},
      {"linker_output_flag",
       {OT_STRING,
       "Linker flag to denote shared library output. Default: '-o '"}},
      {"extra_suffixes",
       {OT_STRINGVECTOR,
       "List of suffixes for extra files that the compiler may generate. Default: None"}},
      {"name",
       {OT_STRING,
        "The file name used to write out compiled objects/libraries. "
        "The actual file names used depend on 'temp_suffix' and include extensions. "
        "Default: 'tmp_casadi_compiler_shell'"}},
      {"temp_suffix",
       {OT_BOOL,
        "Use a temporary (seemingly random) filename suffix for file names. "
        "This is desired for thread-safety. "
        "This behaviour may defeat caching compiler wrappers. "
        "Default: true"}},
     }
  };

  void ShellCompiler::init(const Dict& opts) {
    // Base class
    ImporterInternal::init(opts);

    // Default options

    cleanup_ = true;
    bool temp_suffix = true;
    std::string bare_name = "tmp_casadi_compiler_shell";
    std::string directory = FunctionInternal::get_jit_directory(opts);

    std::vector<std::string> compiler_flags;
    std::vector<std::string> linker_flags;
    std::vector<std::string> include_dirs;
    std::string suffix = OBJECT_FILE_SUFFIX;

#ifdef _WIN32
    std::string compiler = "cl.exe";
    std::string linker = "link.exe";
    std::string compiler_setup = "/c";
    std::string linker_setup = "/DLL";
    std::string compiler_output_flag = "/Fo";
    std::string linker_output_flag = "/out:";
    std::string compiler_include_flag = "/I";
    extra_suffixes_ = {".exp", ".lib"};
#elif defined(__APPLE__)
    std::string compiler = "clang";
    std::string linker = "clang";
    std::string compiler_setup = "-fPIC -c";
    std::string linker_setup = "-shared";
    std::string compiler_output_flag = "-o ";
    std::string linker_output_flag = "-o ";
    std::string compiler_include_flag = "-I";
#else
    std::string compiler = "gcc";
    std::string linker = "gcc";
    std::string compiler_setup = "-fPIC -c";
    std::string linker_setup = "-shared";
    std::string compiler_output_flag = "-o ";
    std::string linker_output_flag = "-o ";
    std::string compiler_include_flag = "-I";
#endif

    // Read options
    for (auto&& op : opts) {
      if (op.first=="compiler") {
        compiler = op.second.to_string();
      } else if (op.first=="linker") {
        linker = op.second.to_string();
      } else if (op.first=="compiler_setup") {
        compiler_setup = op.second.to_string();
      } else if (op.first=="cleanup") {
        cleanup_ = op.second;
      } else if (op.first=="linker_setup") {
        linker_setup = op.second.to_string();
      } else if (op.first=="compiler_flags" || op.first=="flags") {
        compiler_flags = op.second;
      } else if (op.first=="linker_flags") {
        linker_flags = op.second;
      } else if (op.first=="compiler_output_flag") {
        compiler_output_flag = op.second.to_string();
      } else if (op.first=="linker_output_flag") {
        linker_output_flag = op.second.to_string();
      } else if (op.first=="compiler_include_flag") {
        compiler_include_flag = op.second.to_string();
      } else if (op.first=="include_dirs") {
        include_dirs = op.second.to_string_vector();
      } else if (op.first=="extra_suffixes") {
        extra_suffixes_ = op.second.to_string_vector();
      } else if (op.first=="name") {
        bare_name = op.second.to_string();
      } else if (op.first=="temp_suffix") {
        temp_suffix = op.second;
      }
    }

    // Name of temporary file
    if (temp_suffix) {
      obj_name_ = temporary_file(bare_name, suffix, directory);
    } else {
      obj_name_ = directory + bare_name + suffix;
    }
    base_name_ = std::string(obj_name_.begin(), obj_name_.begin()+obj_name_.size()-suffix.size());
    bin_name_ = base_name_+SHARED_LIBRARY_SUFFIX;

    // Compiler command
    std::vector<std::string> cc = expand(compiler);
    for (const std::string& f : compiler_flags) {
      std::vector<std::string> e = expand(f);
      cc.insert(cc.end(), e.begin(), e.end());
    }
    for (const std::string& dir : include_dirs) {
      std::vector<std::string> e = expand(compiler_include_flag + dir);
      cc.insert(cc.end(), e.begin(), e.end());
    }
    std::vector<std::string> setup = expand(compiler_setup);
    cc.insert(cc.end(), setup.begin(), setup.end());
    cc.push_back(name_);
    add_flag(cc, compiler_output_flag, obj_name_);

    // Linker command
    std::vector<std::string> ld = expand(linker);
    ld.push_back(obj_name_);
    add_flag(ld, linker_output_flag, bin_name_);
    for (const std::string& f : linker_flags) {
      std::vector<std::string> e = expand(f);
      ld.insert(ld.end(), e.begin(), e.end());
    }
    setup = expand(linker_setup);
    ld.insert(ld.end(), setup.begin(), setup.end());
    casadi_assert(!cc.empty() && !ld.empty(), "Compiler and linker must not be empty");

#ifdef _WIN32
    // Environment of compiler and linker
    const WinEnv* env = nullptr;
    std::string full;
    if ((cc[0]=="cl" || cc[0]=="cl.exe") && !win_find(cc[0], nullptr, full)) {
      std::string bat = vcvarsall();
#if defined(_M_ARM64) || defined(__aarch64__)
      if (!bat.empty()) env = win_env(bat, "arm64");
#elif defined(_WIN64)
      if (!bat.empty()) env = win_env(bat, "x64");
#else
      if (!bat.empty()) env = win_env(bat, "x86");
#endif
    }
    if (env) {
      win_resolve(cc, env);
      win_resolve(ld, env);
    }
#endif // _WIN32

    // Compile into an object
    if (verbose_) casadi_message("calling \"" + cmdline(cc) + "\"");
#ifdef _WIN32
    int ret = run(cc, env);
#else // _WIN32
    int ret = run(cc);
#endif // _WIN32
    casadi_assert(ret!=-1, "Could not start compiler '" + cc[0] + "'");
    casadi_assert(ret==0, "Compilation failed. Tried \"" + cmdline(cc) + "\"");

    // Link into a shared library
    if (verbose_) casadi_message("calling \"" + cmdline(ld) + "\"");
#ifdef _WIN32
    ret = run(ld, env);
#else // _WIN32
    ret = run(ld);
#endif // _WIN32
    casadi_assert(ret!=-1, "Could not start linker '" + ld[0] + "'");
    casadi_assert(ret==0, "Linking failed. Tried \"" + cmdline(ld) + "\"");

    std::vector<std::string> search_paths = get_search_paths();
    handle_ = open_shared_library(bin_name_, search_paths, "ShellCompiler::init");

  }

  std::string ShellCompiler::library() const {
    return bin_name_;
  }

  signal_t ShellCompiler::get_function(const std::string& symname) {
#ifdef _WIN32
    return (signal_t)GetProcAddress(handle_, TEXT(symname.c_str()));
#else // _WIN32
    signal_t fcnPtr = reinterpret_cast<signal_t>(dlsym(handle_, symname.c_str()));
    if (dlerror()) {
      fcnPtr=nullptr;
      dlerror(); // Reset error flags
    }
    return fcnPtr;
#endif // _WIN32
  }

} // namespace casadi
