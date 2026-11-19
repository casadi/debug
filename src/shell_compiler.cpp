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
#include <cctype>
#include <sstream>

#ifdef _WIN32
#include <windows.h>
#include <map>
#ifdef CASADI_WITH_THREAD
#include <mutex>
#endif // CASADI_WITH_THREAD
#else // _WIN32
#include <unistd.h>
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
    int win_run(const std::string& cmd, const WinEnv* env, std::string* out=nullptr,
        const char* dir=nullptr) {
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
      DWORD flags = GetConsoleCP() ? 0 : CREATE_NO_WINDOW;
      if (env) flags |= CREATE_UNICODE_ENVIRONMENT;
      std::vector<char> buf(cmd.begin(), cmd.end());
      buf.push_back('\0');
      PROCESS_INFORMATION pi;
      BOOL ok = CreateProcessA(nullptr, buf.data(), nullptr, nullptr, TRUE, flags,
        env ? const_cast<wchar_t*>(env->block.data()) : nullptr, dir, &si, &pi);
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

    // Expand variables like the shell did, unless the shell is used
    std::string expand(const std::string& s, bool shell=false) {
#ifdef _WIN32
      if (shell) return s;
      // %VAR% as cmd.exe
      std::string e = s;
      DWORD n = ExpandEnvironmentStringsA(s.c_str(), nullptr, 0);
      if (n) {
        e.resize(n+1);
        DWORD m = ExpandEnvironmentStringsA(s.c_str(), &e[0], n+1);
        if (m && m<=n+1) {
          e.resize(strlen(e.c_str()));
        } else {
          e = s;
        }
      }
      return e;
#else // _WIN32
      // Still interpreted by the shell
      return s;
#endif // _WIN32
    }

    // Quote a path
    std::string quote(const std::string& s, bool shell) {
      // Not to be taken for an option
      if (!s.empty() && (s[0]=='-' || s[0]=='@')) return quote("./" + s, shell);
#ifdef _WIN32
      // A multibyte code page could take a backslash or quote as a trail byte
      CPINFO cp;
      if (GetACP()!=CP_UTF8 && GetCPInfo(CP_ACP, &cp) && cp.MaxCharSize>1) {
        for (char c : s) {
          casadi_assert(static_cast<unsigned char>(c)<0x80, "Unsupported character in '" + s + "'");
        }
      }
      // MSVC runtime rules
      if (shell) {
        // cmd.exe expands % and ! also within quotes
        casadi_assert(s.find_first_of("%!\"")==std::string::npos,
          "Unsupported character in '" + s + "' for a command run through cmd.exe");
        if (!s.empty() && s.find_first_not_of("abcdefghijklmnopqrstuvwxyz"
            "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-+/.:@\\~#$")==std::string::npos) return s;
      } else if (!s.empty() && s.find_first_of(" \t\n\"")==std::string::npos) {
        return s;
      }
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
      // POSIX shell
      if (!s.empty() && s.find_first_not_of("abcdefghijklmnopqrstuvwxyz"
          "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-+=/.,:@%")==std::string::npos) return s;
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

#ifdef _WIN32
    // Does an option need cmd.exe?
    bool cmd_syntax(const std::string& s) {
      // Operators, also after expansion
      for (const std::string& e : {s, expand(s)}) {
        bool in_quotes = false;
        for (char c : e) {
          if (c=='"') {
            in_quotes = !in_quotes;
          } else if (!in_quotes && (c=='&' || c=='|' || c=='<' || c=='>' || c=='^')) {
            return true;
          }
        }
      }
      // Variables only cmd.exe expands, such as %CD%
      std::string e = expand(s);
      size_t i = e.find('%');
      return i!=std::string::npos && e.find('%', i+2)!=std::string::npos;
    }

    // Does cmd.exe end up outside quotes after an option?
    bool cmd_balanced(const std::string& s) {
      bool in_quotes = false;
      for (size_t i=0; i<s.size(); ++i) {
        if (s[i]=='"') {
          in_quotes = !in_quotes;
        } else if (!in_quotes && s[i]=='^' && ++i==s.size()) {
          return false;
        }
      }
      return !in_quotes;
    }

    // Current directory
    std::string win_cwd() {
      std::vector<char> buf(GetCurrentDirectoryA(0, nullptr) + 1);
      DWORD n = GetCurrentDirectoryA(static_cast<DWORD>(buf.size()), buf.data());
      return n && n<buf.size() ? std::string(buf.data(), n) : "";
    }

    // Program of a command line, found as argv[0] by the MSVC runtime
    std::string program(const std::string& cmd) {
      size_t b = cmd.find_first_not_of(" \t");
      if (b==std::string::npos) return "";
      if (cmd[b]=='"') return cmd.substr(b+1, cmd.find('"', b+1)-b-1);
      return cmd.substr(b, cmd.find_first_of(" \t", b)-b);
    }

    // Extensions cmd.exe tries
    std::vector<std::string> win_pathext() {
      const char* pathext = getenv("PATHEXT");
      std::stringstream ss(pathext ? pathext : ".COM;.EXE;.BAT;.CMD");
      std::vector<std::string> r;
      std::string e;
      while (std::getline(ss, e, ';')) if (!e.empty()) r.push_back(e);
      return r;
    }

    // Find a program on the PATH of an environment
    bool win_find(const std::string& prog, const WinEnv* env, std::string& full,
        const std::string& ext=".exe") {
      wchar_t wprog[MAX_PATH], wext[MAX_PATH], wfull[MAX_PATH];
      char afull[MAX_PATH];
      if (!MultiByteToWideChar(CP_ACP, 0, prog.c_str(), -1, wprog, MAX_PATH)) return false;
      if (!MultiByteToWideChar(CP_ACP, 0, ext.c_str(), -1, wext, MAX_PATH)) return false;
      DWORD n = SearchPathW(env ? env->path.c_str() : nullptr, wprog, wext,
        MAX_PATH, wfull, nullptr);
      if (n==0 || n>=MAX_PATH) return false;
      if (!WideCharToMultiByte(CP_ACP, 0, wfull, -1, afull, MAX_PATH, nullptr, nullptr)) {
        return false;
      }
      full = afull;
      return true;
    }

    // Is a file an executable that CreateProcess starts without cmd.exe?
    bool win_exe(std::string f) {
      for (char& c : f) c = static_cast<char>(tolower(static_cast<unsigned char>(c)));
      if (f.size()>=8 && f.compare(f.size()-8, 8, "\\cmd.exe")==0) return false;
      return f.size()>=4 && (f.compare(f.size()-4, 4, ".exe")==0
        || f.compare(f.size()-4, 4, ".com")==0);
    }

    // Can a program be started without cmd.exe, which would find the same?
    bool win_direct(const std::string& prog, const WinEnv* env) {
      std::string f;
      if (!win_find(prog, env, f) || !win_exe(f)) return false;
      // cmd.exe could find a script with the same name
      if (prog.find('.', prog.find_last_of("\\/:")+1)==std::string::npos) {
        for (const std::string& e : win_pathext()) {
          if (!win_exe(e) && win_find(prog, env, f, e)) return false;
        }
      }
      return true;
    }
#endif // _WIN32

    // Include directory, with enclosing quotes removed as by the shell
    std::string include_dir(const std::string& s) {
#ifdef _WIN32
      std::string e = expand(s);
      if (e.size()>=2 && e.front()=='"' && e.back()=='"') return e.substr(1, e.size()-2);
      return e;
#else // _WIN32
      if (s.size()>=2 && s.front()=='\'' && s.back()=='\'') return s.substr(1, s.size()-2);
      std::string e = s;
      if (e.size()>=2 && e.front()=='"' && e.back()=='"') e = e.substr(1, e.size()-2);
      // $VAR and ${VAR}
      std::string r;
      for (size_t i=0; i<e.size(); ++i) {
        if (e[i]=='$') {
          bool brace = i+1<e.size() && e[i+1]=='{';
          size_t j = brace ? i+2 : i+1, k = j;
          while (k<e.size() && (isalnum(static_cast<unsigned char>(e[k])) || e[k]=='_')) k++;
          if (k>j && (!brace || (k<e.size() && e[k]=='}'))) {
            std::string name = e.substr(j, k-j);
            char buf[4096];
            const char* v = name=="PWD" ? getcwd(buf, sizeof(buf)) : getenv(name.c_str());
            if (v) r += v;
            i = brace ? k : k-1;
            continue;
          }
        }
        r += e[i];
      }
      return r;
#endif // _WIN32
    }

    // Run a command
#ifdef _WIN32
    int run(const std::string& cmd, bool shell, const WinEnv* env) {
      if (shell) {
        const char* comspec = getenv("ComSpec");
        return win_run(quote(comspec ? comspec : "cmd.exe", false) + " /s /c \"" + cmd + "\"",
          env);
      }
      return win_run(cmd, env);
#else // _WIN32
    int run(const std::string& cmd, bool shell) {
      return system(cmd.c_str());
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
      if (it!=cache.end()) return it->second.block.empty() ? nullptr : &it->second;

      const char* comspec = getenv("ComSpec");
      std::string cmd = quote(comspec ? comspec : "cmd.exe", false)
        + " /u /d /s /c \"set VSCMD_SKIP_SENDTELEMETRY=1&& " + quote(bat, true) + " " + args
        + " >nul && set\"";
      std::string out;
      std::string dir = bat.substr(0, bat.find_last_of("\\/"));
      if (win_run(cmd, nullptr, &out, dir.c_str())) out.clear();

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
      if (!env.block.empty()) env.block += L'\0';
      const WinEnv& r = cache[key] = env;
      return r.block.empty() ? nullptr : &r;
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
      std::string dir = vswhere.substr(0, vswhere.find_last_of('\\'));
      if (win_run(quote(vswhere, false) + " -latest -products * -requires " + component
          + " -property installationPath", nullptr, &out, dir.c_str())) return "";
      size_t end = out.find_first_of("\r\n");
      if (end!=std::string::npos) out.resize(end);
      if (out.empty()) return "";
      return out + "\\VC\\Auxiliary\\Build\\vcvarsall.bat";
    }

    // Use the full path of the program
    void win_resolve(std::string& cmd, const WinEnv* env) {
      std::string prog = program(cmd), full;
      cmd.erase(0, cmd.find_first_not_of(" \t"));
      size_t end = cmd[0]=='"' ? prog.size()+2 : prog.size();
      if (win_find(prog, env, full)) {
        cmd = quote(full, false) + (end<cmd.size() ? cmd.substr(end) : "");
      }
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

    // Run through cmd.exe as before unless the compiler and linker can be started directly
    bool shell = false;
#ifdef _WIN32
    std::vector<std::string> all = {compiler, compiler_setup, compiler_output_flag,
      compiler_include_flag, linker, linker_setup, linker_output_flag};
    all.insert(all.end(), compiler_flags.begin(), compiler_flags.end());
    all.insert(all.end(), linker_flags.begin(), linker_flags.end());
    for (const std::string& s : all) {
      if (cmd_syntax(s)) shell = true;
    }
    // Set up Visual Studio if cl.exe is not on the PATH
    const WinEnv* env = nullptr;
    std::string prog = program(expand(compiler)), full;
    bool found = false;
    for (const std::string& e : win_pathext()) {
      if (win_find(prog, nullptr, full, e)) found = true;
    }
    if ((_stricmp(prog.c_str(), "cl")==0 || _stricmp(prog.c_str(), "cl.exe")==0) && !found) {
      std::string bat = vcvarsall();
      if (!bat.empty()) {
#if defined(_M_ARM64) || defined(__aarch64__)
        env = win_env(bat, "arm64");
#elif defined(_WIN64)
        env = win_env(bat, "x64");
#else
        env = win_env(bat, "x86");
#endif
        // Otherwise as before
        if (env && !win_find(prog, env, full)) env = nullptr;
      }
    }
    for (const std::string& s : {compiler, linker}) {
      if (!win_direct(program(expand(s)), env)) shell = true;
    }
    if (shell) {
      for (const std::string& s : all) {
        // CasADi's paths must not end up between quotes
        casadi_assert(cmd_balanced(s) && cmd_balanced(expand(s)),
          "Unbalanced quotes or trailing ^ in '" + s + "'");
        // Nor the current directory, which cmd.exe would expand unquoted
        std::string u = s;
        for (char& c : u) c = static_cast<char>(toupper(static_cast<unsigned char>(c)));
        if (u.find("%CD%")!=std::string::npos || u.find("%__CD__%")!=std::string::npos
            || u.find("%=")!=std::string::npos) {
          casadi_assert(win_cwd().find_first_of("&|<>^()%!\"")==std::string::npos,
            "Unsupported character in the current directory for '" + s + "'");
        }
      }
    }
#endif // _WIN32

    // Compiler and linker programs
    std::string cc = expand(compiler, shell), ld = expand(linker, shell);
#ifdef _WIN32
    // CreateProcess searches the PATH of the parent
    if (env && !shell) {
      win_resolve(cc, env);
      win_resolve(ld, env);
    }
#endif // _WIN32

    // Construct the compiler command
    std::stringstream cccmd;
    cccmd << cc;
    for (auto i=compiler_flags.begin(); i!=compiler_flags.end(); ++i) {
      cccmd << " " << expand(*i, shell);
    }
    // Include directories (OS-agnostic: prefix with the compiler's include flag)
    for (const std::string& dir : include_dirs) {
      cccmd << " " << expand(compiler_include_flag, shell)
        << quote(include_dir(dir), shell);
    }
    cccmd << " " << expand(compiler_setup, shell);

    // C/C++ source file
    cccmd << " " << quote(name_, shell);

    // Temporary object file
    cccmd << " " + expand(compiler_output_flag, shell) << quote(obj_name_, shell);

    // Compile into an object
    if (verbose_) casadi_message("calling \"" + cccmd.str() + "\"");
#ifdef _WIN32
    int ret = run(cccmd.str(), shell, env);
#else // _WIN32
    int ret = run(cccmd.str(), shell);
#endif // _WIN32
    if (ret) {
      casadi_error("Compilation failed. Tried \"" + cccmd.str() + "\"");
    }

    // Link step
    std::stringstream ldcmd;
    ldcmd << ld;

    // Temporary file
    ldcmd << " " << quote(obj_name_, shell) << " " + expand(linker_output_flag, shell)
      + quote(bin_name_, shell);

    // Add flags
    for (auto i=linker_flags.begin(); i!=linker_flags.end(); ++i) {
      ldcmd << " " << expand(*i, shell);
    }
    ldcmd << " " << expand(linker_setup, shell);

    // Compile into a shared library
    if (verbose_) casadi_message("calling \"" + ldcmd.str() + "\"");
#ifdef _WIN32
    ret = run(ldcmd.str(), shell, env);
#else // _WIN32
    ret = run(ldcmd.str(), shell);
#endif // _WIN32
    if (ret) {
      casadi_error("Linking failed. Tried \"" + ldcmd.str() + "\"");
    }

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
