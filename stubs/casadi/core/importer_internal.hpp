#pragma once
// Minimal stand-in for the harness
#include "casadi_misc.hpp"
#ifdef _WIN32
// As casadi_os.hpp
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <dlfcn.h>
#endif

#ifdef _WIN32
#define SHARED_LIBRARY_SUFFIX ".dll"
#elif defined(__APPLE__)
#define SHARED_LIBRARY_SUFFIX ".dylib"
#else
#define SHARED_LIBRARY_SUFFIX ".so"
#endif

namespace casadi {

enum TypeID {OT_STRING, OT_STRINGVECTOR, OT_BOOL};

class GenericType {
 public:
  GenericType() : t_(T_NONE), b_(false) {}
  GenericType(const char* s) : t_(T_STR), s_(s), b_(false) {}
  GenericType(const std::string& s) : t_(T_STR), s_(s), b_(false) {}
  GenericType(bool b) : t_(T_BOOL), b_(b) {}
  GenericType(const std::vector<std::string>& v) : t_(T_VEC), b_(false), v_(v) {}
  std::string to_string() const { check(T_STR); return s_; }
  std::vector<std::string> to_string_vector() const { check(T_VEC); return v_; }
  operator bool() const { check(T_BOOL); return b_; }
  operator std::vector<std::string>() const { return to_string_vector(); }
 private:
  enum T {T_NONE, T_STR, T_BOOL, T_VEC};
  void check(T t) const { if (t!=t_) throw CasadiException("GenericType: wrong type"); }
  T t_;
  std::string s_;
  bool b_;
  std::vector<std::string> v_;
};

typedef std::map<std::string, GenericType> Dict;

struct Options {
  struct Entry {
    TypeID type;
    std::string description;
  };
  std::vector<const Options*> bases;
  std::map<std::string, Entry> entries;
};

typedef void (*signal_t)(void);
#ifdef _WIN32
typedef HINSTANCE handle_t;
#else
typedef void* handle_t;
#endif

inline handle_t open_shared_library(const std::string& lib,
    const std::vector<std::string>& /*search_paths*/, const std::string& f) {
#ifdef _WIN32
  handle_t h = LoadLibraryA(lib.c_str());
  if (!h) throw CasadiException(f + ": LoadLibrary failed for '" + lib + "', error "
    + std::to_string(GetLastError()));
#else
  handle_t h = dlopen(lib.c_str(), RTLD_LAZY | RTLD_LOCAL);
  if (!h) throw CasadiException(f + ": dlopen failed for '" + lib + "': " + dlerror());
#endif
  return h;
}

inline void close_shared_library(handle_t h) {
#ifdef _WIN32
  FreeLibrary(h);
#else
  dlclose(h);
#endif
}

struct FunctionInternal {
  // Directory from the options; must end with a separator
  static std::string get_jit_directory(const Dict& opts) {
    auto it = opts.find("directory");
    return it==opts.end() ? std::string() : it->second.to_string();
  }
};

class ImporterInternal {
 public:
  explicit ImporterInternal(const std::string& name) : name_(name), verbose_(false) {}
  virtual ~ImporterInternal() {}
  struct Plugin {
    ImporterInternal* (*creator)(const std::string&);
    const char* name;
    const char* doc;
    int version;
    const Options* options;
  };
  static void registerPlugin(int (*)(Plugin*)) {}
  static const Options options_;
  virtual const Options& get_options() const { return options_; }
  virtual void init(const Dict& opts) {
    auto it = opts.find("verbose");
    if (it!=opts.end()) verbose_ = it->second;
  }
  virtual const char* plugin_name() const = 0;
  virtual std::string class_name() const = 0;
  virtual signal_t get_function(const std::string& symname) = 0;
  virtual std::string library() const { return ""; }
  std::vector<std::string> get_search_paths() const { return {}; }
 protected:
  std::string name_;
  bool verbose_;
};

} // namespace casadi
