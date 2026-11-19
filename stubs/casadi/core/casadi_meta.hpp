#pragma once
// Minimal stand-in for the harness
#define CASADI_VERSION 0
namespace casadi {
struct CasadiMeta {
  static const char* object_file_suffix() {
#ifdef _WIN32
    return ".obj";
#else
    return ".o";
#endif
  }
};
} // namespace casadi
