#pragma once
// Minimal stand-in for the harness
#include <string>
#include <vector>
#include <map>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <atomic>
#ifdef _WIN32
#include <process.h>
#else
#include <unistd.h>
#endif
namespace casadi {
struct CasadiException : std::runtime_error {
  explicit CasadiException(const std::string& s) : std::runtime_error(s) {}
};
inline std::string temporary_file(const std::string& prefix, const std::string& suffix,
    const std::string& directory) {
  static std::atomic<int> counter(0);
#ifdef _WIN32
  int pid = _getpid();
#else
  int pid = getpid();
#endif
  return directory + prefix + "_" + std::to_string(pid) + "_" + std::to_string(counter++) + suffix;
}
} // namespace casadi
#define casadi_assert(cond, msg) \
  do { if (!(cond)) throw casadi::CasadiException(std::string("assertion failed: ") + (msg)); } while (0)
#define casadi_error(msg) throw casadi::CasadiException(msg)
#define casadi_warning(msg) (std::cerr << "warning: " << (msg) << std::endl)
#define casadi_message(msg) (std::cout << (msg) << std::endl)
