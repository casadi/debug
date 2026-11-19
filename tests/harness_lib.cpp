// The plugin source in a shared library, for loading with RTLD_DEEPBIND
#include "../src/shell_compiler.cpp"

namespace casadi {
const Options ImporterInternal::options_ = {};
const std::string ShellCompiler::meta_doc = "";
}  // namespace casadi

extern "C" int harness_run(const char* const* argv) {
  std::vector<std::string> args;
  for (; *argv; ++argv) args.push_back(*argv);
  return casadi::run(args);
}

extern "C" int harness_environ_is_live(void) {
  void* h = dlopen(nullptr, RTLD_LAZY);
  char*** p = reinterpret_cast<char***>(dlsym(h, "environ"));
  int live = p && *p==environ;
  dlclose(h);
  return live;
}
