// Tests for the command line handling of the shell compiler plugin.
// The plugin source is included verbatim; stubs/ provides the few casadi bits it needs.
#ifdef HARNESS_CANDIDATE
#include "../candidate/shell_compiler.cpp"
#else
#include "../src/shell_compiler.cpp"
#endif

namespace casadi {
const Options ImporterInternal::options_ = {};
const std::string ShellCompiler::meta_doc = "";
}  // namespace casadi

#include <cstdio>
#include <fstream>
#include <random>
#include <thread>
#include <set>
#include <sys/stat.h>
#ifdef _WIN32
#include <direct.h>
#include <shellapi.h>
#else
#include <unistd.h>
#endif

using casadi::Dict;
using casadi::ShellCompiler;

#ifdef _WIN32
// Expose the protected static helper
struct Peek : ShellCompiler {
  using ShellCompiler::vcvarsall;
};
#endif

// ---------------------------------------------------------------- reporting
static int n_pass = 0, n_fail = 0, n_note = 0;
static std::string g_section;

static std::string esc(const std::string& s, size_t maxlen = 160) {
  std::string r = "\"";
  for (size_t i=0; i<s.size(); ++i) {
    if (r.size() > maxlen) {
      r += "...(" + std::to_string(s.size()) + " bytes)";
      break;
    }
    unsigned char c = s[i];
    if (c=='\\') r += "\\\\";
    else if (c=='"') r += "\\\"";
    else if (c=='\n') r += "\\n";
    else if (c=='\t') r += "\\t";
    else if (c=='\r') r += "\\r";
    else if (c<32 || c>=127) {
      char b[8];
      snprintf(b, sizeof(b), "\\x%02X", c);
      r += b;
    } else {
      r += static_cast<char>(c);
    }
  }
  return r + "\"";
}

static std::string vesc(const std::vector<std::string>& v) {
  std::string r = "[";
  for (size_t i=0; i<v.size(); ++i) r += (i ? ", " : "") + esc(v[i], 80);
  return r + "]";
}

static void pass() { n_pass++; }
static void fail(const std::string& what) {
  n_fail++;
  std::cout << "FAIL [" << g_section << "] " << what << std::endl;
}
static void note(const std::string& what) {
  n_note++;
  std::cout << "NOTE [" << g_section << "] " << what << std::endl;
}
static void check(bool ok, const std::string& what) {
  if (ok) pass(); else fail(what);
}
static void section(const std::string& s) {
  g_section = s;
  std::cout << "== " << s << std::endl;
}

// ---------------------------------------------------------------- utilities
static std::string g_work;          // scratch root
static std::string g_outfile;       // argv_dump output file
static std::string g_jit = "none";  // cl, gcc, cc, none
struct DumpExe { std::string label, path; };
static std::vector<DumpExe> g_dumps;

#ifdef _WIN32
static const char SEP = '\\';
#else
static const char SEP = '/';
#endif

static std::string join(const std::string& a, const std::string& b) {
  if (!a.empty() && (a.back()=='/' || a.back()=='\\')) return a + b;
  return a + SEP + b;
}

static bool exists(const std::string& p) {
  struct stat st;
  return stat(p.c_str(), &st)==0;
}

static bool make_dir(const std::string& p) {
#ifdef _WIN32
  return _mkdir(p.c_str())==0 || exists(p);
#else
  return mkdir(p.c_str(), 0755)==0 || exists(p);
#endif
}

static void make_dirs(const std::string& p) {
  for (size_t i=1; i<=p.size(); ++i) {
    if (i==p.size() || p[i]=='/' || p[i]=='\\') {
      std::string q = p.substr(0, i);
      if (!q.empty() && q.back()!=':') make_dir(q);
    }
  }
}

static bool change_dir(const std::string& p) {
#ifdef _WIN32
  return _chdir(p.c_str())==0;
#else
  return chdir(p.c_str())==0;
#endif
}

static std::string current_dir() {
  char b[4096];
#ifdef _WIN32
  return _getcwd(b, sizeof(b)) ? b : "";
#else
  return getcwd(b, sizeof(b)) ? b : "";
#endif
}

static std::string read_file(const std::string& p) {
  std::ifstream f(p.c_str(), std::ios::binary);
  std::stringstream ss;
  ss << f.rdbuf();
  return ss.str();
}

static void write_file(const std::string& p, const std::string& s) {
  std::ofstream f(p.c_str(), std::ios::binary);
  f << s;
}

static void set_env(const std::string& k, const std::string& v) {
#ifdef _WIN32
  SetEnvironmentVariableA(k.c_str(), v.c_str());
  _putenv_s(k.c_str(), v.c_str());
#else
  setenv(k.c_str(), v.c_str(), 1);
#endif
}

static void unset_env(const std::string& k) {
#ifdef _WIN32
  SetEnvironmentVariableA(k.c_str(), nullptr);
  _putenv_s(k.c_str(), "");
#else
  unsetenv(k.c_str());
#endif
}

static std::string unhex(const std::string& h) {
  if (h=="-") return "";
  std::string r;
  for (size_t i=0; i+1<h.size(); i+=2) r += static_cast<char>(std::stoi(h.substr(i, 2), nullptr, 16));
  return r;
}

struct Dump {
  int ret = -999;
  bool ok = false;
  std::vector<std::string> crt, w;
  std::string raw, env;
  bool env_set = false;
  int console = -1;
};

static Dump parse_dump(const std::string& text, int ret) {
  Dump d;
  d.ret = ret;
  std::istringstream ss(text);
  std::string line;
  while (std::getline(ss, line)) {
    size_t sp = line.find(' ');
    std::string tag = line.substr(0, sp), val = sp==std::string::npos ? "" : line.substr(sp+1);
    if (tag=="crt") d.ok = true;
    else if (tag=="a") d.crt.push_back(unhex(val));
    else if (tag=="b") d.w.push_back(unhex(val));
    else if (tag=="raw") d.raw = unhex(val);
    else if (tag=="console") d.console = std::stoi(val);
    else if (tag=="env") {
      if (val!="(unset)") {
        d.env_set = true;
        d.env = unhex(val);
      }
    }
  }
  return d;
}

// Run a command line through the plugin's run(); argv_dump writes to g_outfile
static Dump dump_line(const std::string& line, bool shell) {
  remove(g_outfile.c_str());
  set_env("ADUMP_OUT", g_outfile);
#ifdef _WIN32
  int ret = casadi::run(line, shell, nullptr);
#else
  int ret = casadi::run(line, shell);
#endif
  return parse_dump(read_file(g_outfile), ret);
}

static bool non_ascii(const std::string& s) {
  for (char c : s) if (static_cast<unsigned char>(c)>=128) return true;
  return false;
}

// Does quote() refuse the string?
static bool refuses(const std::string& s, bool shell, std::string* msg=nullptr) {
  try {
    casadi::quote(s, shell);
    return false;
  } catch (std::exception& e) {
    if (msg) *msg = e.what();
    return true;
  }
}

// ---------------------------------------------------------------- corpus
static std::vector<std::string> corpus() {
  std::vector<std::string> c = {
    "", " ", "  ", "\t", "a", "a b", "a\tb", "a  b", " leading", "trailing ",
    "\"", "\"\"", "\"\"\"", "a\"b", "\"a b\"", "a\"\"b", "a\" \"b",
    "\\", "\\\\", "\\\\\\", "a\\", "a\\\\", "a b\\", "a b\\\\", "a b\\\\\\",
    "\\\"", "\\\\\"", "\\\\\\\"", "a\\\"b", "a\\\\\"b", "a\\\\\\\"b", "\\\\a", "a\\b\\c",
    "C:\\Program Files\\", "C:\\Program Files (x86)\\x\\", "\\\\server\\share\\a b\\",
    "&|<>^()%!'", "a&b", "a|b", "a^b", "a^&b", "a&&b", "a||b", "a>b", "a<b", "(a)", "a)b",
    "%PATH%", "%%", "100%", "%", "!PATH!", "'single'", "'a b'", "`tick`", "$HOME", "$(x)",
    ";", "a;b", "*", "?", "~", "{a,b}", "#hash", "a=b", "@file", "-", "--", "a,b", "a;b c",
    "\xC3\xA9t\xC3\xA9", "caf\xE9", "\xE4\xB8\xAD", "\x81\x8D\x8F\x90\x9D",
    "line\nbreak", "cr\rx", "vt\vx", "ff\fx",
    "R&D (x86)", "a&mkdir sfxdir&b", "c^&x", "O'Brien", "R&md,sfxdir&b", "x)&md,sfxdir&(y",
  };
  std::string long1(8000, 'a');
  long1 += " & ^ \\ end\\";
  c.push_back(long1);
  return c;
}

// ---------------------------------------------------------------- quote round trip
static void test_quote_roundtrip() {
  for (int shell=0; shell<2; ++shell) {
    section(std::string("quote: argv round trip, shell=") + (shell ? "true" : "false"));
    for (const DumpExe& d : g_dumps) {
      for (const std::string& s : corpus()) {
        std::string msg;
        bool refused = refuses(s, shell!=0, &msg);
#ifdef _WIN32
        if (shell) {
          bool should = s.find_first_of("%!\"")!=std::string::npos;
          check(refused==should, "quote(" + esc(s) + ", true) refused=" + std::to_string(refused));
          if (refused) continue;
          // cmd.exe cannot pass these through a /c line
          if (s.find_first_of("\r\n")!=std::string::npos) continue;
        }
#endif
        if (refused) {
          fail("quote(" + esc(s) + ") refused: " + msg);
          continue;
        }
        std::string q = casadi::quote(s, shell!=0);
        // cmd.exe lines are limited to 8191 characters
        bool once = s.size()>3000;
        std::string line = casadi::quote(d.path, false) + " " + q + (once ? "" : " sep " + q);
        Dump r = dump_line(line, shell!=0);
        std::vector<std::string> want = {s, "sep", s};
        if (once) want = {s};
        std::string what = d.label + " input=" + esc(s) + " quoted=" + esc(q);
        if (!r.ok) {
          fail(what + " child did not run, ret=" + std::to_string(r.ret));
          continue;
        }
        if (r.crt==want) {
          pass();
        } else if (non_ascii(s)) {
          note("non-ASCII (code page) " + what + " crt=" + vesc(r.crt));
        } else {
          fail(what + " expected=" + vesc(want) + " crt=" + vesc(r.crt));
        }
      }
    }
  }
}

#ifdef _WIN32
// ---------------------------------------------------------------- Windows unit tests
static void test_units() {
  section("cmd_syntax");
  struct C { std::string s; bool want; };
  std::vector<C> cs = {
    {"", false}, {"cl.exe", false}, {">nul", true}, {"2>nul", true}, {"a|b", true}, {"a<b", true},
    {"a^b", true}, {"a&b", true}, {"\"a&b\"", false}, {"\"a\"&b", true}, {"\"a&b", false},
    {"-DX=\"a b\"", false}, {"%X%", false}, {"!X!", false}, {"(a)", false},
    {"call x.bat && cl", true}, {"/DTAG=a^&b", true}, {"\"\"&", true},
  };
  for (const C& c : cs) {
    bool got = casadi::cmd_syntax(c.s);
    check(got==c.want, "cmd_syntax(" + esc(c.s) + ")=" + std::to_string(got));
  }

  section("program");
  struct P { std::string cmd, prog, rest; };
  std::vector<P> ps = {
    {"cl.exe /c", "cl.exe", " /c"},
    {"cl", "cl", ""},
    {"\"C:\\Program Files\\VS\\cl.exe\" /c", "C:\\Program Files\\VS\\cl.exe", " /c"},
    {"\"C:\\Program Files\\VS\\cl.exe\"", "C:\\Program Files\\VS\\cl.exe", ""},
    {"\"cl.exe\"", "cl.exe", ""},
    {"cl\t/c", "cl", "\t/c"},
    {"\"unterminated", "unterminated", ""},
    {"\"\"", "", ""},
  };
  for (const P& p : ps) {
    size_t end = 0;
    std::string prog = casadi::program(p.cmd, end);
    std::string rest = end<=p.cmd.size() ? p.cmd.substr(end) : "?";
    check(prog==p.prog && rest==p.rest, "program(" + esc(p.cmd) + ") gave " + esc(prog)
      + " rest " + esc(rest) + ", expected " + esc(p.prog) + " rest " + esc(p.rest));
  }

  section("expand: %VAR%");
  set_env("HARNESS_V", "hello world");
  unset_env("HARNESS_UNDEF");
  struct E { std::string in, want; };
  std::vector<E> es = {
    {"%HARNESS_V%", "hello world"}, {"%harness_v%", "hello world"},
    {"%HARNESS_UNDEF%", "%HARNESS_UNDEF%"}, {"%%", "%%"}, {"%", "%"}, {"100%", "100%"},
    {"a%HARNESS_V%b", "ahello worldb"}, {"", ""}, {"!HARNESS_V!", "!HARNESS_V!"},
  };
  for (const E& e : es) {
    std::string got = casadi::expand(e.in);
    check(got==e.want, "expand(" + esc(e.in) + ")=" + esc(got) + ", expected " + esc(e.want));
  }
}

// ---------------------------------------------------------------- win_run
static void test_win_run() {
  section("win_run");
  if (g_dumps.empty()) return;
  const std::string& exe = g_dumps[0].path;
  std::string qexe = casadi::quote(exe, false);
  for (int code : {0, 1, 3, 42, 255, 1000}) {
    set_env("ADUMP_EXIT", std::to_string(code));
    Dump r = dump_line(qexe, false);
    check(r.ret==code, "exit code " + std::to_string(code) + " gave " + std::to_string(r.ret));
    r = dump_line(qexe, true);
    check(r.ret==code, "shell: exit code " + std::to_string(code) + " gave " + std::to_string(r.ret));
  }
  unset_env("ADUMP_EXIT");
  {
    int ret = casadi::run("no_such_program_harness.exe x", false, nullptr);
    check(ret==-1, "nonexistent program gave " + std::to_string(ret));
    ret = casadi::run("no_such_program_harness.exe x", true, nullptr);
    check(ret!=0 && ret!=-1, "shell: nonexistent program gave " + std::to_string(ret));
  }
  // Environment block
  {
    remove(g_outfile.c_str());
    casadi::WinEnv env;
    const char* sr = getenv("SystemRoot");
    std::vector<std::string> vals = {"ADUMP_ENV=HARNESS_BLOCK", "ADUMP_OUT=" + g_outfile,
      "HARNESS_BLOCK=from block", std::string("SystemRoot=") + (sr ? sr : "C:\\Windows")};
    for (const std::string& v : vals) {
      env.block += std::wstring(v.begin(), v.end());
      env.block += L'\0';
    }
    env.block += L'\0';
    int ret = casadi::win_run(qexe, &env);
    Dump r = parse_dump(read_file(g_outfile), ret);
    check(ret==0 && r.env=="from block", "env block: ret=" + std::to_string(ret)
      + " env=" + esc(r.env));
  }
  // Stdout capture
  {
    unset_env("ADUMP_OUT");
    set_env("ADUMP_SPAM", "3000000");
    set_env("ADUMP_EXIT", "7");
    std::string out;
    int ret = casadi::win_run(qexe + " cap", nullptr, &out);
    size_t nx = 0;
    for (char c : out) if (c=='x') nx++;
    check(ret==7 && nx==3000000, "capture: ret=" + std::to_string(ret) + " x count="
      + std::to_string(nx));
    unset_env("ADUMP_SPAM");
    unset_env("ADUMP_EXIT");
  }
  // Program in a directory with spaces and punctuation
  {
    std::string d = join(g_work, "dir with space & (x) ^ %HARNESS_V% !x!");
    make_dir(d);
    std::string copy = join(d, "argv dump.exe");
    CopyFileA(exe.c_str(), copy.c_str(), FALSE);
    Dump r = dump_line(casadi::quote(copy, false) + " \"a b\"", false);
    check(r.ok && r.crt==std::vector<std::string>{"a b"}, "program path " + esc(copy)
      + " ret=" + std::to_string(r.ret) + " crt=" + vesc(r.crt));
  }
}

static const char* jit_source_probe();
static bool jit_try_probe(const std::string& src, const std::string& outdir, std::string& err);

// Child mode: started without a console
static int console_probe(const std::string& out) {
  std::ostringstream ss;
  ss << "self_console " << (GetConsoleWindow() ? 1 : 0) << "\n";
  Dump r = dump_line(casadi::quote(g_dumps[0].path, false) + " probe", false);
  ss << "child_ret " << r.ret << "\nchild_console " << r.console << "\n";
  r = dump_line(casadi::quote(g_dumps[0].path, false) + " probe", true);
  ss << "shell_child_ret " << r.ret << "\nshell_child_console " << r.console << "\n";
  if (g_jit!="none") {
    std::string d = join(g_work, "probe jit & (x)");
    make_dir(d);
    std::string src = join(d, "jit src.c");
    write_file(src, jit_source_probe());
    std::string err;
    bool ok = jit_try_probe(src, d + SEP, err);
    ss << "probe_jit " << (ok ? "ok" : "failed: " + err) << "\n";
  }
  write_file(out, ss.str());
  return 0;
}

static void test_no_console() {
  section("win_run: parent without console");
  char self[MAX_PATH];
  GetModuleFileNameA(nullptr, self, MAX_PATH);
  std::string out = join(g_work, "probe.txt");
  remove(out.c_str());
  std::string cmd = casadi::quote(self, false) + " --work " + casadi::quote(g_work, false)
    + " --dump " + casadi::quote(g_dumps[0].label + "=" + g_dumps[0].path, false) + " --jit "
    + g_jit + " --console-probe " + casadi::quote(out, false);
  std::vector<char> buf(cmd.begin(), cmd.end());
  buf.push_back('\0');
  STARTUPINFOA si;
  ZeroMemory(&si, sizeof(si));
  si.cb = sizeof(si);
  PROCESS_INFORMATION pi;
  BOOL ok = CreateProcessA(nullptr, buf.data(), nullptr, nullptr, FALSE, DETACHED_PROCESS,
    nullptr, nullptr, &si, &pi);
  check(ok, "start detached probe");
  if (!ok) return;
  WaitForSingleObject(pi.hProcess, 120000);
  CloseHandle(pi.hProcess);
  CloseHandle(pi.hThread);
  std::string res = read_file(out);
  std::cout << res;
  check(res.find("self_console 0")!=std::string::npos, "probe has no console");
  check(res.find("\nchild_ret 0")!=std::string::npos, "child ran from a console-less parent");
  check(res.find("\nshell_child_ret 0")!=std::string::npos, "shell child ran");
  check(res.find("child_console 2")==std::string::npos, "child got a visible console window");
  if (g_jit!="none") check(res.find("probe_jit ok")!=std::string::npos, "JIT from a console-less parent");
}

static void test_vs_setup() {
  section("Visual Studio setup");
  std::string full;
  bool on_path = casadi::win_find("cl.exe", nullptr, full);
  std::cout << "   cl.exe on PATH: " << (on_path ? full : "no") << std::endl;
  std::string bat = Peek::vcvarsall();
  std::cout << "   vcvarsall: " << bat << std::endl;
  check(!bat.empty() && exists(bat), "vcvarsall.bat located");
  if (bat.empty()) return;
#if defined(_M_ARM64) || defined(__aarch64__)
  std::vector<std::string> archs = {"arm64", "arm64_x64"};
#elif defined(_WIN64)
  std::vector<std::string> archs = {"x64", "x64_x86", "x64_arm64"};
#else
  std::vector<std::string> archs = {"x86", "x86_x64"};
#endif
  for (const std::string& a : archs) {
    const casadi::WinEnv* env = nullptr;
    try {
      env = casadi::win_env(bat, a);
    } catch (std::exception& e) {
      fail("win_env(" + a + ") threw: " + e.what());
      continue;
    }
    std::string cl, link;
    bool f1 = casadi::win_find("cl.exe", env, cl), f2 = casadi::win_find("link.exe", env, link);
    std::cout << "   " << a << ": cl=" << cl << " link=" << link << std::endl;
    check(f1 && cl.find("\\MSVC\\")!=std::string::npos, a + ": cl.exe resolved: " + cl);
    check(f2 && link.find("\\MSVC\\")!=std::string::npos, a + ": link.exe resolved: " + link);
    check(casadi::win_env(bat, a)==env, a + ": cached");
  }
  // vcvarsall copied to a directory with cmd.exe syntax in its name
  std::string d = join(g_work, "vs copy & (x) ^ y");
  make_dir(d);
  std::string bat2 = join(d, "setup env.bat");
  write_file(bat2, "@call \"" + bat + "\" %*\r\n");
  try {
    const casadi::WinEnv* env = casadi::win_env(bat2, archs[0]);
    std::string cl;
    check(casadi::win_find("cl.exe", env, cl), "batch file in " + esc(d) + ": cl.exe resolved");
  } catch (std::exception& e) {
    fail("batch file in " + esc(d) + " threw: " + e.what());
  }
  std::string d3 = join(g_work, "vs copy 100%");
  make_dir(d3);
  try {
    casadi::win_env(join(d3, "x.bat"), archs[0]);
    fail("batch file in " + esc(d3) + " was not refused");
  } catch (std::exception& e) {
    check(std::string(e.what()).find("Unsupported character")!=std::string::npos,
      "batch file in " + esc(d3) + ": " + e.what());
  }
}
#endif

// ---------------------------------------------------------------- JIT
static const char* jit_source =
  "#ifdef _WIN32\n#define EXPORT __declspec(dllexport)\n#else\n#define EXPORT\n#endif\n"
  "#define XSTR(x) STR(x)\n#define STR(x) #x\n"
  "#ifdef USE_INC\n#include \"harness_inc.h\"\n#else\n#define INC_VALUE 0\n#endif\n"
  "EXPORT int square(int x) { return x*x + INC_VALUE; }\n"
  "#ifdef TAG\nEXPORT const char* tag(void) { return XSTR(TAG); }\n#endif\n"
  "#ifdef ARR\nstatic int arr[] = ARR;\nEXPORT int arr_len(void) { return sizeof(arr)/sizeof(arr[0]); }\n#endif\n";

static Dict jit_base() {
  Dict d;
  if (g_jit=="gcc") {
    d["compiler"] = "gcc";
    d["linker"] = "gcc";
    d["compiler_setup"] = "-c";
    d["linker_setup"] = "-shared";
    d["compiler_output_flag"] = "-o ";
    d["linker_output_flag"] = "-o ";
    d["compiler_include_flag"] = "-I";
    d["extra_suffixes"] = std::vector<std::string>{};
  }
  return d;
}

static std::string def_flag() { return g_jit=="cl" ? "/D" : "-D"; }
static std::string inc_flag() { return g_jit=="cl" ? "/I" : "-I"; }

struct Want {
  int square = 49;
  const char* tag = nullptr;
  int arr_len = -1;
};

static bool jit_try(const std::string& src, const std::string& outdir, Dict opts, const Want& want,
    std::string& err) {
  opts["directory"] = outdir;
  opts["verbose"] = true;
  try {
    ShellCompiler sc(src);
    sc.init(opts);
    typedef int (*sq_t)(int);
    sq_t f = reinterpret_cast<sq_t>(sc.get_function("square"));
    if (!f) {
      err = "symbol 'square' not found";
      return false;
    }
    if (f(7)!=want.square) {
      err = "square(7) gave " + std::to_string(f(7));
      return false;
    }
    if (want.tag) {
      typedef const char* (*tag_t)(void);
      tag_t t = reinterpret_cast<tag_t>(sc.get_function("tag"));
      if (!t) {
        err = "symbol 'tag' not found";
        return false;
      }
      if (std::string(t())!=want.tag) {
        err = "tag() gave " + esc(t()) + ", expected " + esc(want.tag);
        return false;
      }
    }
    if (want.arr_len>=0) {
      sq_t a = reinterpret_cast<sq_t>(sc.get_function("arr_len"));
      int n = a ? reinterpret_cast<int (*)(void)>(a)() : -1;
      if (n!=want.arr_len) {
        err = "arr_len() gave " + std::to_string(n);
        return false;
      }
    }
    return true;
  } catch (std::exception& e) {
    err = e.what();
    return false;
  }
}

static const char* jit_source_probe() { return jit_source; }
static bool jit_try_probe(const std::string& src, const std::string& outdir, std::string& err) {
  return jit_try(src, outdir, jit_base(), Want(), err);
}

// Markers that only appear if part of a name was run as a command
static std::vector<std::string> sentinels(const std::string& d) {
  std::vector<std::string> r;
  for (const char* n : {"sfxdir", "sfxdir2", "sfxdir3"}) {
    for (const std::string& base : {d, g_work, current_dir(), join(g_work, "jitdirs"),
        join(g_work, "shelldirs")}) {
      if (exists(join(base, n))) r.push_back(join(base, n));
    }
  }
  return r;
}

static std::vector<std::string> awkward_names() {
  return {
    "plain", "with space", "O'Brien", "R&D (x86)", "a&mkdir sfxdir&b", "c^d", "c^&x",
    "c^&mkdir sfxdir2", "R&md,sfxdir3&b", "x)&md,sfxdir3&(y", "100%", "pct%HARNESS_V%dir",
    "(paren)", "semi;colon", "bang!x!", "comma,eq=", "tilde~1", "\xC3\xA9t\xC3\xA9", "  two  spaces",
#ifndef _WIN32
    "x|y", "a;mkdir sfxdir;b", "dollar$HOME", "back`mkdir sfxdir`", "$(mkdir sfxdir2)",
    "quote\"dq", "star*", "brace{a,b}", "lt<gt>", "a\\b", "new\nline", "a&&mkdir sfxdir",
    "x'$(mkdir sfxdir3)'y", "-dash",
#endif
  };
}

// Compile in each awkward directory: cwd, source, output and include dir
static void jit_dirs(const std::string& label, const std::string& sub, const Dict& extra,
    bool shell) {
  section(label);
  set_env("HARNESS_V", "hello world");
  std::string root = join(g_work, sub);
  make_dir(root);
  std::string cwd = current_dir();
  for (const std::string& n : awkward_names()) {
    std::string d = join(root, n);
    if (!make_dir(d)) {
      note("cannot create directory " + esc(d));
      continue;
    }
    change_dir(d);
    std::string src = join(d, "jit src.c");
    write_file(src, jit_source);
    write_file(join(d, "harness_inc.h"), "#define INC_VALUE 1000\n");
    Dict o = jit_base();
    for (auto&& kv : extra) o[kv.first] = kv.second;
    std::vector<std::string> fl = {def_flag() + "USE_INC"};
    if (o.count("flags")) {
      std::vector<std::string> e = o["flags"];
      fl.insert(fl.end(), e.begin(), e.end());
    }
    o["flags"] = fl;
    o["include_dirs"] = std::vector<std::string>{d};
    std::string err;
    Want w;
    w.square = 1049;
    // Source compiled from cwd=d, so the quoted include also works without include_dirs
    write_file(join(d, "harness_inc.h"), "#define INC_VALUE 1000\n");
    bool ok = jit_try(src, d + SEP, o, w, err);
    std::vector<std::string> s = sentinels(d);
    check(s.empty(), "dir=" + esc(n) + ": side effect " + vesc(s));
#ifdef _WIN32
    bool expect_refusal = shell && n.find_first_of("%!")!=std::string::npos;
#else
    bool expect_refusal = false;
#endif
    if (expect_refusal) {
      check(!ok && err.find("Unsupported character")!=std::string::npos,
        "dir=" + esc(n) + ": expected 'Unsupported character', got " + (ok ? "success" : err));
    } else {
      check(ok, "dir=" + esc(n) + ": " + err);
    }
    change_dir(cwd);
  }
}

static void test_jit_options() {
  section("JIT options (" + g_jit + ")");
  std::string d = join(g_work, "jitopts dir");
  make_dir(d);
  std::string src = join(d, "jit src.c");
  write_file(src, jit_source);
  std::string err;
  set_env("HARNESS_V", "hello world");
  Want w;
  auto run = [&](const std::string& what, Dict o, const Want& want) {
    bool ok = jit_try(src, d + SEP, o, want, err);
    check(ok, what + ": " + err);
    std::vector<std::string> s = sentinels(d);
    check(s.empty(), what + ": side effect " + vesc(s));
  };

  // Include directories with spaces and punctuation
  for (const std::string& n : {"inc space", "inc dir & (x)", "inc 100%", "inc O'Brien"}) {
    std::string inc = join(d, n);
    make_dir(inc);
    write_file(join(inc, "harness_inc.h"), "#define INC_VALUE 1000\n");
    Dict o = jit_base();
    o["flags"] = std::vector<std::string>{def_flag() + "USE_INC"};
    o["include_dirs"] = std::vector<std::string>{inc};
    Want wi;
    wi.square = 1049;
    run("include_dirs " + esc(n), o, wi);
  }

  Dict o = jit_base();
  o["flags"] = std::vector<std::string>{def_flag() + "TAG=\"a b\""};
  Want wt;
  wt.tag = "a b";
  run("flags " + esc(def_flag() + "TAG=\"a b\""), o, wt);

#ifdef _WIN32
  o = jit_base();
  o["flags"] = std::vector<std::string>{"\"" + def_flag() + "TAG=%HARNESS_V%\""};
  wt.tag = "hello world";
  run("flags with %VAR% in quotes", o, wt);

  // cmd.exe syntax in options
  o = jit_base();
  o["flags"] = std::vector<std::string>{">nul"};
  run("flags [\">nul\"]", o, w);

  o = jit_base();
  o["flags"] = std::vector<std::string>{def_flag() + "TAG=a^&b"};
  wt.tag = "a&b";
  run("flags with ^& " + esc(def_flag() + "TAG=a^&b"), o, wt);

  o = jit_base();
  o["linker_flags"] = std::vector<std::string>{">nul", "2>&1"};
  run("linker_flags [\">nul\", \"2>&1\"]", o, w);

  if (g_jit=="cl") {
    o = jit_base();
    o["compiler"] = "cl.exe 2>nul";
    run("compiler \"cl.exe 2>nul\"", o, w);

    o = jit_base();
    o["compiler"] = "cl /nologo";
    run("compiler \"cl /nologo\"", o, w);

    set_env("HARNESS_CC", "cl.exe");
    o = jit_base();
    o["compiler"] = "%HARNESS_CC%";
    run("compiler \"%HARNESS_CC%\"", o, w);

    o = jit_base();
    o["compiler"] = " cl.exe";
    run("compiler \" cl.exe\" (leading space)", o, w);

    std::string full;
    std::string bat = Peek::vcvarsall();
    if (!bat.empty()) {
#if defined(_M_ARM64) || defined(__aarch64__)
      std::string arch = "arm64";
#elif defined(_WIN64)
      std::string arch = "x64";
#else
      std::string arch = "x86";
#endif
      const casadi::WinEnv* env = casadi::win_env(bat, arch);
      std::string clfull, linkfull;
      casadi::win_find("cl.exe", env, clfull);
      casadi::win_find("link.exe", env, linkfull);
      // Quoted full path with spaces
      o = jit_base();
      o["compiler"] = "\"" + clfull + "\"";
      run("compiler quoted full path (no Visual Studio setup)", o, w);
      // A setup batch file in the compiler command
      std::string sd = join(d, "setup & (x)");
      make_dir(sd);
      std::string sb = join(sd, "vcsetup.bat");
      write_file(sb, "@call \"" + bat + "\" " + arch + " >nul\r\n");
      o = jit_base();
      o["compiler"] = "call \"" + sb + "\" && cl";
      o["linker"] = "call \"" + sb + "\" && link";
      run("compiler \"call <bat> && cl\"", o, w);
    }
  }
#else
  o = jit_base();
  o["flags"] = std::vector<std::string>{"\"" + def_flag() + "TAG=$HARNESS_V\""};
  wt.tag = "hello world";
  run("flags with $VAR in quotes", o, wt);

  o = jit_base();
  o["flags"] = std::vector<std::string>{"2>/dev/null", ">/dev/null"};
  run("flags with redirection", o, w);

  o = jit_base();
  o["flags"] = std::vector<std::string>{"-DARR={1,2,3}"};
  Want wa;
  wa.arr_len = 3;
  bool ok = jit_try(src, d + SEP, o, wa, err);
  // Reference: what /bin/sh does with the word
  int r = system("test \"$(printf '%s ' -DARR={1,2,3})\" = '-DARR={1,2,3} '");
  note("flags [\"-DARR={1,2,3}\"]: " + std::string(ok ? "ok" : err.substr(0, 120))
    + "; /bin/sh brace expansion: " + (r==0 ? "no" : "yes"));

  const char* ccenv = getenv("CC");
  set_env("HARNESS_CC", g_jit=="gcc" ? "gcc" : (ccenv ? ccenv : "cc"));
  o = jit_base();
  o["compiler"] = "$HARNESS_CC";
  o["linker"] = "$HARNESS_CC";
  run("compiler $HARNESS_CC", o, w);
#endif

  // A compiler that cannot be started
  o = jit_base();
  o["compiler"] = "no_such_compiler_harness";
  bool ok2 = jit_try(src, d + SEP, o, w, err);
  check(!ok2, "missing compiler fails: " + err);
}

#ifdef _WIN32
// Compiler started through a batch file wrapper, without cmd.exe syntax in options
static void test_bat_wrapper() {
  section("JIT with a batch file as compiler (" + g_jit + ")");
  std::string wd = join(g_work, "wrap");
  make_dir(wd);
  std::string ccw = join(wd, "ccwrap.bat"), ldw = join(wd, "ldwrap.bat");
  if (g_jit=="cl") {
    std::string bat = Peek::vcvarsall();
#if defined(_M_ARM64) || defined(__aarch64__)
    std::string arch = "arm64";
#elif defined(_WIN64)
    std::string arch = "x64";
#else
    std::string arch = "x86";
#endif
    write_file(ccw, "@call \"" + bat + "\" " + arch + " >nul\r\n@cl %*\r\n");
    write_file(ldw, "@call \"" + bat + "\" " + arch + " >nul\r\n@link %*\r\n");
  } else {
    write_file(ccw, "@gcc %*\r\n");
    write_file(ldw, "@gcc %*\r\n");
  }
  // A root without spaces, so that CasADi's paths are not quoted in direct mode
  const char* rt = getenv("RUNNER_TEMP");
  std::string root = join(rt ? rt : g_work, "wrapdirs_" + g_jit);
  make_dir(root);
  std::string cwd = current_dir();
  for (const std::string& n : {"plain", "with space", "R&md,sfxdir3&b", "x)&md,sfxdir3&(y",
      "c^&md,sfxdir2", "pct%HARNESS_V%dir", "100%"}) {
    std::string d = join(root, n);
    make_dir(d);
    change_dir(d);
    std::string src = join(d, "jitsrc.c");
    write_file(src, jit_source);
    Dict o = jit_base();
    o["compiler"] = ccw;
    o["linker"] = ldw;
    std::string err;
    bool ok = jit_try(src, d + SEP, o, Want(), err);
    std::vector<std::string> s = sentinels(d);
    for (const char* m : {"sfxdir", "sfxdir2", "sfxdir3"}) {
      if (exists(join(root, m))) s.push_back(join(root, m));
    }
    check(s.empty(), "batch compiler, dir=" + esc(n) + ": side effect " + vesc(s));
    if (std::string(n).find('%')!=std::string::npos) {
      // cmd.exe expands %VAR% in the paths it is given
      if (!ok) note("batch compiler, dir=" + esc(n) + ": " + err.substr(0, 200));
      change_dir(cwd);
      continue;
    }
    check(ok, "batch compiler, dir=" + esc(n) + ": " + err);
    change_dir(cwd);
  }
}
#endif

static void test_jit_threads() {
  section("JIT from several threads (" + g_jit + ")");
#ifdef CASADI_WITH_THREAD
  const int n = 6;
  std::vector<std::string> errs(n);
  std::vector<int> oks(n, 0);
  std::vector<std::thread> th;
  for (int i=0; i<n; ++i) {
    std::string d = join(g_work, "thread " + std::to_string(i) + " & x");
    make_dir(d);
    std::string src = join(d, "jit src.c");
    write_file(src, jit_source);
    th.emplace_back([=, &errs, &oks]() {
      std::string err;
      oks[i] = jit_try(src, d + SEP, jit_base(), Want(), err);
      errs[i] = err;
    });
  }
  for (auto& t : th) t.join();
  for (int i=0; i<n; ++i) check(oks[i], "thread " + std::to_string(i) + ": " + errs[i]);
#else
  note("built without CASADI_WITH_THREAD");
#endif
}

// ---------------------------------------------------------------- main
int main(int argc, char** argv) {
  std::vector<std::string> only;
  std::string probe;
  for (int i=1; i<argc; ++i) {
    std::string a = argv[i];
    if (a=="--work" && i+1<argc) {
      g_work = argv[++i];
    } else if (a=="--dump" && i+1<argc) {
      std::string v = argv[++i];
      size_t eq = v.find('=');
      g_dumps.push_back({v.substr(0, eq), v.substr(eq+1)});
    } else if (a=="--jit" && i+1<argc) {
      g_jit = argv[++i];
    } else if (a=="--only" && i+1<argc) {
      only.push_back(argv[++i]);
    } else if (a=="--console-probe" && i+1<argc) {
      probe = argv[++i];
    } else {
      std::cerr << "unknown argument " << a << std::endl;
      return 2;
    }
  }
  if (g_work.empty()) g_work = current_dir();
  make_dirs(g_work);
  g_outfile = join(g_work, probe.empty() ? "adump.txt" : "adump_probe.txt");
#ifdef _WIN32
  if (!probe.empty()) return console_probe(probe);
#endif
  auto want = [&](const std::string& s) {
    if (only.empty()) return s!="threads";
    for (const std::string& o : only) if (o==s) return true;
    return false;
  };
  std::cout << "work=" << g_work << " jit=" << g_jit;
#ifdef CASADI_WITH_THREAD
  std::cout << " CASADI_WITH_THREAD";
#endif
#ifdef _WIN32
  std::cout << " ACP=" << GetACP();
#endif
  std::cout << std::endl;

  if (want("threads") && g_jit!="none") test_jit_threads();
  if (want("quote")) test_quote_roundtrip();
#ifdef _WIN32
  if (want("units")) test_units();
  if (want("winrun")) test_win_run();
  if (want("console")) test_no_console();
  if (want("vs") && g_jit=="cl") test_vs_setup();
#endif
  if (want("jit") && g_jit!="none") {
    jit_dirs("JIT in awkward directories, direct (" + g_jit + ")", "jitdirs", Dict(), false);
#ifdef _WIN32
    Dict sh;
    sh["flags"] = std::vector<std::string>{">nul"};
    jit_dirs("JIT in awkward directories, cmd.exe (" + g_jit + ")", "shelldirs", sh, true);
#endif
    test_jit_options();
  }
#ifdef _WIN32
  if (want("wrap") && g_jit!="none") test_bat_wrapper();
#endif
  std::cout << "SUMMARY pass=" << n_pass << " fail=" << n_fail << " note=" << n_note << std::endl;
  return n_fail ? 1 : 0;
}
