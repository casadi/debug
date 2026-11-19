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
  using ShellCompiler::vcvarsall_;
  using ShellCompiler::vcvarsall_searched_;
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
static Dump dump_line(const std::string& line) {
  remove(g_outfile.c_str());
  set_env("ADUMP_OUT", g_outfile);
#ifdef _WIN32
  int ret = casadi::run(line, nullptr);
#else
  int ret = casadi::run(line);
#endif
  return parse_dump(read_file(g_outfile), ret);
}

// Run a command line through system(), as the plugin did before
static Dump dump_system(const std::string& line) {
  remove(g_outfile.c_str());
  set_env("ADUMP_OUT", g_outfile);
  int ret = system(line.c_str());
  return parse_dump(read_file(g_outfile), ret);
}

// Quote for a CreateProcess command line (MSVC runtime rules)
static std::string aquote(const std::string& s) {
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
}

#ifdef _WIN32
// Full path of a program on the PATH of an environment
static bool find_full(const std::string& prog, const casadi::WinEnv* env, std::string& full) {
  char buf[MAX_PATH];
  std::string path;
  if (env) path = std::string(env->path.begin(), env->path.end());
  DWORD n = SearchPathA(env ? path.c_str() : nullptr, prog.c_str(), ".exe", MAX_PATH, buf, nullptr);
  if (n==0 || n>=MAX_PATH) return false;
  full = buf;
  return true;
}
#endif

static bool non_ascii(const std::string& s) {
  for (char c : s) if (static_cast<unsigned char>(c)>=128) return true;
  return false;
}

// Does quote() refuse the string?
static bool refuses(const std::string& s, std::string* msg=nullptr) {
  try {
    casadi::quote(s);
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
  section("quote: argv round trip through the shell");
  for (const DumpExe& d : g_dumps) {
    for (const std::string& s : corpus()) {
      std::string msg;
      bool refused = refuses(s, &msg);
#ifdef _WIN32
      bool should = s.find_first_of("%!\"")!=std::string::npos;
      for (char c : s) if (static_cast<unsigned char>(c)<0x20) should = true;
      if (non_ascii(s) && !casadi::win_bytes_safe()) should = true;
      check(refused==should, "quote(" + esc(s) + ") refused=" + std::to_string(refused));
      if (refused) continue;
      // cmd.exe cannot pass these through a /c line
      if (s.find_first_of("\r\n")!=std::string::npos) continue;
#endif
      if (refused) {
        fail("quote(" + esc(s) + ") refused: " + msg);
        continue;
      }
      std::string q = casadi::quote(s);
      // cmd.exe lines are limited to 8191 characters
      bool once = s.size()>3000;
      std::string line = casadi::quote(d.path) + " " + q + (once ? "" : " sep " + q);
      Dump r = dump_line(line);
      std::string ws = !s.empty() && (s[0]=='-' || s[0]=='@') ? "./" + s : s;
      std::vector<std::string> want = {ws, "sep", ws};
      if (once) want = {ws};
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

#ifdef _WIN32
// ---------------------------------------------------------------- Windows unit tests
static void test_units() {
  set_env("HARNESS_V", "hello world");
  set_env("HARNESS_OP", "a&b");
  unset_env("HARNESS_UNDEF");
  section("cmd_balanced");
  struct C { std::string s; bool want; };
  std::vector<C> cs = {
    {"", true}, {"cl.exe", true}, {">nul", true}, {"a^&b", true}, {"-DX=\"a b\"", true},
    {"\"", false}, {"a\"b", false}, {"^", false}, {"a^", false}, {"\"a^\"", true},
    {"-DX=^\"y\"", false}, {"-DS=^\"v1^\"", true}, {"/Fo^", false}, {"^^", true}, {"\"^\"^", false},
    {"\"a\"\"b\"", true}, {"^\"", true},
  };
  for (const C& c : cs) {
    bool got = casadi::cmd_balanced(casadi::win_wide(c.s, false));
    check(got==c.want, "cmd_balanced(" + esc(c.s) + ")=" + std::to_string(got));
  }

  section("program");
  struct P { std::string cmd, prog; };
  std::vector<P> ps = {
    {"cl.exe /c", "cl.exe"}, {"cl", "cl"},
    {"\"C:\\Program Files\\VS\\cl.exe\" /c", "C:\\Program Files\\VS\\cl.exe"},
    {"\"C:\\Program Files\\VS\\cl.exe\"", "C:\\Program Files\\VS\\cl.exe"},
    {"\"cl.exe\"", "cl.exe"}, {"cl\t/c", "cl"}, {"\"unterminated", "unterminated"}, {"\"\"", ""},
    {" cl.exe /c", "cl.exe"}, {"\t \"cl.exe\" /c", "cl.exe"}, {"C:\\a b\\cc.bat x", "C:\\a"},
  };
  for (const P& p : ps) {
    std::string prog = casadi::program(p.cmd);
    check(prog==p.prog, "program(" + esc(p.cmd) + ") gave " + esc(prog) + ", expected " + esc(p.prog));
  }

  section("include_dir");
  struct I { std::string in, want; };
  std::vector<I> is = {
    {"C:\\a b", "C:\\a b"}, {"\"C:\\a b\"", "C:\\a b"}, {"%HARNESS_V%", "hello world"},
    {"\"%HARNESS_V%\"", "hello world"}, {"\"a\"b\"", "a\"b"}, {"\"", "\""}, {"\"\"", ""},
  };
  for (const I& c : is) {
    std::string got = casadi::include_dir(c.in);
    check(got==c.want, "include_dir(" + esc(c.in) + ")=" + esc(got) + ", expected " + esc(c.want));
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
    std::string got;
    if (!casadi::win_narrow(casadi::win_wide(e.in, true), got)) got = "(lossy)";
    check(got==e.want, "expand(" + esc(e.in) + ")=" + esc(got) + ", expected " + esc(e.want));
  }
}

// ---------------------------------------------------------------- win_run
static void test_win_run() {
  section("win_run");
  if (g_dumps.empty()) return;
  const std::string& exe = g_dumps[0].path;
  std::string qexe = casadi::quote(exe);
  for (int code : {0, 1, 3, 42, 255, 1000}) {
    set_env("ADUMP_EXIT", std::to_string(code));
    Dump r = dump_line(qexe);
    check(r.ret==code, "exit code " + std::to_string(code) + " gave " + std::to_string(r.ret));
    Dump v = dump_system(qexe);
    check(v.ret==r.ret, "exit code " + std::to_string(code) + ": system() gave " + std::to_string(v.ret));
  }
  unset_env("ADUMP_EXIT");
  {
    int ret = casadi::win_run("no_such_program_harness.exe x", nullptr);
    check(ret==-1, "win_run: nonexistent program gave " + std::to_string(ret));
    ret = casadi::run("no_such_program_harness.exe x", nullptr);
    int v = system("no_such_program_harness.exe x");
    check(ret!=0 && ret==v, "run: nonexistent program gave " + std::to_string(ret)
      + ", system() " + std::to_string(v));
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
    std::string d = join(g_work, "dir with space & (x) ^ y");
    make_dir(d);
    std::string copy = join(d, "argv dump.exe");
    CopyFileA(exe.c_str(), copy.c_str(), FALSE);
    Dump r = dump_line(casadi::quote(copy) + " \"a b\"");
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
  Dump r = dump_line(casadi::quote(g_dumps[0].path) + " probe");
  ss << "child_ret " << r.ret << "\nchild_console " << r.console << "\n";
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

// Child mode: hidden console, stdout redirected by the parent
static int console_probe2() {
  std::cout << "probe2 console_cp " << GetConsoleCP() << " window " << (GetConsoleWindow() ? 1 : 0)
            << std::endl;
  unset_env("ADUMP_OUT");
  std::cout.flush();
  int ret = casadi::run(casadi::quote(g_dumps[0].path) + " fromshell", nullptr);
  std::cout << "probe2 shell_ret " << ret << std::endl;
  if (g_jit!="none") {
    std::string d = join(g_work, "probe2 jit");
    make_dir(d);
    std::string src = join(d, "jit src.c");
    write_file(src, std::string("#pragma message(\"harnessmsg\")\n") + jit_source_probe());
    std::string err;
    std::cout.flush();
    bool ok = jit_try_probe(src, d + SEP, err);
    std::cout << "probe2 jit " << (ok ? "ok" : err) << std::endl;
  }
  std::cout << "probe2_done" << std::endl;
  return 0;
}

static void test_no_console() {
  section("win_run: parent without console");
  char self[MAX_PATH];
  GetModuleFileNameA(nullptr, self, MAX_PATH);
  std::string out = join(g_work, "probe.txt");
  remove(out.c_str());
  std::string cmd = aquote(self) + " --work " + aquote(g_work)
    + " --dump " + aquote(g_dumps[0].label + "=" + g_dumps[0].path) + " --jit "
    + g_jit + " --console-probe " + aquote(out);
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
  check(res.find("child_console 2")==std::string::npos, "child got a visible console window");
  if (g_jit!="none") check(res.find("probe_jit ok")!=std::string::npos, "JIT from a console-less parent");

  // Parent with a hidden console, or the inherited console, and redirected stdout
  for (int hidden=1; hidden>=0; --hidden) {
    std::string kind = hidden ? "hidden console" : "inherited console";
    section("win_run: parent with " + kind);
    std::string log = join(g_work, "probe2 log " + std::to_string(hidden) + ".txt");
    SECURITY_ATTRIBUTES sa = {sizeof(SECURITY_ATTRIBUTES), nullptr, TRUE};
    HANDLE f = CreateFileA(log.c_str(), GENERIC_WRITE, FILE_SHARE_READ, &sa, CREATE_ALWAYS,
      FILE_ATTRIBUTE_NORMAL, nullptr);
    std::string cmd2 = aquote(self) + " --work " + aquote(join(g_work, "p2_" + std::to_string(hidden)))
      + " --dump " + aquote(g_dumps[0].label + "=" + g_dumps[0].path)
      + " --jit " + g_jit + " --console-probe2 x";
    std::vector<char> buf2(cmd2.begin(), cmd2.end());
    buf2.push_back('\0');
    ZeroMemory(&si, sizeof(si));
    si.cb = sizeof(si);
    si.dwFlags = STARTF_USESTDHANDLES;
    si.hStdInput = GetStdHandle(STD_INPUT_HANDLE);
    si.hStdOutput = f;
    si.hStdError = f;
    ok = CreateProcessA(nullptr, buf2.data(), nullptr, nullptr, TRUE,
      hidden ? CREATE_NO_WINDOW : 0, nullptr, nullptr, &si, &pi);
    check(ok, "start probe, " + kind);
    if (!ok) continue;
    WaitForSingleObject(pi.hProcess, 180000);
    CloseHandle(pi.hProcess);
    CloseHandle(pi.hThread);
    CloseHandle(f);
    res = read_file(log);
    std::cout << res.substr(0, 3000);
    check(res.find("probe2_done")!=std::string::npos, kind + ": probe finished");
    check(res.find("a 66726f6d7368656c6c")!=std::string::npos,
      kind + ": child output reaches the parent's redirected stdout");
    if (g_jit!="none") {
      check(res.find("probe2 jit ok")!=std::string::npos, kind + ": JIT");
      check(res.find("harnessmsg")!=std::string::npos, kind + ": compiler output reaches the parent's stdout");
    }
  }
}

static void test_vs_setup() {
  section("Visual Studio setup");
  std::string full;
  bool on_path = find_full("cl.exe", nullptr, full);
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
    if (!env) {
      fail("win_env(" + a + ") gave no environment");
      continue;
    }
    std::string cl, link;
    bool f1 = find_full("cl.exe", env, cl), f2 = find_full("link.exe", env, link);
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
    check(env && find_full("cl.exe", env, cl), "batch file in " + esc(d) + ": cl.exe resolved");
  } catch (std::exception& e) {
    fail("batch file in " + esc(d) + " threw: " + e.what());
  }
  std::string d3 = join(g_work, "vs copy 100%");
  make_dir(d3);
  try {
    const casadi::WinEnv* e3 = casadi::win_env(join(d3, "x.bat"), archs[0]);
    check(e3==nullptr, "batch file in " + esc(d3) + ": expected no environment");
  } catch (std::exception& e) {
    fail("batch file in " + esc(d3) + " threw: " + e.what());
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
    "fw\xEF\xBC\x82q", "fw\xEF\xBC\x85p\xEF\xBC\x85", "tab\tx", "a\xCA\xBA" "b",
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
    // Trailing separator in the second pass
    o["include_dirs"] = std::vector<std::string>{shell ? d + SEP : d};
    std::string err;
    Want w;
    w.square = 1049;
    // Source compiled from cwd=d, so the quoted include also works without include_dirs
    write_file(join(d, "harness_inc.h"), "#define INC_VALUE 1000\n");
    bool ok = jit_try(src, d + SEP, o, w, err);
    std::vector<std::string> s = sentinels(d);
    check(s.empty(), "dir=" + esc(n) + ": side effect " + vesc(s));
#ifdef _WIN32
    bool expect_refusal = n.find_first_of("%!\t")!=std::string::npos
      || (non_ascii(n) && !casadi::win_bytes_safe());
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
#ifdef _WIN32
    if (std::string(n).find('%')!=std::string::npos) {
      bool ok = jit_try(src, d + SEP, o, wi, err);
      check(!ok && err.find("Unsupported character")!=std::string::npos, "include_dirs " + esc(n)
        + ": expected 'Unsupported character', got " + (ok ? std::string("success") : err));
      continue;
    }
#endif
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
      find_full("cl.exe", env, clfull);
      find_full("link.exe", env, linkfull);
      // Quoted full path with spaces
      o = jit_base();
      o["compiler"] = "\"" + clfull + "\"";
      if (!jit_try(src, d + SEP, o, w, err)) {
        note("compiler quoted full path (no Visual Studio setup): " + err.substr(0, 160));
      }
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


// Mode choice, include directories, unbalanced quotes
static void test_modes() {
  section("JIT mode choice and include_dirs (" + g_jit + ")");
  std::string d = join(g_work, "modes dir");
  make_dir(d);
  std::string src = join(d, "jit src.c");
  write_file(src, jit_source);
  std::string cwd = current_dir();
  std::string err;
  Want w, wi;
  wi.square = 1049;
  auto run = [&](const std::string& what, Dict o, const Want& want) {
    bool ok = jit_try(src, d + SEP, o, want, err);
    check(ok, what + ": " + err);
    std::vector<std::string> s = sentinels(d);
    check(s.empty(), what + ": side effect " + vesc(s));
  };
  auto refused = [&](const std::string& what, Dict o, const std::string& msg) {
    bool ok = jit_try(src, d + SEP, o, w, err);
    check(!ok && err.find(msg)!=std::string::npos, what + ": expected '" + msg + "', got "
      + (ok ? std::string("success") : err));
    std::vector<std::string> s = sentinels(d);
    check(s.empty(), what + ": side effect " + vesc(s));
  };

  // Include directories
  std::vector<std::string> names = {"inc space", "inc & (x)", "inc O'Brien", "inc 100%",
    "inc ^caret", "inc semi;colon", "inc \xC3\xA9t\xC3\xA9", "inc\xE2\x88\x96"};
#ifndef _WIN32
  names.push_back("inc $(mkdir sfxdir)");
  names.push_back("inc `mkdir sfxdir2`");
  names.push_back("inc \"dq\"");
  names.push_back("inc $HOME");
#endif
  for (const std::string& n : names) {
    std::string inc = join(d, n);
    make_dir(inc);
    write_file(join(inc, "harness_inc.h"), "#define INC_VALUE 1000\n");
    set_env("HARNESS_INC", inc);
    std::vector<std::string> forms;
    // include_dirs expand variables, so a literal $ or %VAR% in a name needs quoting
    if (n.find('$')==std::string::npos) forms.push_back(inc);
#ifdef _WIN32
    forms.push_back("\"" + inc + "\"");
    forms.push_back("%HARNESS_INC%");
    forms.push_back("\"%HARNESS_INC%\"");
#else
    if (inc.find('\'')==std::string::npos) forms.push_back("'" + inc + "'");
    forms.push_back("$HARNESS_INC");
    forms.push_back("${HARNESS_INC}");
    forms.push_back("\"$HARNESS_INC\"");
#endif
    for (int shell=0; shell<2; ++shell) {
#ifndef _WIN32
      if (shell) continue;
#endif
      for (const std::string& f : forms) {
        Dict o = jit_base();
        std::vector<std::string> fl = {def_flag() + "USE_INC"};
        if (shell) fl.push_back(">nul");
        o["flags"] = fl;
        o["include_dirs"] = std::vector<std::string>{f};
        std::string what = std::string(shell ? "cmd.exe" : "direct") + " include_dirs " + esc(f)
          + " (" + esc(n) + ")";
#ifdef _WIN32
        if (inc.find_first_of("%!")!=std::string::npos
            || (non_ascii(inc) && !casadi::win_bytes_safe())) {
          refused(what, o, "Unsupported character");
          continue;
        }
#endif
        run(what, o, wi);
      }
    }
  }

#ifdef _WIN32
  // Undefined variable as include directory
  {
    unset_env("HARNESS_UNDEF_INC");
    Dict o = jit_base();
    o["include_dirs"] = std::vector<std::string>{"%HARNESS_UNDEF_INC%"};
    refused("include_dirs [\"%HARNESS_UNDEF_INC%\"] (undefined)", o, "Unsupported character");
  }
  // Quoted non-ASCII option, in a directory with &
  {
    std::string ud = join(d, "u8&md,sfxdir&b");
    make_dir(ud);
    std::string usrc = join(ud, "jitsrc.c");
    write_file(usrc, jit_source);
    change_dir(ud);
    Dict o = jit_base();
    o["flags"] = std::vector<std::string>{def_flag() + "T=\"\xE3\x83\x86\xE3\x82\xB9\xE3\x83\x88\" "
      + def_flag() + "U=\"x\""};
    bool ok = jit_try(usrc, ud + SEP, o, w, err);
    std::vector<std::string> x = sentinels(ud);
    if (exists(join(ud, "sfxdir")) || exists(join(d, "sfxdir"))) x.push_back("sfxdir");
    check(x.empty(), "quoted non-ASCII option: side effect " + vesc(x));
    note("quoted non-ASCII option (ACP " + std::to_string(GetACP()) + "): "
      + (ok ? std::string("ok") : err.substr(0, 160)));
    change_dir(cwd);
  }

  // Fallback triggers
  Dict o = jit_base();
  o["flags"] = std::vector<std::string>{">nul"};
  run("flags [\">nul\"]", o, w);

  // %CD% is not quoted, so use a directory without spaces
  const char* rt = getenv("RUNNER_TEMP");
  std::string cdd = join(rt ? rt : g_work, "cdtest_" + g_jit);
  make_dir(cdd);
  std::string inc = join(cdd, "cdinc");
  make_dir(inc);
  write_file(join(inc, "harness_inc.h"), "#define INC_VALUE 1000\n");
  change_dir(cdd);
  o = jit_base();
  o["flags"] = std::vector<std::string>{def_flag() + "USE_INC", inc_flag() + "%CD%\\cdinc"};
  run("flags [\"" + inc_flag() + "%CD%\\cdinc\"]", o, wi);
  o["flags"] = std::vector<std::string>{def_flag() + "USE_INC", "\"" + inc_flag() + "%CD%\\cdinc\""};
  run("flags [\"\\\"" + inc_flag() + "%CD%\\cdinc\\\"\"]", o, wi);
  change_dir(cwd);

  // Unbalanced quotes with cmd.exe
  o = jit_base();
  o["flags"] = std::vector<std::string>{def_flag() + "Q=\\\"", ">nul"};
  refused("flags [" + esc(def_flag() + "Q=\\\"") + ", \">nul\"]", o, "Unbalanced quotes");
  // No operator outside quotes: started directly, nothing interprets the rest
  o["flags"] = std::vector<std::string>{"-DQ=\"a&mkdir sfxdir&b"};
  {
    bool ok = jit_try(src, d + SEP, o, w, err);
    std::vector<std::string> sx = sentinels(d);
    check(sx.empty(), "flags [" + esc("-DQ=\"a&mkdir sfxdir&b") + "]: side effect " + vesc(sx));
    note("flags [" + esc("-DQ=\"a&mkdir sfxdir&b") + "]: " + (ok ? "ok" : err.substr(0, 160)));
  }

  // Compiler through a variable holding cmd.exe syntax
  std::string sd = join(d, "setup & (x)");
  make_dir(sd);
  std::string sb = join(sd, "setup.bat");
  if (g_jit=="cl") {
    std::string bat = Peek::vcvarsall();
#if defined(_M_ARM64) || defined(__aarch64__)
    std::string arch = "arm64";
#elif defined(_WIN64)
    std::string arch = "x64";
#else
    std::string arch = "x86";
#endif
    write_file(sb, "@call \"" + bat + "\" " + arch + " >nul\r\n");
    set_env("HARNESS_CC", "call \"" + sb + "\" && cl");
    set_env("HARNESS_LD", "call \"" + sb + "\" && link");
  } else {
    write_file(sb, "@set HARNESS_FROM_BAT=1\r\n");
    set_env("HARNESS_CC", "call \"" + sb + "\" && gcc");
    set_env("HARNESS_LD", "call \"" + sb + "\" && gcc");
  }
  o = jit_base();
  o["compiler"] = "%HARNESS_CC%";
  o["linker"] = "%HARNESS_LD%";
  run("compiler \"%HARNESS_CC%\" = " + esc(getenv("HARNESS_CC")), o, w);

  // Program names with a dot, found by cmd.exe through PATHEXT
  if (g_jit=="gcc") {
    std::string gcc;
    if (find_full("gcc", nullptr, gcc)) {
      std::string dir = gcc.substr(0, gcc.find_last_of('\\'));
      WIN32_FIND_DATAA fd;
      HANDLE h = FindFirstFileA((dir + "\\*-gcc-*.exe").c_str(), &fd);
      if (h!=INVALID_HANDLE_VALUE) {
        std::string name = fd.cFileName;
        FindClose(h);
        name = name.substr(0, name.size()-4);
        std::cout << "   dotted gcc: " << name << std::endl;
        o = jit_base();
        o["compiler"] = name;
        o["linker"] = name;
        run("compiler " + esc(name), o, w);
      } else {
        note("no *-gcc-*.exe next to " + gcc);
      }
    }
  }
#else
  Dict o = jit_base();
  o["flags"] = std::vector<std::string>{"2>/dev/null"};
  run("flags [\"2>/dev/null\"]", o, w);
#endif
}

#ifdef _WIN32
static std::string native_arch() {
#if defined(_M_ARM64) || defined(__aarch64__)
  return "arm64";
#elif defined(_WIN64)
  return "x64";
#else
  return "x86";
#endif
}

// Caret handling, PATHEXT, plain paths, %CD%
static void test_cmd_details() {
  section("cmd.exe details (" + g_jit + ")");
  const char* rt = getenv("RUNNER_TEMP");
  std::string root = join(rt ? rt : g_work, "details_" + g_jit);
  make_dir(root);
  std::string cwd = current_dir();
  std::string err;
  Want w;
  auto in_dir = [&](const std::string& n) {
    std::string d = join(root, n);
    make_dir(d);
    change_dir(d);
    std::string src = join(d, "jitsrc.c");
    write_file(src, jit_source);
    return d;
  };
  auto sfx = [&](const std::string& d) {
    std::vector<std::string> r = sentinels(d);
    for (const char* m : {"sfxdirQ", "sfxdirC", "sfxdir"}) {
      if (exists(join(d, m))) r.push_back(join(d, m));
      if (exists(join(root, m))) r.push_back(join(root, m));
    }
    return r;
  };
  auto expect_refusal = [&](const std::string& what, const std::string& d, Dict o,
      const std::string& msg) {
    bool ok = jit_try(join(d, "jitsrc.c"), d + SEP, o, w, err);
    check(!ok && err.find(msg)!=std::string::npos, what + ": expected '" + msg + "', got "
      + (ok ? std::string("success") : err.substr(0, 200)));
    std::vector<std::string> x = sfx(d);
    check(x.empty(), what + ": side effect " + vesc(x));
  };
  auto expect_ok = [&](const std::string& what, const std::string& d, Dict o, const Want& want) {
    bool ok = jit_try(join(d, "jitsrc.c"), d + SEP, o, want, err);
    check(ok, what + ": " + err.substr(0, 300));
    std::vector<std::string> x = sfx(d);
    check(x.empty(), what + ": side effect " + vesc(x));
  };

  // 1. Carets and quotes
  set_env("HARNESS_CQ", "a^\"b");
  set_env("HARNESS_OQ", "a\"b");
  for (const char* n : {"R&md,sfxdirQ&b x", "a&mkdir sfxdirQ&b"}) {
    std::string d = in_dir(n);
    std::string dn = esc(n);
    Dict o = jit_base();
    o["flags"] = std::vector<std::string>{def_flag() + "X=^\"y\""};
    expect_refusal("cwd " + dn + " flags [" + esc(def_flag() + "X=^\"y\"") + "]", d, o,
      "Unbalanced quotes or trailing ^");
    o = jit_base();
    o["flags"] = std::vector<std::string>{">nul"};
    o["compiler_output_flag"] = g_jit=="cl" ? "/Fo^" : "-o^";
    expect_refusal("cwd " + dn + " compiler_output_flag ending in ^", d, o,
      "Unbalanced quotes or trailing ^");
    o = jit_base();
    o["flags"] = std::vector<std::string>{"%HARNESS_CQ%", ">nul"};
    {
      // cmd.exe sees an escaped quote: accepted
      bool ok = jit_try(join(d, "jitsrc.c"), d + SEP, o, w, err);
      std::vector<std::string> x = sfx(d);
      check(x.empty(), "cwd " + dn + " flag %HARNESS_CQ%=a^\"b: side effect " + vesc(x));
      note("cwd " + dn + " flag %HARNESS_CQ%=a^\"b: " + (ok ? std::string("ok") : err.substr(0, 120)));
    }
    o = jit_base();
    o["flags"] = std::vector<std::string>{"%HARNESS_OQ%", ">nul"};
    expect_refusal("cwd " + dn + " flag %HARNESS_OQ%=a\"b", d, o, "Unbalanced quotes");
    o = jit_base();
    o["flags"] = std::vector<std::string>{def_flag() + "TAG=^\"v1^\""};
    Want wt;
    wt.tag = "v1";
    expect_ok("cwd " + dn + " flags [" + esc(def_flag() + "TAG=^\"v1^\"") + "]", d, o, wt);
    o = jit_base();
    o["flags"] = std::vector<std::string>{def_flag() + "TAG=a^&b"};
    wt.tag = "a&b";
    expect_ok("cwd " + dn + " flags [" + esc(def_flag() + "TAG=a^&b") + "]", d, o, wt);
    change_dir(cwd);
  }

  // 2. PATHEXT: a script with the compiler's name
  if (g_jit=="gcc") {
    std::string wrap = join(root, "pathext wrap");
    make_dir(wrap);
    std::string marker = join(root, "pathext_marker.txt");
    write_file(join(wrap, "gcc.bat"), "@echo ran>>\"" + marker + "\"\r\n@gcc.exe %*\r\n");
    const char* path = getenv("PATH");
    std::string old_path = path ? path : "";
    set_env("PATH", wrap + ";" + old_path);
    std::string d = in_dir("pathext");
    remove(marker.c_str());
    Dict o = jit_base();
    expect_ok("gcc.bat first on PATH", d, o, w);
    check(exists(marker), "gcc.bat first on PATH: wrapper ran");
    set_env("PATH", old_path);
    // In the current directory
    std::string d2 = in_dir("pathext cwd");
    write_file(join(d2, "gcc.bat"), "@echo ran>>\"" + marker + "\"\r\n@gcc.exe %*\r\n");
    remove(marker.c_str());
    expect_ok("gcc.bat in the current directory", d2, jit_base(), w);
    check(exists(marker), "gcc.bat in the current directory: wrapper ran");
    change_dir(cwd);
  }
  if (g_jit=="cl") {
    std::string wrap = join(root, "clwrap");
    make_dir(wrap);
    std::string marker = join(root, "cl_marker.txt");
    std::string bat = Peek::vcvarsall();
    write_file(join(wrap, "cl.bat"), "@echo ran>>\"" + marker + "\"\r\n@call \"" + bat + "\" "
      + native_arch() + " >nul\r\n@cl.exe %*\r\n");
    write_file(join(wrap, "lk.bat"), "@call \"" + bat + "\" " + native_arch()
      + " >nul\r\n@link.exe %*\r\n");
    const char* path = getenv("PATH");
    std::string old_path = path ? path : "";
    set_env("PATH", wrap + ";" + old_path);
    std::string d = in_dir("clbat");
    remove(marker.c_str());
    Dict o = jit_base();
    o["compiler"] = "cl";
    o["linker"] = "lk";
    expect_ok("user cl.bat on PATH", d, o, w);
    check(exists(marker), "user cl.bat on PATH: wrapper ran");
    set_env("PATH", old_path);
    change_dir(cwd);
  }

  // 4. Plain paths stay unquoted with cmd.exe
  {
    std::string wrap = join(root, "argwrap");
    make_dir(wrap);
    std::string args = join(root, "args.txt");
    std::string real = g_jit=="cl" ? "@call \"" + Peek::vcvarsall() + "\" " + native_arch()
      + " >nul\r\n@cl.exe %*\r\n" : "@gcc.exe %*\r\n";
    write_file(join(wrap, "argcc.bat"), "@echo %*>\"" + args + "\"\r\n@set SRC=%1\r\n" + real);
    std::string link = g_jit=="cl" ? "@call \"" + Peek::vcvarsall() + "\" " + native_arch()
      + " >nul\r\n@link.exe %*\r\n" : "@gcc.exe %*\r\n";
    write_file(join(wrap, "argld.bat"), link);
    for (const char* n : {"plainpath", "with space", "amp&x"}) {
      std::string d = in_dir(std::string("q_") + n);
      Dict o = jit_base();
      o["compiler"] = join(wrap, "argcc.bat");
      o["linker"] = join(wrap, "argld.bat");
      o["name"] = "plainobj";
      o["temp_suffix"] = false;
      expect_ok(std::string("argument capture, dir ") + esc(n), d, o, w);
      std::string got = read_file(args);
      while (!got.empty() && (got.back()=='\n' || got.back()=='\r' || got.back()==' ')) got.pop_back();
      std::string src = join(d, "jitsrc.c"), obj = join(d, "plainobj.obj");
      bool plain = std::string(n)=="plainpath";
      std::string q = plain ? "" : "\"";
      std::string want = (g_jit=="cl" ? "/c " : "-c ") + q + src + q + " "
        + (g_jit=="cl" ? "/Fo" : "-o ") + q + obj + q;
      check(got==want, std::string("argument capture, dir ") + esc(n) + ": got " + esc(got)
        + ", expected " + esc(want));
      if (plain) {
        // Same line through system()
        remove(args.c_str());
        int r = system((join(wrap, "argcc.bat") + " " + want).c_str());
        std::string v = read_file(args);
        while (!v.empty() && (v.back()=='\n' || v.back()=='\r' || v.back()==' ')) v.pop_back();
        check(r==0 && v==got, "argument capture: system() gave " + esc(v) + ", plugin " + esc(got));
      }
    }
    change_dir(cwd);
  }

  // 3. The compiler sees the same as with system()
  if (!g_dumps.empty()) {
    for (const char* n : {"vanilla", "vanilla space"}) {
      std::string d = in_dir(n);
      std::string src = join(d, "jitsrc.c");
      Dict o = jit_base();
      o["compiler"] = casadi::quote(g_dumps[0].path);
      o["linker"] = "cmd /c exit 1";
      o["name"] = "vobj";
      o["temp_suffix"] = false;
      remove(g_outfile.c_str());
      set_env("ADUMP_OUT", g_outfile);
      jit_try(src, d + SEP, o, w, err);
      Dump a = parse_dump(read_file(g_outfile), 0);
      // What the plugin used to run: compiler, flags, setup, source, output flag and object
      std::string setup = g_jit=="cl" ? "/c" : "-c", outf = g_jit=="cl" ? "/Fo" : "-o ";
      std::string obj = join(d, "vobj.obj");
      std::string vanilla = casadi::quote(g_dumps[0].path) + " " + setup + " " + src + " " + outf + obj;
      Dump b = dump_system(vanilla);
      if (std::string(n)=="vanilla") {
        check(a.ok && a.crt==b.crt && a.raw==b.raw, std::string("dir ") + esc(n)
          + ": argv " + vesc(a.crt) + " raw " + esc(a.raw) + " vs system() argv " + vesc(b.crt)
          + " raw " + esc(b.raw));
      } else {
        // system() splits the unquoted path; the plugin keeps it together
        check(a.ok && a.crt.size()>=3 && a.crt[1]==src, std::string("dir ") + esc(n) + ": argv "
          + vesc(a.crt));
        note(std::string("dir ") + esc(n) + ": system() gave " + vesc(b.crt));
      }
    }
    change_dir(cwd);
  }

  // Compile errors give a nonzero exit code
  {
    std::string d = in_dir("syntaxerr");
    std::string src = join(d, "jitsrc.c");
    write_file(src, "int f( {\n");
    bool ok = jit_try(src, d + SEP, jit_base(), w, err);
    check(!ok && err.find("Compilation failed")!=std::string::npos, "syntax error: "
      + (ok ? std::string("success") : err.substr(0, 120)));
    change_dir(cwd);
  }

  // CMDCMDLINE in options
  {
    std::string d = in_dir("cmdline&md,sfxdirC&x");
    for (const char* f : {"-DX=%CMDCMDLINE%", "-DX=%cmdcmdline:~0,5%", "-DX=!CMDCMDLINE!",
        "/I%=C:%inc", "-I%=c:%inc"}) {
      Dict o = jit_base();
      o["flags"] = std::vector<std::string>{f};
      expect_refusal("flags [" + esc(f) + "]", d, o, "Unsupported variable");
    }
    change_dir(cwd);
  }

  // 5. %CD% in options
  {
    for (const char* n : {"cdplain", "cd&md,sfxdirC&x", "cd(x)", "cd^x", "cd%x", "cd!x", "cd space"}) {
      std::string d = in_dir(n);
      make_dir(join(d, "inc"));
      write_file(join(join(d, "inc"), "harness_inc.h"), "#define INC_VALUE 1000\n");
      Want wi;
      wi.square = 1049;
      for (const char* form : {"%CD%", "%cd%", "%CD:~0%", "%CD:x=x%", "%__CD__%", "%cd:~0,200%",
          "%CD:\\=/%"}) {
        RemoveDirectoryA(join(d, "sfxdirC").c_str());
        RemoveDirectoryA(join(root, "sfxdirC").c_str());
        Dict o = jit_base();
        std::string f = form;
        std::string fl = inc_flag() + f + (f=="%__CD__%" ? "inc" : f=="%CD:\\=/%" ? "/inc" : "\\inc");
        o["flags"] = std::vector<std::string>{def_flag() + "USE_INC", fl};
        std::string what = "cwd " + esc(n) + " flags [" + esc(fl) + "]";
        bool awkward = std::string(n).find_first_of("&()^%!")!=std::string::npos;
        if (std::string(n)=="cd space") {
          bool ok = jit_try(join(d, "jitsrc.c"), d + SEP, o, wi, err);
          note(what + ": " + (ok ? "ok" : err.substr(0, 120)));
          continue;
        }
        if (awkward) {
          bool ok = jit_try(join(d, "jitsrc.c"), d + SEP, o, wi, err);
          std::vector<std::string> x = sfx(d);
          check(x.empty(), what + ": side effect " + vesc(x));
          check(!ok && err.find("current directory")!=std::string::npos, what
            + ": expected 'Unsupported character in the current directory', got "
            + (ok ? std::string("success") : err.substr(0, 200)));
        } else {
          expect_ok(what, d, o, wi);
        }
      }
    }
    change_dir(cwd);
  }
}

// Visual Studio setup that fails: JIT as without the setup
static void test_vs_fail() {
  section("Visual Studio setup failures (" + g_jit + ")");
  std::string real = Peek::vcvarsall();
  std::string d = join(g_work, "vsfail");
  make_dir(d);
  std::string src = join(d, "jit src.c");
  write_file(src, jit_source);
  struct B { std::string name, body; };
  std::vector<B> bs = {
    {"exit1.bat", "@exit /b 1\r\n"},
    {"nocl.bat", "@set HARNESS_FROM_BAT=1\r\n"},
    {"bogus.bat", "@call \"" + real + "\" bogus_arch\r\n"},
    {"empty.bat", ""},
  };
  for (const B& b : bs) {
    std::string bat = join(d, b.name);
    write_file(bat, b.body);
    Peek::vcvarsall_ = bat;
    Peek::vcvarsall_searched_ = true;
    std::string err;
    bool ok = jit_try(src, d + SEP, jit_base(), Want(), err);
    check(!ok && err.find("Environment setup failed")==std::string::npos
      && err.find("Compilation failed")!=std::string::npos,
      "vcvarsall " + b.name + ": expected 'Compilation failed', got " + (ok ? "success" : err.substr(0, 250)));
  }
  Peek::vcvarsall_ = real;
}

static void test_vs_time() {
  section("first JIT timing (" + g_jit + ")");
  std::string d = join(g_work, "vstime");
  make_dir(d);
  std::string src = join(d, "jit src.c");
  write_file(src, jit_source);
  std::string err;
  for (int i=0; i<2; ++i) {
    DWORD t0 = GetTickCount();
    bool ok = jit_try(src, d + SEP, jit_base(), Want(), err);
    DWORD t1 = GetTickCount();
    std::cout << "   JIT " << (i ? "second" : "first") << ": " << (t1-t0) << " ms" << std::endl;
    check(ok, "timed JIT: " + err);
  }
  DWORD t0 = GetTickCount();
  casadi::win_env(Peek::vcvarsall(), "x64_x86");
  std::cout << "   win_env x64_x86 alone: " << (GetTickCount()-t0) << " ms" << std::endl;
}
#endif

#ifdef _WIN32
// %CD% with a current directory that only wide APIs can name
static void test_cd_wide() {
  section("%CD% with non-ASCII current directories (" + g_jit + ", ACP " + std::to_string(GetACP()) + ")");
  const char* rt = getenv("RUNNER_TEMP");
  std::string root = join(rt ? rt : g_work, "cdwide_" + g_jit + "_" + std::to_string(GetACP()));
  make_dir(root);
  std::string src = join(root, "jitsrc.c");
  write_file(src, jit_source);
  std::string cwd = current_dir();
  std::wstring wroot(root.begin(), root.end());
  struct C { std::wstring name; const char* label; bool ok; };
  bool cp_single = GetACP()!=CP_UTF8;
  std::vector<C> cs = {
    {L"a\uFF02 -wrapper x \uFF02b", "a\\uFF02 -wrapper x \\uFF02b", false},
    {L"\u00E9t\u00E9", "\\u00E9t\\u00E9", cp_single},
    {L"plainwide", "plainwide", true},
  };
  for (const C& c : cs) {
    std::wstring d = wroot + L"\\" + c.name;
    CreateDirectoryW(d.c_str(), nullptr);
    CreateDirectoryW((d + L"\\inc").c_str(), nullptr);
    FILE* f = _wfopen((d + L"\\inc\\harness_inc.h").c_str(), L"wb");
    if (f) {
      fputs("#define INC_VALUE 1000\n", f);
      fclose(f);
    }
    if (!SetCurrentDirectoryW(d.c_str())) {
      fail(std::string("cannot enter ") + c.label);
      continue;
    }
    for (const std::string& form : {inc_flag() + "\"%CD%\\inc\"", inc_flag() + "%CD%\\inc"}) {
      Dict o = jit_base();
      o["flags"] = std::vector<std::string>{def_flag() + "USE_INC", form};
      Want wi;
      wi.square = 1049;
      std::string err;
      bool ok = jit_try(src, root + SEP, o, wi, err);
      std::string what = std::string("cwd ") + c.label + " flags [" + esc(form) + "]";
      if (c.ok) {
        check(ok, what + ": " + err.substr(0, 200));
      } else {
        check(!ok && err.find("current directory")!=std::string::npos, what
          + ": expected 'Unsupported character in the current directory', got "
          + (ok ? std::string("success") : err.substr(0, 200)));
      }
    }
    change_dir(cwd);
  }
  // What a narrow program receives for U+FF02 in this process
  if (!g_dumps.empty()) {
    std::wstring w = L"x\uFF02y";
    char nb[16] = {0};
    BOOL used = FALSE;
    WideCharToMultiByte(CP_ACP, 0, w.c_str(), -1, nb, sizeof(nb), nullptr,
      GetACP()==CP_UTF8 ? nullptr : &used);
    std::cout << "   U+FF02 with best fit in ACP " << GetACP() << ": " << esc(nb) << std::endl;
  }
}
#endif

#ifdef _WIN32
// Compiler started through a batch file wrapper, without cmd.exe syntax in options
static void test_bat_wrapper() {
  section("JIT with a batch file as compiler (" + g_jit + ")");
  std::string wd = join(g_work, "wrap");
  make_dir(wd);
  std::string arch;
#if defined(_M_ARM64) || defined(__aarch64__)
  arch = "arm64";
#elif defined(_WIN64)
  arch = "x64";
#else
  arch = "x86";
#endif
  std::string cc_body, ld_body;
  if (g_jit=="cl") {
    std::string bat = Peek::vcvarsall();
    cc_body = "@call \"" + bat + "\" " + arch + " >nul\r\n@cl %*\r\n";
    ld_body = "@call \"" + bat + "\" " + arch + " >nul\r\n@link %*\r\n";
  } else {
    cc_body = "@gcc %*\r\n";
    ld_body = "@gcc %*\r\n";
  }
  for (const char* ext : {".bat", ".cmd"}) {
    write_file(join(wd, std::string("ccwrap") + ext), cc_body);
    write_file(join(wd, std::string("ldwrap") + ext), ld_body);
  }
  const char* path = getenv("PATH");
  std::string old_path = path ? path : "";
  set_env("PATH", wd + ";" + old_path);
  struct V { std::string label, cc, ld; };
  std::vector<V> vs = {
    {"quoted .bat path", "\"" + join(wd, "ccwrap.bat") + "\"", "\"" + join(wd, "ldwrap.bat") + "\""},
    {".cmd path", "\"" + join(wd, "ccwrap.cmd") + "\"", "\"" + join(wd, "ldwrap.cmd") + "\""},
    {"no extension (PATHEXT)", "ccwrap", "ldwrap"},
    {".bat on PATH", "ccwrap.bat", "ldwrap.bat"},
  };
  std::string full;
  if (g_jit=="gcc") {
    vs.push_back({"cmd /c gcc", "cmd /c gcc", "cmd /c gcc"});
  } else if (find_full("cl.exe", nullptr, full)) {
    vs.push_back({"cmd /c cl", "cmd /c cl", "cmd /c link"});
  } else {
    vs.push_back({"cmd /c ccwrap", "cmd /c ccwrap", "cmd /c ldwrap"});
  }
  // A root without spaces, so that CasADi's paths are not quoted for being spaced
  const char* rt = getenv("RUNNER_TEMP");
  std::string root = join(rt ? rt : g_work, "wrapdirs_" + g_jit);
  make_dir(root);
  std::string cwd = current_dir();
  int k = 0;
  for (const V& v : vs) {
    std::string vroot = join(root, std::to_string(k++));
    make_dir(vroot);
    for (const char* n : {"plain", "with space", "R&D", "R&md,sfxdir3&b", "x)&md,sfxdir3&(y",
        "c^&x", "c^&md,sfxdir2", "(x86)", "O'Brien", "pct%HARNESS_V%dir", "100%", "bang!x!"}) {
      std::string d = join(vroot, n);
      make_dir(d);
      change_dir(d);
      std::string src = join(d, "jitsrc.c");
      write_file(src, jit_source);
      Dict o = jit_base();
      o["compiler"] = v.cc;
      o["linker"] = v.ld;
      std::string err;
      bool ok = jit_try(src, d + SEP, o, Want(), err);
      std::vector<std::string> s = sentinels(d);
      for (const char* m : {"sfxdir", "sfxdir2", "sfxdir3"}) {
        if (exists(join(vroot, m))) s.push_back(join(vroot, m));
      }
      std::string what = v.label + " (" + esc(v.cc) + "), dir=" + esc(n);
      check(s.empty(), what + ": side effect " + vesc(s));
      if (std::string(n).find_first_of("%!")!=std::string::npos) {
        check(!ok && err.find("Unsupported character")!=std::string::npos,
          what + ": expected 'Unsupported character', got " + (ok ? std::string("success") : err));
      } else {
        check(ok, what + ": " + err);
      }
      change_dir(cwd);
    }
  }
  set_env("PATH", old_path);
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
    } else if (a=="--console-probe2" && i+1<argc) {
      ++i;
      probe = "\x02";
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
  if (probe=="\x02") return console_probe2();
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
#ifndef _WIN32
  if (want("units")) {
    section("include_dir");
    set_env("HARNESS_V", "hello world");
    unset_env("HARNESS_UNDEF");
    struct I { std::string in, want; };
    std::vector<I> is = {
      {"/a b", "/a b"}, {"'/a b'", "/a b"}, {"'$HARNESS_V'", "$HARNESS_V"},
      {"\"$HARNESS_V\"", "hello world"}, {"$HARNESS_V", "hello world"},
      {"${HARNESS_V}x", "hello worldx"}, {"$HARNESS_Vx", ""}, {"$HARNESS_UNDEF", ""},
      {"a$", "a$"}, {"${HARNESS_V", "${HARNESS_V"}, {"$(x)", "$(x)"}, {"`x`", "`x`"},
      {"$1", ""}, {"'", "'"}, {"\"", "\""}, {"\"'a'\"", "'a'"},
    };
    for (const I& c : is) {
      std::string got = casadi::include_dir(c.in);
      check(got==c.want, "include_dir(" + esc(c.in) + ")=" + esc(got) + ", expected " + esc(c.want));
    }
  }
#endif
#ifdef _WIN32
  if (want("units")) test_units();
  if (want("winrun")) test_win_run();
  if (want("console")) test_no_console();
  if (want("vs") && g_jit=="cl") test_vs_setup();
#endif
  if (want("jit") && g_jit!="none") {
    jit_dirs("JIT in awkward directories (" + g_jit + ")", "jitdirs", Dict(), false);
    Dict sh;
#ifdef _WIN32
    sh["flags"] = std::vector<std::string>{">nul"};
#else
    sh["flags"] = std::vector<std::string>{">/dev/null"};
#endif
    jit_dirs("JIT in awkward directories, redirected, include dir with trailing separator ("
      + g_jit + ")", "shelldirs", sh, true);
    test_jit_options();
    test_modes();
#ifdef _WIN32
    test_cd_wide();
#endif
  }
#ifdef _WIN32
  if (want("wrap") && g_jit!="none") test_bat_wrapper();
  if (want("details") && g_jit!="none") test_cmd_details();
  if (only.size() && want("vsfail") && g_jit=="cl") test_vs_fail();
  if (only.size() && want("vstime") && g_jit=="cl") test_vs_time();
#endif
  std::cout << "SUMMARY pass=" << n_pass << " fail=" << n_fail << " note=" << n_note << std::endl;
  return n_fail ? 1 : 0;
}
