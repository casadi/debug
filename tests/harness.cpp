// Tests for the argument handling helpers of the shell compiler plugin.
// The plugin source is included verbatim; stubs/ provides the few casadi bits it needs.
#include "../src/shell_compiler.cpp"

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

// Expose the protected static helper
#ifdef _WIN32
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

// Start argv_dump through the plugin's run()
static Dump dump_args(const std::string& exe, const std::vector<std::string>& extra) {
  remove(g_outfile.c_str());
  set_env("ADUMP_OUT", g_outfile);
  std::vector<std::string> args = {exe};
  args.insert(args.end(), extra.begin(), extra.end());
#ifdef _WIN32
  int ret = casadi::run(args, nullptr);
#else
  int ret = casadi::run(args);
#endif
  return parse_dump(read_file(g_outfile), ret);
}

#ifdef _WIN32
// Start argv_dump with a raw command line tail
static Dump dump_raw(const std::string& exe, const std::string& tail) {
  remove(g_outfile.c_str());
  set_env("ADUMP_OUT", g_outfile);
  int ret = casadi::win_run(casadi::quote(exe) + " " + tail, nullptr);
  return parse_dump(read_file(g_outfile), ret);
}

// Start argv_dump through cmd.exe, as system() did
static Dump dump_cmd(const std::string& exe, const std::string& tail) {
  remove(g_outfile.c_str());
  set_env("ADUMP_OUT", g_outfile);
  const char* comspec = getenv("ComSpec");
  std::string cmd = casadi::quote(comspec ? comspec : "cmd.exe") + " /d /s /c \""
    + casadi::quote(exe) + " " + tail + "\"";
  int ret = casadi::win_run(cmd, nullptr);
  return parse_dump(read_file(g_outfile), ret);
}

static std::wstring widen(const std::string& s) {
  int n = MultiByteToWideChar(CP_ACP, 0, s.c_str(), -1, nullptr, 0);
  std::wstring w(n, L'\0');
  MultiByteToWideChar(CP_ACP, 0, s.c_str(), -1, &w[0], n);
  w.resize(n ? n-1 : 0);
  return w;
}

static std::string narrow(const std::wstring& w) {
  int n = WideCharToMultiByte(CP_ACP, 0, w.c_str(), -1, nullptr, 0, nullptr, nullptr);
  std::string s(n, '\0');
  WideCharToMultiByte(CP_ACP, 0, w.c_str(), -1, &s[0], n, nullptr, nullptr);
  s.resize(n ? n-1 : 0);
  return s;
}

// In-process reference: CommandLineToArgvW on "x <tail>"
static std::vector<std::string> ref_w(const std::string& tail) {
  int n = 0;
  std::wstring cl = L"x " + widen(tail);
  LPWSTR* w = CommandLineToArgvW(cl.c_str(), &n);
  std::vector<std::string> r;
  for (int i=1; i<n; ++i) r.push_back(narrow(w[i]));
  LocalFree(w);
  return r;
}
#endif

static bool non_ascii(const std::string& s) {
  for (char c : s) if (static_cast<unsigned char>(c)>=128) return true;
  return false;
}

// ---------------------------------------------------------------- corpus
static std::vector<std::string> corpus() {
  std::vector<std::string> c = {
    "", " ", "  ", "\t", "a", "a b", "a\tb", "a  b", " leading", "trailing ",
    "\"", "\"\"", "\"\"\"", "a\"b", "\"a b\"", "a\"\"b", "a\" \"b",
    "\\", "\\\\", "\\\\\\", "a\\", "a\\\\", "a b\\", "a b\\\\", "a b\\\\\\",
    "\\\"", "\\\\\"", "\\\\\\\"", "a\\\"b", "a\\\\\"b", "a\\\\\\\"b", "\\\\a", "a\\b\\c",
    "C:\\Program Files\\", "C:\\Program Files (x86)\\x\\", "\\\\server\\share\\a b\\",
    "&|<>^()%!'", "a&b", "a|b", "a^b", "a^&b", "a&&b", "a||b", "a>b", "a<b", "(a)",
    "%PATH%", "%%", "100%", "%", "!PATH!", "'single'", "'a b'", "`tick`", "$HOME", "$(x)",
    ";", "a;b", "*", "?", "~", "{a,b}", "#hash", "a=b", "@file", "-", "--",
    "\xC3\xA9t\xC3\xA9", "caf\xE9", "\xE4\xB8\xAD", "\x81\x8D\x8F\x90\x9D",
    "line\nbreak", "cr\rx", "vt\vx", "ff\fx",
    "-DTAG=\"a b\"", "/FoC:\\a b\\x.obj", "-I\"C:\\a b\\\"",
    "\"C:\\a b\\\\\" rest", "a\\\\\\\\\"\"b",
  };
  std::string long1(8000, 'a');
  long1 += " \" \\ \\\" end\\";
  c.push_back(long1);
  std::string long2;
  for (int i=0; i<2000; ++i) long2 += "x\\\" ";
  c.push_back(long2);
  return c;
}

// ---------------------------------------------------------------- quote round trip
static void test_quote_roundtrip() {
  section("quote: argv round trip through run()");
  for (const DumpExe& d : g_dumps) {
    for (const std::string& s : corpus()) {
      Dump r = dump_args(d.path, {s, "sep", s});
      std::vector<std::string> want = {s, "sep", s};
      std::string what = d.label + " input=" + esc(s) + " quoted=" + esc(casadi::quote(s));
      if (!r.ok) {
        fail(what + " child did not run, ret=" + std::to_string(r.ret));
        continue;
      }
      bool acp = non_ascii(s);
      if (r.crt==want) {
        pass();
      } else if (acp) {
        note("non-ASCII (code page) " + what + " crt=" + vesc(r.crt));
      } else {
        fail(what + " expected=" + vesc(want) + " crt=" + vesc(r.crt));
      }
#ifdef _WIN32
      if (r.w==want) {
        pass();
      } else if (acp) {
        note("non-ASCII (code page) CommandLineToArgvW " + what + " w=" + vesc(r.w));
      } else {
        fail("CommandLineToArgvW " + what + " expected=" + vesc(want) + " w=" + vesc(r.w));
      }
#endif
    }
  }
#ifdef _WIN32
  // Command line beyond the CreateProcess limit
  if (!g_dumps.empty()) {
    std::string big(33000, 'a');
    Dump r = dump_args(g_dumps[0].path, {big});
    if (r.ret==-1) {
      note("command line of 33000 chars: run() returns -1 (reported as 'Could not start compiler')");
    } else {
      note("command line of 33000 chars: run() returned " + std::to_string(r.ret));
    }
  }
#else
  // Display quoting is valid for sh
  section("quote: display quoting is valid sh");
  for (const DumpExe& d : g_dumps) {
    for (const std::string& s : corpus()) {
      if (s.size()>4000) continue;
      Dump r = dump_args("/bin/sh", {"-c", casadi::quote(d.path) + " " + casadi::quote(s)});
      check(r.ok && r.crt==std::vector<std::string>{s},
        d.label + " input=" + esc(s) + " quoted=" + esc(casadi::quote(s)) + " sh gave " + vesc(r.crt));
    }
  }
#endif
}

// ---------------------------------------------------------------- expand: splitting
#ifdef _WIN32
static std::string fuzz_string(std::mt19937& rng, int maxlen) {
  static const char alphabet[] = {'a', 'b', ' ', '\t', '"', '"', '\\', '\\', '\\'};
  std::uniform_int_distribution<int> len(0, maxlen), pick(0, sizeof(alphabet)-1);
  int n = len(rng);
  std::string s;
  for (int i=0; i<n; ++i) s += alphabet[pick(rng)];
  return s;
}

static bool has_cmd_meta(const std::string& s) {
  return s.find_first_of("&|<>^()%!\n\r") != std::string::npos;
}

static void test_expand_split() {
  section("expand: splitter vs child CRT / CommandLineToArgvW / cmd.exe");
  std::vector<std::string> inputs;
  for (const std::string& s : corpus()) if (s.find('%')==std::string::npos) inputs.push_back(s);
  std::mt19937 rng(12345);
  for (int i=0; i<400; ++i) inputs.push_back(fuzz_string(rng, 10));

  std::set<std::string> seen;
  for (const DumpExe& d : g_dumps) {
    int div = 0, div_cmd = 0;
    for (const std::string& s : inputs) {
      if (!seen.insert(d.label + "\x01" + s).second) continue;
      std::vector<std::string> e = casadi::expand(s);
      Dump r = dump_raw(d.path, s);
      bool acp = non_ascii(s);
      if (!r.ok) {
        fail(d.label + " child did not run for " + esc(s));
        continue;
      }
      if (r.crt==e) {
        pass();
      } else if (acp || d.label!="msvc") {
        note(d.label + " CRT differs from expand(): input=" + esc(s) + " expand=" + vesc(e)
          + " crt=" + vesc(r.crt));
        div++;
      } else {
        fail(d.label + " CRT differs from expand(): input=" + esc(s) + " expand=" + vesc(e)
          + " crt=" + vesc(r.crt));
        div++;
      }
      // What system() gave: cmd.exe, then the CRT
      if (!has_cmd_meta(s) && s.size()<4000) {
        Dump c = dump_cmd(d.path, s);
        if (c.crt==e) {
          pass();
        } else if (acp || d.label!="msvc") {
          note(d.label + " cmd.exe+CRT differs from expand(): input=" + esc(s) + " expand="
            + vesc(e) + " cmd=" + vesc(c.crt));
          div_cmd++;
        } else {
          fail(d.label + " cmd.exe+CRT differs from expand(): input=" + esc(s) + " expand="
            + vesc(e) + " cmd=" + vesc(c.crt));
          div_cmd++;
        }
      }
    }
    std::cout << "   " << d.label << ": " << div << " CRT and " << div_cmd
              << " cmd.exe divergences" << std::endl;
  }

  // Large in-process fuzz: expand() vs CommandLineToArgvW, and quote() round trip
  int div_w = 0, bad_rt = 0;
  std::mt19937 rng2(777);
  for (int i=0; i<100000; ++i) {
    std::string s = fuzz_string(rng2, 12);
    std::vector<std::string> e = casadi::expand(s);
    std::vector<std::string> w = ref_w(s);
    if (w!=e) {
      if (div_w<15) note("CommandLineToArgvW differs from expand(): input=" + esc(s)
        + " expand=" + vesc(e) + " w=" + vesc(w));
      div_w++;
    }
    std::vector<std::string> rt = ref_w(casadi::cmdline(e));
    std::vector<std::string> rt2 = casadi::expand(casadi::cmdline(e));
    if (rt!=e || rt2!=e) {
      if (bad_rt<15) fail("cmdline(expand(s)) does not round trip: input=" + esc(s)
        + " expand=" + vesc(e) + " cmdline=" + esc(casadi::cmdline(e)) + " w=" + vesc(rt)
        + " expand2=" + vesc(rt2));
      bad_rt++;
    } else {
      pass();
    }
  }
  std::cout << "   in-process fuzz: " << div_w << " CommandLineToArgvW divergences, "
            << bad_rt << " round trip failures out of 100000" << std::endl;
}

static void test_expand_vars() {
  section("expand: %VAR%");
  set_env("HARNESS_V", "hello world");
  set_env("HARNESS_Q", "say \"hi there\"");
  set_env("HARNESS_E", "");
  unset_env("HARNESS_UNDEF");
  struct Case { std::string in; std::vector<std::string> want; };
  std::vector<Case> cases = {
    {"%HARNESS_V%", {"hello", "world"}},
    {"\"%HARNESS_V%\"", {"hello world"}},
    {"%harness_v%", {"hello", "world"}},
    {"a%HARNESS_V%b", {"ahello", "worldb"}},
    {"%HARNESS_Q%", {"say", "hi there"}},
    {"%HARNESS_UNDEF%", {"%HARNESS_UNDEF%"}},
    {"%HARNESS_UNDEF", {"%HARNESS_UNDEF"}},
    {"%%", {"%%"}},
    {"%", {"%"}},
    {"100%", {"100%"}},
    {"50% off", {"50%", "off"}},
    {"%%HARNESS_V%%", {"%hello", "world%"}},
    {"%HARNESS_V%%HARNESS_V%", {"hello", "worldhello", "world"}},
    {"x%HARNESS_E%y", {"xy"}},
    {"%HARNESS_E%", {}},
    {"", {}},
    {"   ", {}},
  };
  for (const Case& c : cases) {
    std::vector<std::string> e = casadi::expand(c.in);
    check(e==c.want, "input=" + esc(c.in) + " expected=" + vesc(c.want) + " expand=" + vesc(e));
    // Reference: what cmd.exe /c and the CRT produce
    for (const DumpExe& d : g_dumps) {
      if (d.label!="msvc") continue;
      Dump r = dump_cmd(d.path, c.in);
      if (r.crt!=e) note("cmd.exe differs: input=" + esc(c.in) + " expand=" + vesc(e)
        + " cmd=" + vesc(r.crt));
      else pass();
    }
  }
  // Long expansion
  std::string big(20000, 'v');
  set_env("HARNESS_BIG", big);
  std::vector<std::string> e = casadi::expand("-D%HARNESS_BIG%");
  check(e.size()==1 && e[0]=="-D" + big, "20000-char expansion: got " + std::to_string(e.size())
    + " args, first size " + std::to_string(e.empty() ? 0 : e[0].size()));
}
#else
static void test_expand_posix() {
  section("expand: wordexp");
  set_env("HARNESS_V", "hello world");
  set_env("HARNESS_W", "$(mkdir sfxdir_w)");
  unset_env("HARNESS_UNDEF");
  std::string home = getenv("HOME") ? getenv("HOME") : "";
  struct Case { std::string in; std::vector<std::string> want; };
  std::vector<Case> ok = {
    {"a b", {"a", "b"}},
    {"'a b'", {"a b"}},
    {"\"a b\"", {"a b"}},
    {"a\\ b", {"a b"}},
    {"$HARNESS_V", {"hello", "world"}},
    {"\"$HARNESS_V\"", {"hello world"}},
    {"${HARNESS_V}x", {"hello", "worldx"}},
    {"$HARNESS_UNDEF", {}},
    {"x${HARNESS_UNDEF:-dflt}", {"xdflt"}},
    {"~", {home}},
    {"-DTAG=\"a b\"", {"-DTAG=a b"}},
    {"'a;b'", {"a;b"}},
    {"\"a|b\"", {"a|b"}},
    {"'$(x)'", {"$(x)"}},
    {"$HARNESS_W", {"$(mkdir", "sfxdir_w)"}},
    {"", {}},
  };
  for (const Case& c : ok) {
    try {
      std::vector<std::string> e = casadi::expand(c.in);
      check(e==c.want, "input=" + esc(c.in) + " expected=" + vesc(c.want) + " expand=" + vesc(e));
    } catch (std::exception& ex) {
      fail("input=" + esc(c.in) + " threw: " + ex.what());
    }
  }
  // Must be refused, with no side effect
  std::string dir = join(g_work, "wordexp");
  make_dir(dir);
  std::string cwd = current_dir();
  change_dir(dir);
  std::vector<std::string> bad = {
    "$(mkdir sfxdir1)", "`mkdir sfxdir2`", "\"$(mkdir sfxdir3)\"", "\"`mkdir sfxdir4`\"",
    "${HARNESS_UNDEF:-$(mkdir sfxdir5)}", "${HARNESS_UNDEF:-`mkdir sfxdir6`}",
    "$((`mkdir sfxdir7`))", "$(( $(mkdir sfxdir8) ))", "${HARNESS_V:+$(mkdir sfxdir9)}",
    "${HARNESS_UNDEF=$(mkdir sfxdir10)}", "x$(mkdir sfxdir11)y", "\"${HARNESS_UNDEF:-$(mkdir sfxdir12)}\"",
    "a;mkdir sfxdir13", "a|mkdir sfxdir14", "a&mkdir sfxdir15", "a>sfxdir16", "a<b", "(mkdir sfxdir17)",
    "{a}", "a\nmkdir sfxdir18", "a&&mkdir sfxdir19",
  };
  for (const std::string& s : bad) {
    std::string what;
    try {
      std::vector<std::string> e = casadi::expand(s);
      what = "accepted -> " + vesc(e);
    } catch (std::exception& ex) {
      what = std::string("refused: ") + ex.what();
    }
    bool refused = what.compare(0, 7, "refused")==0;
    bool sidefx = false;
    for (int i=1; i<=19; ++i) if (exists("sfxdir" + std::to_string(i))) sidefx = true;
    if (exists("sfxdir_w")) sidefx = true;
    if (sidefx) fail("side effect for input=" + esc(s) + " (" + what + ")");
    else pass();
    if (!refused) note("not refused: input=" + esc(s) + " " + what);
    std::cout << "   " << esc(s) << ": " << what.substr(0, 140) << std::endl;
  }
  change_dir(cwd);
}
#endif

// ---------------------------------------------------------------- win_run
#ifdef _WIN32
static void test_win_run() {
  section("win_run");
  if (g_dumps.empty()) return;
  const std::string& exe = g_dumps[0].path;
  for (int code : {0, 1, 3, 42, 255, 1000}) {
    set_env("ADUMP_EXIT", std::to_string(code));
    Dump r = dump_args(exe, {});
    check(r.ret==code, "exit code " + std::to_string(code) + " gave " + std::to_string(r.ret));
  }
  unset_env("ADUMP_EXIT");
  {
    std::vector<std::string> a = {"no_such_program_harness.exe", "x"};
    int ret = casadi::run(a, nullptr);
    check(ret==-1, "nonexistent program gave " + std::to_string(ret));
    a = {"C:\\no such dir\\x.exe"};
    ret = casadi::run(a, nullptr);
    check(ret==-1, "nonexistent path gave " + std::to_string(ret));
    a = {""};
    ret = casadi::run(a, nullptr);
    check(ret==-1, "empty program gave " + std::to_string(ret));
  }
  // Environment block
  {
    remove(g_outfile.c_str());
    casadi::WinEnv env;
    const char* sr = getenv("SystemRoot");
    std::wstring vals[] = {
      L"ADUMP_ENV=HARNESS_BLOCK",
      L"ADUMP_OUT=" + widen(g_outfile),
      L"HARNESS_BLOCK=from block \u00e9",
      L"SystemRoot=" + widen(sr ? sr : "C:\\Windows"),
    };
    for (const std::wstring& v : vals) {
      env.block += v;
      env.block += L'\0';
    }
    env.block += L'\0';
    int ret = casadi::win_run(casadi::quote(exe), &env);
    Dump r = parse_dump(read_file(g_outfile), ret);
    check(ret==0 && r.env_set && r.env=="from block \xE9", "env block: ret=" + std::to_string(ret)
      + " env=" + esc(r.env));
    // Parent variable must not leak into the child
    set_env("HARNESS_PARENT_ONLY", "1");
    env.block = L"ADUMP_ENV=HARNESS_PARENT_ONLY";
    env.block += L'\0';
    env.block += L"ADUMP_OUT=" + widen(g_outfile);
    env.block += L'\0';
    env.block += L'\0';
    remove(g_outfile.c_str());
    ret = casadi::win_run(casadi::quote(exe), &env);
    r = parse_dump(read_file(g_outfile), ret);
    check(ret==0 && !r.env_set, "env block replaces parent environment");
  }
  // Stdout capture, large output, exit code
  {
    unset_env("ADUMP_OUT");
    set_env("ADUMP_SPAM", "3000000");
    set_env("ADUMP_EXIT", "7");
    std::string out;
    int ret = casadi::win_run(casadi::quote(exe) + " cap", nullptr, &out);
    size_t nx = 0;
    for (char c : out) if (c=='x') nx++;
    check(ret==7 && nx==3000000, "capture: ret=" + std::to_string(ret) + " x count="
      + std::to_string(nx) + " size=" + std::to_string(out.size()));
    unset_env("ADUMP_SPAM");
    unset_env("ADUMP_EXIT");
    out.clear();
    ret = casadi::win_run("no_such_program_harness.exe", nullptr, &out);
    check(ret==-1 && out.empty(), "capture with nonexistent program gave " + std::to_string(ret));
  }
  // Program in a directory with spaces and punctuation
  {
    std::string d = join(g_work, "dir with space & (x) ^ %HARNESS_V% !x!");
    make_dir(d);
    std::string copy = join(d, "argv dump.exe");
    CopyFileA(exe.c_str(), copy.c_str(), FALSE);
    Dump r = dump_args(copy, {"a b"});
    check(r.ok && r.ret==0 && r.crt==std::vector<std::string>{"a b"}, "program path "
      + esc(copy) + " ret=" + std::to_string(r.ret) + " crt=" + vesc(r.crt));
  }
  // Console of the child under a console parent
  {
    Dump r = dump_args(exe, {});
    std::cout << "   parent console window: " << (GetConsoleWindow() ? "yes" : "no")
              << ", child console state: " << r.console << std::endl;
  }
}

static const char* jit_source_probe();
static bool jit_try_probe(const std::string& src, const std::string& outdir, std::string& err);

// Child mode: started without a console
static int console_probe(const std::string& out) {
  std::ostringstream ss;
  ss << "self_console " << (GetConsoleWindow() ? 1 : 0) << "\n";
  Dump r = dump_args(g_dumps[0].path, {"probe"});
  ss << "child_ret " << r.ret << "\nchild_console " << r.console << "\n";
  // Baseline: plain CreateProcess without CREATE_NO_WINDOW
  {
    remove(g_outfile.c_str());
    STARTUPINFOA si;
    ZeroMemory(&si, sizeof(si));
    si.cb = sizeof(si);
    PROCESS_INFORMATION pi;
    std::string cmd = casadi::quote(g_dumps[0].path);
    std::vector<char> buf(cmd.begin(), cmd.end());
    buf.push_back('\0');
    if (CreateProcessA(nullptr, buf.data(), nullptr, nullptr, FALSE, 0, nullptr, nullptr,
        &si, &pi)) {
      WaitForSingleObject(pi.hProcess, INFINITE);
      CloseHandle(pi.hProcess);
      CloseHandle(pi.hThread);
    }
    Dump b = parse_dump(read_file(g_outfile), 0);
    ss << "baseline_child_console " << b.console << "\n";
  }
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
  std::string cmd = casadi::quote(self) + " --work " + casadi::quote(g_work) + " --dump "
    + casadi::quote(g_dumps[0].label + "=" + g_dumps[0].path) + " --jit " + g_jit + " --console-probe "
    + casadi::quote(out);
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
  WaitForSingleObject(pi.hProcess, 60000);
  CloseHandle(pi.hProcess);
  CloseHandle(pi.hThread);
  std::string res = read_file(out);
  std::cout << res;
  check(res.find("self_console 0")!=std::string::npos, "probe has no console");
  check(res.find("child_ret 0")!=std::string::npos, "child ran from a console-less parent");
  check(res.find("\nchild_console 2")==std::string::npos, "child got a visible console window");
  if (g_jit!="none") check(res.find("probe_jit ok")!=std::string::npos, "JIT from a console-less parent");
}

static void test_batch_status() {
  section("win_env: batch file exit status");
  struct Case { std::string body; bool fails; };
  std::vector<Case> cases = {
    {"@exit /b 0\r\n", false},
    {"@exit /b 1\r\n", true},
    {"@echo off\r\nexit /b 1\r\n", true},
    {"@echo off\r\nsetlocal\r\nexit /b 1\r\n", true},
    {"@echo off\r\ngoto :end\r\n:end\r\nexit /B 1\r\n", true},
    {"@echo off\r\ncall :sub\r\nexit /b 1\r\n:sub\r\nexit /b 0\r\n", true},
    {"@echo off\r\nset __S=1\r\nendlocal & exit /b %__S%\r\n", true},
    {"@echo off\r\nno_such_command_harness\r\n", true},
    {"@echo off\r\nset HARNESS_FROM_BAT=1\r\n", false},
  };
  std::string d = join(g_work, "bat dir & (x)");
  make_dir(d);
  int k = 0;
  for (const Case& c : cases) {
    std::string bat = join(d, "case" + std::to_string(k++) + ".bat");
    write_file(bat, c.body);
    bool threw = false;
    try {
      casadi::win_env(bat, "x64");
    } catch (std::exception&) {
      threw = true;
    }
    // Reference: with call
    const char* comspec = getenv("ComSpec");
    std::string cmd = casadi::quote(comspec ? comspec : "cmd.exe") + " /d /s /c \"call "
      + casadi::quote(bat) + " x64 >nul && exit 0 || exit 1\"";
    int ref = casadi::win_run(cmd, nullptr);
    check(threw==c.fails, "batch " + esc(c.body) + ": win_env threw=" + std::to_string(threw)
      + " expected=" + std::to_string(c.fails) + "; with 'call' the status is " + std::to_string(ref));
  }
}

static void test_vs_setup() {
  section("Visual Studio setup");
  std::string full;
  bool on_path = casadi::win_find("cl.exe", nullptr, full);
  std::cout << "   cl.exe on PATH: " << (on_path ? full : "no") << std::endl;
  SetEnvironmentVariableW(L"HARNESS_UNI", L"\u00e9\u4e2d \U0001F600 x=y");
  std::string bat = Peek::vcvarsall();
  std::cout << "   vcvarsall: " << bat << std::endl;
  check(!bat.empty() && exists(bat), "vcvarsall.bat located");
  if (bat.empty()) return;
  check(Peek::vcvarsall()==bat, "vcvarsall cached");
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
    std::wstring blk(env->block.data(), env->block.size());
    auto has = [&](const wchar_t* k) {
      std::wstring key = std::wstring(L"\0", 1) + k;
      return blk.compare(0, wcslen(k), k)==0 || blk.find(key)!=std::wstring::npos;
    };
    check(has(L"INCLUDE="), a + ": INCLUDE in environment");
    check(has(L"LIB="), a + ": LIB in environment");
    check(has(L"Path=") || has(L"PATH="), a + ": PATH in environment");
    check(!env->path.empty(), a + ": PATH captured");
    check(blk.find(std::wstring(L"\0HARNESS_UNI=\u00e9\u4e2d \U0001F600 x=y\0", 23))!=std::wstring::npos,
      a + ": non-ASCII variable preserved");
    check(env->block.size()>=2 && env->block[env->block.size()-1]==L'\0'
      && env->block[env->block.size()-2]==L'\0', a + ": block double-NUL terminated");
    std::string cl, link;
    bool f1 = casadi::win_find("cl.exe", env, cl), f2 = casadi::win_find("link.exe", env, link);
    std::cout << "   " << a << ": cl=" << cl << " link=" << link << " block="
              << env->block.size() << " wchars" << std::endl;
    check(f1 && cl.find("\\MSVC\\")!=std::string::npos, a + ": cl.exe resolved: " + cl);
    check(f2 && link.find("\\MSVC\\")!=std::string::npos, a + ": link.exe resolved: " + link);
    check(casadi::win_find("cl", env, cl), a + ": 'cl' without extension resolved");
    check(casadi::win_env(bat, a)==env, a + ": cached");
    // The environment runs the tools
    std::string out;
    int ret = casadi::win_run(casadi::quote(cl), env, &out);
    check(ret!=-1, a + ": cl.exe starts with the environment, ret=" + std::to_string(ret));
  }
  try {
    const casadi::WinEnv* env = casadi::win_env(bat, "bogus_arch");
    std::wstring blk(env->block.data(), env->block.size());
    std::string cl;
    fail("win_env with an invalid architecture did not throw; INCLUDE in block: "
      + std::to_string(blk.find(std::wstring(L"\0INCLUDE=", 9))!=std::wstring::npos)
      + ", cl.exe found: " + std::to_string(casadi::win_find("cl.exe", env, cl)) + " " + cl);
  } catch (std::exception& e) {
    pass();
    std::cout << "   bogus_arch: " << std::string(e.what()).substr(0, 200) << std::endl;
  }
}
#endif

// ---------------------------------------------------------------- JIT
static const char* jit_source =
  "#ifdef _WIN32\n#define EXPORT __declspec(dllexport)\n#else\n#define EXPORT\n#endif\n"
  "#define XSTR(x) STR(x)\n#define STR(x) #x\n"
  "#ifdef USE_INC\n#include \"harness_inc.h\"\n#else\n#define INC_VALUE 0\n#endif\n"
  "EXPORT int square(int x) { return x*x + INC_VALUE; }\n"
  "#ifdef TAG\nEXPORT const char* tag(void) { return XSTR(TAG); }\n#endif\n";

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

static bool jit_try(const std::string& src, const std::string& outdir, Dict opts, int want,
    const char* want_tag, std::string& err) {
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
    if (f(7)!=want) {
      err = "square(7) gave " + std::to_string(f(7));
      return false;
    }
    if (want_tag) {
      typedef const char* (*tag_t)(void);
      tag_t t = reinterpret_cast<tag_t>(sc.get_function("tag"));
      if (!t) {
        err = "symbol 'tag' not found";
        return false;
      }
      if (std::string(t())!=want_tag) {
        err = "tag() gave " + esc(t()) + ", expected " + esc(want_tag);
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
  return jit_try(src, outdir, jit_base(), 49, nullptr, err);
}

static std::vector<std::string> sentinels(const std::string& d) {
  std::vector<std::string> r;
  for (const char* n : {"sfxdir", "sfxdir2"}) {
    for (const std::string& base : {d, g_work, current_dir()}) {
      if (exists(join(base, n))) r.push_back(join(base, n));
    }
  }
  return r;
}

static void test_jit_dirs() {
  section("JIT in directories with spaces and punctuation (" + g_jit + ")");
  set_env("HARNESS_V", "hello world");
  std::vector<std::string> names = {
    "plain", "with space", "O'Brien", "R&D (x86)", "a&mkdir sfxdir&b", "c^d", "c^&mkdir sfxdir2",
    "100%", "pct%HARNESS_V%dir", "(paren)", "semi;colon", "bang!x!", "comma,eq=", "tilde~1",
    "trail.", "\xC3\xA9t\xC3\xA9", "  two  spaces",
#ifndef _WIN32
    "x|y", "a;mkdir sfxdir;b", "dollar$HOME", "back`mkdir sfxdir`", "$(mkdir sfxdir2)",
    "quote\"dq", "star*", "brace{a,b}", "lt<gt>", "a\\b", "new\nline", "a&&mkdir sfxdir",
#endif
  };
  std::string root = join(g_work, "jitdirs");
  make_dir(root);
  std::string cwd = current_dir();
  for (const std::string& n : names) {
#ifdef _WIN32
    if (!n.empty() && n.back()=='.') continue;  // Windows strips trailing dots
#endif
    std::string d = join(root, n);
    if (!make_dir(d)) {
      note("cannot create directory " + esc(d));
      continue;
    }
    change_dir(d);
    std::string src = join(d, "jit src.c");
    write_file(src, jit_source);
    std::string err;
    bool ok = jit_try(src, d + SEP, jit_base(), 49, nullptr, err);
    check(ok, "dir=" + esc(n) + ": " + err);
    std::vector<std::string> s = sentinels(d);
    check(s.empty(), "dir=" + esc(n) + ": side effect " + vesc(s));
    change_dir(cwd);
  }
}

static void test_jit_options() {
  section("JIT options (" + g_jit + ")");
  std::string d = join(g_work, "jitopts dir");
  make_dir(d);
  std::string src = join(d, "jit src.c");
  write_file(src, jit_source);
  std::string inc = join(d, "inc dir & (x)");
  make_dir(inc);
  write_file(join(inc, "harness_inc.h"), "#define INC_VALUE 1000\n");
  std::string err;
  bool ok;
  set_env("HARNESS_V", "hello world");

  Dict o = jit_base();
  o["flags"] = std::vector<std::string>{def_flag() + "TAG=\"a b\""};
  ok = jit_try(src, d + SEP, o, 49, "a b", err);
  check(ok, "flags " + esc(def_flag() + "TAG=\"a b\"") + ": " + err);

#ifdef _WIN32
  o = jit_base();
  o["flags"] = std::vector<std::string>{"\"" + def_flag() + "TAG=%HARNESS_V%\""};
  ok = jit_try(src, d + SEP, o, 49, "hello world", err);
  check(ok, "flags with %VAR% in quotes: " + err);
#else
  o = jit_base();
  o["flags"] = std::vector<std::string>{"\"" + def_flag() + "TAG=$HARNESS_V\""};
  ok = jit_try(src, d + SEP, o, 49, "hello world", err);
  check(ok, "flags with $VAR in quotes: " + err);
#endif

  // Compiler given through a variable
  if (g_jit!="none") {
#ifdef _WIN32
    std::string cc = g_jit=="cl" ? "cl.exe" : "gcc";
    set_env("HARNESS_CC", cc);
    o = jit_base();
    o["compiler"] = "%HARNESS_CC%";
    if (g_jit=="gcc") o["linker"] = "%HARNESS_CC%";
#else
    const char* ccenv = getenv("CC");
    std::string cc = g_jit=="gcc" ? "gcc" : (ccenv ? ccenv : "cc");
    set_env("HARNESS_CC", cc);
    o = jit_base();
    o["compiler"] = "$HARNESS_CC";
    o["linker"] = "$HARNESS_CC";
#endif
    ok = jit_try(src, d + SEP, o, 49, nullptr, err);
    check(ok, "compiler through a variable: " + err);
  }

  // Include directory with a space
  std::string inc2 = join(d, "inc space");
  make_dir(inc2);
  write_file(join(inc2, "harness_inc.h"), "#define INC_VALUE 1000\n");
  o = jit_base();
  o["flags"] = std::vector<std::string>{def_flag() + "USE_INC"};
  o["include_dirs"] = std::vector<std::string>{inc2};
  ok = jit_try(src, d + SEP, o, 1049, nullptr, err);
  check(ok, "include_dirs " + esc(inc2) + ": " + err);

  // Include directory with spaces and punctuation
  o = jit_base();
  o["flags"] = std::vector<std::string>{def_flag() + "USE_INC"};
  o["include_dirs"] = std::vector<std::string>{inc};
  ok = jit_try(src, d + SEP, o, 1049, nullptr, err);
  check(ok, "include_dirs " + esc(inc) + ": " + err);

  o = jit_base();
  o["flags"] = std::vector<std::string>{def_flag() + "USE_INC", "\"" + inc_flag() + inc + "\""};
  ok = jit_try(src, d + SEP, o, 1049, nullptr, err);
  check(ok, "quoted include flag " + esc(inc) + ": " + err);

  // Output flag variants
  if (g_jit=="gcc" || g_jit=="cc") {
    o = jit_base();
    o["compiler_output_flag"] = "-o";
    o["linker_output_flag"] = "-o";
    ok = jit_try(src, d + SEP, o, 49, nullptr, err);
    check(ok, "glued output flags '-o': " + err);
  }

  // A compiler that cannot be started gives its own error
  o = jit_base();
  o["compiler"] = "no_such_compiler_harness";
  ok = jit_try(src, d + SEP, o, 49, nullptr, err);
  check(!ok && err.find("Could not start compiler")!=std::string::npos,
    "missing compiler error: " + err);
}

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
      oks[i] = jit_try(src, d + SEP, jit_base(), 49, nullptr, err);
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
  for (const DumpExe& d : g_dumps) std::cout << "dump " << d.label << "=" << d.path << std::endl;

  if (want("threads") && g_jit!="none") test_jit_threads();
  if (want("quote")) test_quote_roundtrip();
#ifdef _WIN32
  if (want("split")) test_expand_split();
  if (want("vars")) test_expand_vars();
  if (want("winrun")) test_win_run();
  if (want("console")) test_no_console();
  if (want("vs") && g_jit=="cl") test_vs_setup();
  if (want("bat")) test_batch_status();
#else
  if (want("split")) test_expand_posix();
#endif
  if (want("jit") && g_jit!="none") {
    test_jit_dirs();
    test_jit_options();
  }
  std::cout << "SUMMARY pass=" << n_pass << " fail=" << n_fail << " note=" << n_note << std::endl;
  return n_fail ? 1 : 0;
}
