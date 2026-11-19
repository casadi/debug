/* Prints its arguments, hex encoded, one per line.
   Control through environment variables:
     ADUMP_OUT=<file>   write to file instead of stdout
     ADUMP_ENV=<name>   also print the value of that variable
     ADUMP_SPAM=<n>     print n extra bytes to stdout
     ADUMP_EXIT=<n>     exit code */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#include <windows.h>
#include <shellapi.h>
#endif

static void hex(FILE* f, const char* tag, const char* s) {
  fprintf(f, "%s ", tag);
  if (!*s) fprintf(f, "-");
  for (; *s; ++s) fprintf(f, "%02x", (unsigned char)*s);
  fprintf(f, "\n");
}

int main(int argc, char** argv) {
  FILE* f = stdout;
  const char* out = getenv("ADUMP_OUT");
  const char* e = getenv("ADUMP_ENV");
  const char* spam = getenv("ADUMP_SPAM");
  const char* ex = getenv("ADUMP_EXIT");
  int i;
  if (out && *out) {
    f = fopen(out, "wb");
    if (!f) return 99;
  }
  fprintf(f, "crt %d\n", argc-1);
  for (i=1; i<argc; ++i) hex(f, "a", argv[i]);
#ifdef _WIN32
  {
    int n = 0;
    LPWSTR* w = CommandLineToArgvW(GetCommandLineW(), &n);
    fprintf(f, "w %d\n", n-1);
    for (i=1; i<n; ++i) {
      int len = WideCharToMultiByte(CP_ACP, 0, w[i], -1, NULL, 0, NULL, NULL);
      char* buf = (char*)malloc(len+1);
      WideCharToMultiByte(CP_ACP, 0, w[i], -1, buf, len, NULL, NULL);
      hex(f, "b", buf);
      free(buf);
    }
    LocalFree(w);
    hex(f, "raw", GetCommandLineA());
    fprintf(f, "console %d\n", GetConsoleWindow() ? (IsWindowVisible(GetConsoleWindow()) ? 2 : 1) : 0);
  }
#endif
  if (e) {
    const char* v = getenv(e);
    if (v) hex(f, "env", v); else fprintf(f, "env (unset)\n");
  }
  if (f!=stdout) fclose(f);
  if (spam) {
    long n = atol(spam);
    long k;
    for (k=0; k<n; ++k) putchar('x');
    fflush(stdout);
  }
  return ex ? atoi(ex) : 0;
}
