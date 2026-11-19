/* Loads the plugin code with RTLD_DEEPBIND and starts argv_dump after setenv */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern char** environ;

int main(int argc, char** argv) {
  if (argc<4) return 2;
  /* argv[1]: library, argv[2]: argv_dump, argv[3]: output file */
  printf("environ in executable: %p\n", (void*)environ);
  void* h = dlopen(argv[1], RTLD_NOW | RTLD_LOCAL | RTLD_DEEPBIND);
  if (!h) { printf("dlopen: %s\n", dlerror()); return 3; }
  int (*run)(const char* const*) = (int (*)(const char* const*))dlsym(h, "harness_run");
  int (*live)(void) = (int (*)(void))dlsym(h, "harness_environ_is_live");
  setenv("HARNESS_LATE", "late value", 1);
  setenv("ADUMP_ENV", "HARNESS_LATE", 1);
  setenv("ADUMP_OUT", argv[3], 1);
  printf("plain environ in library is live: %d\n", live());
  const char* a[] = {argv[2], "x", NULL};
  int ret = run(a);
  printf("run returned %d\n", ret);
  FILE* f = fopen(argv[3], "rb");
  char buf[4096] = {0};
  if (f) { if (fread(buf, 1, sizeof(buf)-1, f)) {} fclose(f); }
  printf("%s", buf);
  /* "late value" hex */
  int ok = ret==0 && strstr(buf, "env 6c6174652076616c7565")!=NULL;
  printf("%s deepbind environment\n", ok ? "PASS" : "FAIL");
  return ok ? 0 : 1;
}
