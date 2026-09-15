#define _GNU_SOURCE
#include <mosek.h>
#include <dlfcn.h>
#include <stdio.h>
MSKrescodee MSK_optimizetrm(MSKtask_t task, MSKrescodee *trm) {
  static int count;
  char name[80];
  sprintf(name,"model-%d.ptf",count++);
  MSK_writedata(task,name);
  MSKrescodee (*real)(MSKtask_t,MSKrescodee*) = dlsym(RTLD_NEXT,"MSK_optimizetrm");
  return real(task,trm);
}
