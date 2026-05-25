# Julia 1.12.6 Windows cross-build bootstrap diagnostic

A from-source (USE_BINARYBUILDER=0) cross-build of Julia 1.12.6 in the MXE
dockcross image `ghcr.io/jgillis/windows-shared-x64-posix:production` gets all
the way through `make julia-src-release` — 18 third-party deps + LLVM (static)
+ libjulia.dll + libjulia-codegen.dll + julia.exe all link cleanly. The
sysimage bootstrap step then has to run `julia.exe --output-o basecompiler-o.a`
on the build host. Under wine inside the dockcross image, that invocation
exits silently with rc=1 and zero diagnostic output — even with WINEPATH and
the MXE runtime DLLs (libgcc_s_seh-1.dll, libstdc++-6.dll, etc.) staged into
`usr/bin/`.

This branch's workflow ships the built artifacts to a real Windows runner so
we can see what julia.exe actually does on native Windows.

## How to run

`workflow_dispatch` on the `julia-win-bootstrap` branch:
1. **Linux job** (`build-linux`): pulls
   `ghcr.io/jgillis/windows-shared-x64-posix:production`, applies the
   Dockerfile here, downloads the patched Julia source tarball from this
   repo's release, runs `make -j 8 julia-src-release` inside the container,
   tars `usr/` into a workflow artifact.
2. **Windows job** (`run-windows`): downloads the artifact onto a
   `windows-2022` runner, runs `julia.exe --version` and the sysimage
   bootstrap command (`--output-o basecompiler-o.a Base_compiler.jl`),
   captures stdout/stderr/exit code.

## The 15 fixes already in the source tarball

(Full skill writeup elsewhere; the patched tarball already has all of them.)

Dockerfile shims:
1. Canonical `x86_64-w64-mingw32-*` symlinks for MXE's `.shared.posix-*`
2. Unprefixed exec wrappers (not symlinks) for `windres`/`dlltool`
3. Mixed-case Windows header symlinks (NTSecAPI.h etc.)
4. Stub `afunix.h` (MXE mingw-w64 11.2 predates it; LLVM 18 needs it)
5. Upstream CMake 3.29 (MXE's 3.20.1 has a `try_compile(COPY_FILE)` bug)
6. Custom cmake wrapper: unsets `CMAKE_TOOLCHAIN_FILE`, detects NATIVE
   sub-build, rewrites paths so cmake doesn't see CMAKE_RC_COMPILER
   "changed" between reconfigures and wipe cache

Make.user:
- `XC_HOST=x86_64-w64-mingw32`, `USE_BINARYBUILDER=0`, `OS=WINNT`
- `USE_SYSTEM_P7ZIP=1` (p7zip's Makefile assumes `-ldl`)
- `USE_LLVM_SHLIB=0` (libLLVM-18jl.dll exceeds Windows PE 65535-ordinal
  cap with MXE's bfd ld)
- `MMTK_PLAN=None` (empty default fails `ifneq (,None)` causing bare `-I`
  that swallows `JL_SYSTEM_IMAGE_PATH`)

Source patches:
- `deps/llvm.mk`: unconditional `BUILD_SHARED_LIBS=OFF`; extended
  `CROSS_TOOLCHAIN_FLAGS_NATIVE` with `CMAKE_SYSTEM_NAME=Linux`/`CROSSCOMPILING=FALSE`;
  gate `LLVM_INSTALL` DLL copy on `USE_LLVM_SHLIB`
- `deps/libgit2.mk`: replace Debian-only `CMAKE_FIND_ROOT_PATH=/usr/$(XC_HOST)`
  with MXE sysroot + `$(build_prefix)`, `MODE_INCLUDE=BOTH`
- `src/Makefile`: fix upstream paren typo in `CG_LLVMLINK` `--system-libs`
- `Makefile`: drop `julia-libllvmcalltest` from default targets

Container env:
- `unset CMAKE_TOOLCHAIN_FILE` (dockcross env-var override silently
  wins over -D flags)
- `PKG_CONFIG_PATH=PKG_CONFIG_LIBDIR=/work/usr/lib/pkgconfig` for libgit2
- `WINEPREFIX=/tmp/wineprefix` (root-owned `/tmp` defeats default wine)
- `WINEPATH` pointing at MXE runtime DLL dir
- MXE runtime DLLs copied to `/work/usr/bin/` for julia.exe dlopen
