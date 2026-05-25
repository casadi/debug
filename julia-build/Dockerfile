FROM ghcr.io/jgillis/windows-shared-x64-posix:production

# Canonical x86_64-w64-mingw32-* symlinks for MXE's .shared.posix tools,
# so Julia (which uses XC_HOST=x86_64-w64-mingw32, the triple autoconf
# accepts) finds its CROSS_COMPILE-prefixed tools.
RUN set -eux; \
    cd /usr/src/mxe/usr/bin; \
    for f in x86_64-w64-mingw32.shared.posix-*; do \
        target="${f#x86_64-w64-mingw32.shared.posix-}"; \
        ln -sf "$f" "x86_64-w64-mingw32-$target"; \
    done

# Unprefixed wrappers for build systems (openssl Makefile) that hardcode
# bare tool names like `windres`. These must be exec wrappers (not
# symlinks) — when MXE tools like windres look up their default include
# path, they derive it from argv[0]'s real directory. A bare symlink in
# /usr/local/bin loses that lookup; an exec to the real path keeps it.
RUN set -eux; \
    for t in windres dlltool; do \
        printf '#!/bin/sh\nexec /usr/src/mxe/usr/bin/x86_64-w64-mingw32.shared.posix-%s "$@"\n' "$t" \
            > "/usr/local/bin/$t"; \
        chmod +x "/usr/local/bin/$t"; \
    done

# Case-insensitive header shims for deps that include Windows headers with
# mixed case (e.g. p7zip uses #include <NTSecAPI.h>) but MXE's mingw-w64
# distribution ships them lowercase.
RUN set -eux; \
    inc=/usr/src/mxe/usr/x86_64-w64-mingw32.shared.posix/include; \
    for hdr in NTSecAPI CommCtrl HtmlHelp MAPI OleCtl Psapi ShlGuid ShlObj Shlwapi ShObjIdl TlHelp32 WinUser; do \
        target_lc="$(echo "$hdr" | tr A-Z a-z).h"; \
        if [ -e "$inc/$target_lc" ] && [ ! -e "$inc/${hdr}.h" ]; then \
            ln -sf "$target_lc" "$inc/${hdr}.h"; \
        fi; \
    done

# MXE's mingw-w64 11.2 predates Unix-domain-socket support on Windows
# (afunix.h, added to mingw-w64 v9+). LLVM 18's raw_socket_stream.cpp
# unconditionally includes it. Ship a minimal compatible header.
RUN cat > /usr/src/mxe/usr/x86_64-w64-mingw32.shared.posix/include/afunix.h <<'EOF'
#ifndef _AFUNIX_H_
#define _AFUNIX_H_
#define UNIX_PATH_MAX 108
typedef struct sockaddr_un {
    ADDRESS_FAMILY sun_family;
    char sun_path[UNIX_PATH_MAX];
} SOCKADDR_UN, *PSOCKADDR_UN;
#define SIO_AF_UNIX_GETPEERPID _WSAIOR(IOC_VENDOR, 256)
#endif
EOF

# Fresh upstream cmake (MXE bundles 3.20.1 which has a try_compile(COPY_FILE)
# regression that breaks check_type_size on Windows).
ARG CMAKE_VERSION=3.29.6
RUN set -eux; \
    cd /tmp; \
    curl -fsSLo cmake.tgz "https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-linux-x86_64.tar.gz"; \
    mkdir -p /opt/cmake; \
    tar -xz --strip-components=1 -C /opt/cmake -f cmake.tgz; \
    rm cmake.tgz; \
    /opt/cmake/bin/cmake --version

# Wrapper used by Julia for all cmake invocations of deps. It:
#  - Bypasses cross-compile setup on `--build / --install / -E /
#    --system-information` so build-time helpers (cmake -E copy etc.) work.
#  - Normalizes /usr/src/mxe/usr/bin/x86_64-w64-mingw32-* refs to the
#    .shared.posix-prefixed real path, so cmake doesn't see CMAKE_RC_COMPILER
#    "changed" between reconfigures and wipe the cache (was breaking
#    BLAS_LIBRARIES for SuiteSparse).
#  - Does NOT use the MXE toolchain file: it sets FIND_ROOT_PATH to MXE
#    only, hiding deps Julia installs into /work/usr (LLVM couldn't find
#    zlib). Instead we set just the essentials cross-compile vars and
#    load MXE's TryRunResults cache so try_run probes don't hang.
RUN cat > /usr/local/bin/x86_64-w64-mingw32-cmake <<'EOF'
#!/usr/bin/env bash
# The dockcross image exports CMAKE_TOOLCHAIN_FILE=mxe-conf.cmake via env
# which forces the cross compiler and CMAKE_SYSTEM_NAME=Windows on every
# cmake invocation regardless of -D args. We supply toolchain settings
# ourselves; unset the env var so it doesn't fight us.
unset CMAKE_TOOLCHAIN_FILE
# Pass-through mode for build helpers.
if echo -- "$@" | grep -Ewq -- "(--build|--install|-E|--system-information)"; then
    exec /opt/cmake/bin/cmake "$@"
fi
# When called for LLVM's NATIVE sub-build (PWD ends in /NATIVE), don't
# inject cross-compile flags — that sub-build is a host build.
# Detect LLVM NATIVE sub-build either by PWD or by an explicit
# CMAKE_SYSTEM_NAME=Linux / host-compiler arg from CROSS_TOOLCHAIN_FLAGS_NATIVE.
# In that case strip out any cross-compiler args and force host gcc.
is_native=0
case "$PWD" in */NATIVE|*/NATIVE/*) is_native=1 ;; esac
for a in "$@"; do
    case "$a" in
        -DCMAKE_SYSTEM_NAME=Linux|-DCMAKE_C_COMPILER=/usr/bin/gcc|-DCMAKE_C_COMPILER=gcc)
            is_native=1 ;;
    esac
done
if [ "$is_native" = 1 ]; then
    # Strip cross-compiler -D overrides from args so host gcc wins.
    args=()
    for a in "$@"; do
        case "$a" in
            -DCMAKE_C_COMPILER=*|-DCMAKE_CXX_COMPILER=*|-DCMAKE_RC_COMPILER=*|\
            -DCMAKE_C_COMPILER_ARG1=*|-DCMAKE_CXX_COMPILER_ARG1=*|\
            -DCMAKE_LINKER=*|-DCMAKE_AR=*|-DCMAKE_RANLIB=*|\
            -DCMAKE_SYSTEM_NAME=*|-DCMAKE_TOOLCHAIN_FILE=*)
                continue ;;
        esac
        args+=("$a")
    done
    exec /opt/cmake/bin/cmake \
        -DCMAKE_C_COMPILER=/usr/bin/gcc \
        -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
        -DCMAKE_SYSTEM_NAME=Linux \
        -DCMAKE_CROSSCOMPILING=FALSE \
        "${args[@]}"
fi
# Cross-compile mode: rewrite x86_64-w64-mingw32-* paths to their MXE
# .shared.posix-prefixed real names so cmake doesn't see CMAKE_RC_COMPILER
# "changed" between reconfigures (would wipe the cache and break things
# like SuiteSparse's BLAS detection).
args=()
for a in "$@"; do
    args+=("${a//\/usr\/src\/mxe\/usr\/bin\/x86_64-w64-mingw32-/\/usr\/src\/mxe\/usr\/bin\/x86_64-w64-mingw32.shared.posix-}")
done
exec /opt/cmake/bin/cmake \
    -DCMAKE_SYSTEM_NAME=Windows \
    -DCMAKE_SYSTEM_PROCESSOR=x86_64 \
    -DCMAKE_CROSSCOMPILING=TRUE \
    -DCMAKE_CROSSCOMPILING_EMULATOR=/usr/bin/wine \
    "-DCMAKE_FIND_ROOT_PATH=/usr/src/mxe/usr/x86_64-w64-mingw32.shared.posix;/work/usr" \
    -DCMAKE_FIND_ROOT_PATH_MODE_PROGRAM=NEVER \
    -DCMAKE_FIND_ROOT_PATH_MODE_INCLUDE=BOTH \
    -DCMAKE_FIND_ROOT_PATH_MODE_LIBRARY=BOTH \
    -DCMAKE_FIND_ROOT_PATH_MODE_PACKAGE=BOTH \
    -C /usr/src/mxe/usr/share/cmake/modules/TryRunResults.cmake \
    "${args[@]}"
EOF
RUN chmod +x /usr/local/bin/x86_64-w64-mingw32-cmake
