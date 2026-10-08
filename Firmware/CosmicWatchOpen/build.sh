#!/usr/bin/env bash
# Build the firmware with the tools installed in venv_cosmicwatch.nosync:
#   bin/cmake, bin/ninja                       (pip install cmake ninja)
#   opt/arm-gnu-toolchain-*-arm-none-eabi      (Arm GNU Toolchain)
#   opt/pico-sdk                               (Raspberry Pi pico-sdk + tinyusb)
# Output: build/cosmicwatch_open.uf2
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
VENV="${VENV:-$HERE/../../../venv_cosmicwatch.nosync}"
VENV="$(cd "$VENV" && pwd)"

export PICO_SDK_PATH="${PICO_SDK_PATH:-$VENV/opt/pico-sdk}"
TOOLCHAIN="$(ls -d "$VENV"/opt/arm-gnu-toolchain-*-arm-none-eabi 2>/dev/null | head -1)"
export PICO_TOOLCHAIN_PATH="${PICO_TOOLCHAIN_PATH:-$TOOLCHAIN}"
export PATH="$VENV/bin:$PICO_TOOLCHAIN_PATH/bin:$PATH"

cmake -S "$HERE" -B "$HERE/build" -G Ninja -DCMAKE_BUILD_TYPE=Release

# The SDK builds two host tools (pioasm, picotool) as external projects during
# the build. On some macOS installations the Command Line Tools ship an
# incomplete include/c++/v1 that hides the SDK's libc++ headers: point the host
# compiler to the SDK headers. CXXFLAGS is only set for the build step, so the
# cached flags of the ARM cross-compilation are not affected.
HOST_CXXFLAGS=""
if [[ "$(uname)" == Darwin ]] && ! echo '#include <string>' | c++ -x c++ -fsyntax-only - 2>/dev/null; then
    HOST_CXXFLAGS="-isystem $(xcrun --show-sdk-path)/usr/include/c++/v1"
fi
CXXFLAGS="$HOST_CXXFLAGS" cmake --build "$HERE/build"
ls -l "$HERE/build/cosmicwatch_open.uf2"
