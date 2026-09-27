#!/usr/bin/env bash
# Fetches and builds Cinder for macOS so this project can link against it.
#
#   scripts/setup_cinder.sh [Release|Debug|all]     (default: Release)
#
# Cinder is cloned next to this repo (../Cinder) unless CINDER_PATH is set.
# Requires the Xcode command line tools and CMake (brew install cmake).
set -euo pipefail

CINDER_VERSION="v0.9.3"
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CINDER_PATH="${CINDER_PATH:-$(dirname "$REPO_DIR")/Cinder}"
CONFIGS="${1:-Release}"
[[ "$CONFIGS" == "all" ]] && CONFIGS="Release Debug"

command -v cmake >/dev/null || { echo "cmake not found; install it with: brew install cmake" >&2; exit 1; }

if [[ ! -d "$CINDER_PATH/.git" ]]; then
	echo "Cloning Cinder $CINDER_VERSION into $CINDER_PATH"
	git clone --branch "$CINDER_VERSION" --depth 1 --recursive https://github.com/cinder/Cinder.git "$CINDER_PATH"
fi

# Cinder 0.9.3's bundled FreeType zlib only defines 'Byte' when TARGET_OS_MAC is undefined,
# expecting MacTypes.h to supply it; recent macOS SDKs no longer do. This is the same fix
# Cinder master uses.
ZCONF="$CINDER_PATH/src/freetype/gzip/ftzconf.h"
if grep -q '^#if !defined(MACOS) && !defined(TARGET_OS_MAC)$' "$ZCONF"; then
	echo "Patching $ZCONF"
	sed -i '' 's/^#if !defined(MACOS) \&\& !defined(TARGET_OS_MAC)$/#if !defined(Byte)/' "$ZCONF"
fi

for config in $CONFIGS; do
	echo "Building Cinder ($config)"
	cmake -S "$CINDER_PATH" -B "$CINDER_PATH/build/$config" -DCMAKE_BUILD_TYPE="$config"
	cmake --build "$CINDER_PATH/build/$config" -j "$(sysctl -n hw.ncpu)"
done

echo "Cinder is ready at $CINDER_PATH"
