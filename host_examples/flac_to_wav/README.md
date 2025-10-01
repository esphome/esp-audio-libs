# FLAC to WAV Converter Example

This example demonstrates how to use the esp-audio-libs FLAC decoder to convert FLAC files to WAV format.

## Overview

The `flac_to_wav` program:
- Reads a FLAC file
- Decodes it using the esp-audio-libs FLAC decoder
- Writes the decoded audio as a WAV file
- Supports all bit depths (8, 12, 16, 20, 24, 32 bits)
- Handles various channel configurations (mono to 8 channels)
- Verifies decoded audio against the MD5 signature in the FLAC header
- Displays file information and verification results

## Building

### Prerequisites

- CMake 3.10 or later
- A C++11 compatible compiler (gcc, clang, etc.)
- Make or Ninja build system

### Build Steps

```bash
# From the flac_to_wav directory
cmake -B build
cmake --build build
```

The compiled binary will be placed in the project directory as `flac_to_wav`.

### Build Options

You can specify a different build type:

```bash
# Debug build with symbols
cmake -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build

# Release build with optimizations (default)
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

## Usage

```bash
./flac_to_wav <input.flac> <output.wav>
```

### Example

```bash
# Convert a FLAC file to WAV
./flac_to_wav song.flac song.wav

# The program will display information about the file:
# - Sample rate
# - Number of channels
# - Bit depth
# - Total samples
# - Metadata blocks found
#
# After conversion, it will verify the decoded audio:
# - Expected MD5 (from FLAC header)
# - Computed MD5 (from decoded audio)
# - Verification result (PASS/FAIL)
```

### MD5 Verification

The program automatically verifies the decoded audio by computing an MD5 hash and comparing it with the signature stored in the FLAC file's STREAMINFO header. This ensures bit-perfect decoding.

**Supported for all bit depths:**
- Byte-aligned: 8, 16, 24, 32 bits
- Non-byte-aligned: 12, 15, 20 bits (properly packed per FLAC spec)

Example output:
```
=== MD5 Verification ===
Expected MD5: ac3c581ce17991866b0dcdea3b9dfd43
Computed MD5: ac3c581ce17991866b0dcdea3b9dfd43
Result: PASS - MD5 signatures match!
```

## Testing

A comprehensive test suite is available to validate the FLAC decoder against the official FLAC test files. See [TESTING.md](TESTING.md) for details on:
- Setting up test files
- Running the test suite
- Interpreting test results
- Current decoder capabilities

## Performance

On a typical modern desktop CPU:
- **16-bit stereo @ 44.1kHz**: ~100-200x realtime (decoding a 3-minute song takes <1 second)
- **24-bit stereo @ 96kHz**: ~80-150x realtime

Performance on ESP32-S3 @ 240MHz:
- **16-bit stereo @ 48kHz**: ~30x realtime (~3-4% CPU)
- **24-bit stereo @ 48kHz**: ~20x realtime (~5% CPU)

## Output Format

The converter produces standard WAV files:
- **8-bit, 16-bit**: PCM format (format code 1)
- **12-bit, 20-bit, 24-bit, 32-bit**: WAVE_FORMAT_EXTENSIBLE (format code 0xFFFE)
- **Byte order**: Little-endian (standard for WAV)
- **Bit alignment**: Non-byte-aligned bit depths are zero-padded in the LSB

Example: 12-bit samples are stored in 16-bit containers (2 bytes):
```
Original 12-bit value: 0x123
Stored in WAV:         0x1230  (left-shifted by 4, LSBs = 0)
```

## Further Reading

- [FLAC Format Specification](https://xiph.org/flac/format.html)
- [FLAC Decoder Implementation Details](../../src/decode/flac/README.md)
- [FLAC Test Files Repository](https://github.com/ietf-wg-cellar/flac-test-files)
- [RFC 9639: FLAC Specification](https://www.rfc-editor.org/rfc/rfc9639.html)