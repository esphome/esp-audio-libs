# FLAC Decoder Test Suite

This document describes how to test the esp-audio-libs FLAC decoder using the official FLAC test suite.

## Overview

The test suite validates the FLAC decoder against the official FLAC test files. It performs validation using:
1. **MD5 verification** (primary): Decoded audio MD5 matches the signature in the FLAC header
2. **ffmpeg comparison** (secondary): Cross-validation against ffmpeg's output

## Prerequisites

1. **Python 3** with standard library
2. **ffmpeg** (for reference decoding)
3. **FLAC test files** from the official test suite

## Setting Up Test Files

Download the official FLAC test suite:

```bash
# From the flac_to_wav directory
git clone https://github.com/ietf-wg-cellar/flac-test-files.git
```

This creates a `flac-test-files` directory in `flac_to_wav`.

Expected directory structure:
```
esp-audio-libs/
├── host_examples/
│   └── flac_to_wav/
│   ├── flac-test-files/        # Test suite (cloned from GitHub)
│       │   ├── subset/             # Standard FLAC files
│       │   ├── uncommon/           # Edge cases
│       │   └── faulty/             # Invalid files
│       ├── build/
│       ├── flac_to_wav.cpp
│       ├── test_flac_decoder.py
│       ├── README.md
│       └── TESTING.md (this file)
```

## Running Tests

```bash
# Build the decoder first
cmake -B build
cmake --build build

# Run the test suite
python3 test_flac_decoder.py
```

## Understanding Test Results

The test script will:
1. Decode each FLAC file using our decoder
2. Verify the MD5 hash computed during decoding against the FLAC header signature
3. Cross-validate against ffmpeg output (for additional confidence)
4. Generate a detailed report

**Test output example:**
```
Testing subset files (64 files)...
  [1/64] Testing 01 - blocksize 4096.flac... ✓ PASS
  [2/64] Testing 02 - blocksize 4608.flac... ✓ PASS
  ...
```

**Results are saved to:**
- `test_results/test_report.txt` - Human-readable report
- `test_results/test_report.json` - Machine-readable JSON report
- `test_results/<category>/our_decoder/` - WAV files from our decoder
- `test_results/<category>/ffmpeg/` - WAV files from ffmpeg

## Test Categories

The test suite includes three categories:

### 1. subset (64 files)
Standard FLAC files using streaming subset features (should all pass)
- Various bit depths (8, 12, 16, 20, 24, 32 bits)
  - All bit depths supported with proper MD5 verification
  - Non-byte-aligned depths (12, 15, 20 bits) are correctly packed per FLAC spec for MD5 computation
- Various sample rates (22.05kHz to 384kHz)
- Various block sizes (16 to 16384 samples)
- Various channel configurations (mono to 8 channels)
- Files with large metadata blocks
- Variable block size files

### 2. uncommon (11 files)
Files with uncommon or exotic features of the FLAC format
- Uncommon bit depths (15 or 32 bits) - passes
- Uncommon sample rate (768kHz) - passes
- Uncommon block size (65535) or rice partition orders (15) - passes
- Changes to streaming settings (sample rate, number of channels, bit depth) - fails (not supported)
- Files starting with frame header or unparseable data - fails (requires STREAMINFO)

### 3. faulty (11 files)
Files with invalid data or corruption
- Expected to fail gracefully if unable to decode
- Some files may be accepted if the error is in metadata we don't validate

## Interpreting Test Results

### Primary Test Results (MD5-based)

**✓ PASS - MD5 verified**
- Computed MD5 matches header signature
- ✅ Cryptographically verified bit-perfect decode
- Most reliable test result

**✓ PASS - MD5 verified + matches ffmpeg**
- MD5 matches header signature
- Also matches ffmpeg output
- ✅ Highest confidence result (dual validation)

**✗ FAIL - MD5 mismatch**
- Computed MD5 does not match header signature
- ❌ Decoder bug, needs investigation

**✓ PASS - Matches ffmpeg (no MD5 in file)**
- No MD5 signature available in FLAC header
- Falls back to ffmpeg comparison
- ✅ Correct decode (secondary validation)

### Other Results

**? FAIL - Decoder failed but ffmpeg succeeded**
- Our decoder rejected the file
- May indicate missing feature support or decoder bug
- ⚠️ Needs investigation

**? EXPECTED - Both decoders failed**
- Both our decoder and ffmpeg rejected the file
- ✅ Correct behavior for invalid files

**? UNKNOWN - No MD5 available and ffmpeg failed**
- Cannot validate output
- ⚠️ Manual inspection needed

## Troubleshooting Tests

**"Error: flac_to_wav not found"**
```bash
# Build the decoder first
cmake -B build
cmake --build build
```

**"Error: Test files not found"**
```bash
# Clone the test files repository
cd ..
git clone https://github.com/ietf-wg-cellar/flac-test-files.git
cd flac_to_wav
```

**"Error: ffmpeg not found"**
```bash
# On macOS:
brew install ffmpeg

# On Ubuntu/Debian:
sudo apt-get install ffmpeg

# On Windows:
# Download from https://ffmpeg.org/download.html
```

**Test suite runs but all tests fail**
- Check that `flac_to_wav` is executable
- Try running manually: `./flac_to_wav "flac-test-files/subset/01 - blocksize 4096.flac" test.wav`
- Check console output for error messages

## Current Test Results

As of the latest run, the decoder achieves:
- **73/86 tests passing (84%)**
- **All 64 subset files pass** (100%) - all with MD5 verification
- Known limitations:
  - Files with changing stream parameters (sample rate, channels, bit depth) mid-stream
  - Files without STREAMINFO metadata block

## Customizing Tests

Edit `test_flac_decoder.py` to:
- Test specific file categories: Modify `TEST_CATEGORIES` list
- Change timeout for slow systems: Adjust `timeout=` parameter in `run_command()`
- Test against different reference decoder: Replace ffmpeg commands
- Add custom validation logic

## Further Reading

- [FLAC Format Specification](https://xiph.org/flac/format.html)
- [FLAC Test Files Repository](https://github.com/ietf-wg-cellar/flac-test-files)
- [RFC 9639: FLAC Specification](https://www.rfc-editor.org/rfc/rfc9639.html)
