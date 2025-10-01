# FLAC Decoder

This directory contains a FLAC (Free Lossless Audio Codec) decoder optimized for ESP32 devices.

## Overview

The decoder is based on [Nayuki's Simple FLAC Implementation](https://www.nayuki.io/res/simple-flac-implementation/) and [Mike Hansen's initial C++ FLAC decoder port](https://github.com/synesthesiam/flac-decoder). It has been extensively enhanced and optimized for embedded systems, particularly ESP32 microcontrollers.

**Specification**: [FLAC Format Specification](https://xiph.org/flac/format.html)

## Features

- **Streaming support**: Both header and frame decoding can be done incrementally with small buffers
- **Memory efficient**: Configurable metadata size limits to control memory usage
- **Hardware optimized**: Assembly implementations for ESP32 (Xtensa AE32/AES3 instructions)
- **CRC validation**: Optional frame integrity checking with CRC-8 (frame header) and CRC-16 (frame data)
- **Full metadata support**: Extracts metadata blocks including album art (with size limits)
- **MD5 signature extraction**: Provides access to the 128-bit MD5 signature from STREAMINFO for external validation
- **Automatic overflow detection**: Selects 32-bit or 64-bit accumulator for LPC restoration based on sample depth and coefficients

## File Organization

### Public API
- `flac_decoder.h` (in `include/`) - Main decoder class interface

### Core Decoder
- `flac_decoder.cpp` - Main decoder implementation
  - Header and metadata parsing
  - Frame decoding orchestration
  - Subframe decoding (constant, verbatim, fixed, LPC)
  - Bit stream reading utilities

### Linear Predictive Coding (LPC)
- `flac_lpc.h` / `flac_lpc.cpp` - LPC restoration helpers
  - Overflow detection logic
  - 32-bit arithmetic implementation (fast path)
  - 64-bit arithmetic implementation (safe for high-resolution audio)

### Assembly Optimizations (ESP32 only)
- `flac_lpc_32_asm.S` - 32-bit LPC restoration (Xtensa assembly)
- `flac_lpc_64_asm.S` - 64-bit LPC restoration (Xtensa assembly)
- `flac_lpc_asm.h` - Assembly function declarations
- `flac_lpc_platform.h` - Platform detection for assembly support

### Utilities
- `flac_crc.h` / `flac_crc.cpp` - CRC-8 and CRC-16 lookup tables and functions

## Architecture

The FLACDecoder class is designed as a self-contained, non-inheritable implementation. All methods are either public API or private implementation details


### Optimizations

#### Bitstream Reading

The decoder uses a 32-bit integer bit buffer to avoid unnecessary 64-bit arithmetic. It always attempts to refill the bit buffer with 4 bytes whenever needed. Read functions are inlined.

#### LPC Accumulator Type

The decoder automatically determines whether 32-bit or 64-bit arithmetic is safe:

1. **Overflow Detection** (`can_use_32bit_lpc()`):
   - Analyzes sample depth, LPC coefficients, order, and quantization shift
   - Calculates maximum possible values before and after shift
   - Determines if 32-bit arithmetic will overflow

2. **32-bit Fast Path** (`restore_linear_prediction_32bit()`):
   - Used when overflow is impossible (most 16-bit audio)
   - ~2x faster on ESP32-S3
   - Assembly implementation with loop unrolling for orders 1-12

3. **64-bit Safe Path** (`restore_linear_prediction_64bit()`):
   - Used for high-resolution audio (e.g., 24-bit) with large coefficients
   - Prevents overflow in all valid FLAC streams
   - Assembly implementation uses MULL/MULSH instructions on ESP32 with loop unrolling for orders 1-12

#### Output Packing

The decoder uses optimized packing algorithms for output audio from the internal 32-bit representation. All paths use unrolled loops, with specialized fast paths for the most common audio scenarios: mono 16-bits per sample, stereo 16-bits per sample, and stereo 24-bits per sample.

### Decoder Pipeline

The FLAC decoder processes audio data through the following stages:

```
┌─────────────────┐
│  Input Buffer   │
└────────┬────────┘
         │
         ├──────────────────────────┐
         │                          │
         ▼                          ▼
┌─────────────────┐        ┌──────────────┐
│ Header Parser   │        │  Frame Sync  │
└────────┬────────┘        └──────┬───────┘
         │                        │
         ▼                        ▼
┌─────────────────┐        ┌──────────────┐
│ Metadata Blocks │        │ Frame Header │
│  - STREAMINFO   │        │   (CRC-8)    │
│  - VORBIS_      │        └──────┬───────┘
│    COMMENT      │               │
│  - PICTURE      │               ▼
│  - etc.         │        ┌──────────────┐
└─────────────────┘        │  Subframes   │
                           │  - Constant  │
                           │  - Verbatim  │
                           │  - Fixed     │
                           │  - LPC       │
                           └──────┬───────┘
                                  │
                                  ▼
                           ┌──────────────┐
                           │   Residuals  │
                           │ (Rice coded) │
                           └──────┬───────┘
                                  │
                                  ▼
                           ┌──────────────┐
                           │   Channel    │
                           │ Decorrelation│
                           └──────┬───────┘
                                  │
                                  ▼
                           ┌──────────────┐
                           │ Frame CRC-16 │
                           │  Validation  │
                           └──────┬───────┘
                                  │
                                  ▼
                           ┌──────────────┐
                           │Output Buffer │
                           │(Interleaved  │
                           │ PCM samples) │
                           └──────────────┘
```

### Streaming Support

The decoder is designed for streaming scenarios where the entire file isn't in memory:

- **Header Parsing**: Returns `FLAC_DECODER_HEADER_OUT_OF_DATA` when more data is needed
- **Partial Metadata**: Accumulates large metadata blocks (e.g., album art) across multiple calls
- **Frame Sync**: Searches for frame boundaries in the input buffer
- **Buffer Tracking**: `get_bytes_index()` tells caller how many bytes were consumed

### Memory Management

- **PSRAM-aware**: Prefers PSRAM for large allocations, falls back to internal RAM
- **Customizable allocators**: FLAC_MALLOC/FLAC_FREE macros defined in header, can be overridden at compile time
- **Block samples buffer**: Allocated once based on `max_block_size` from STREAMINFO
- **Configurable metadata limits**: Per-type size limits prevent memory exhaustion (see Configuration section)
- **No dynamic reallocation**: All buffers allocated upfront after reading STREAMINFO

### Metadata Configuration

The decoder provides fine-grained control over metadata storage through configurable size limits:

**Default Limits** (optimized for ESP32):
- STREAMINFO: Always stored (required)
- PADDING: 0 bytes (skipped)
- APPLICATION: 0 bytes (skipped)
- SEEKTABLE: 0 bytes (skipped, no seeking support yet)
- VORBIS_COMMENT: 2KB (typical artist/album/title tags)
- CUESHEET: 0 bytes (skipped)
- PICTURE: 0 bytes (album art disabled by default)
- Unknown types: 0 bytes (skipped)

**API Methods**:
- `set_max_metadata_size(FLACMetadataType type, uint32_t max_size)`: Set limit for specific type
- `get_max_metadata_size(FLACMetadataType type)`: Get current limit
- `set_max_album_art_size(uint32_t max_size)`: Convenience method for PICTURE blocks
- `get_max_album_art_size()`: Get album art limit

**Usage**: Configure limits before calling `read_header()` to control memory usage based on your application's needs.

## Usage Example

```cpp
#include "flac_decoder.h"
#include <cstring>  // for std::memmove

using namespace esp_audio_libs::flac;

FLACDecoder decoder;

// Optional: Configure metadata limits before reading header
decoder.set_max_metadata_size(FLAC_METADATA_TYPE_PICTURE, 50 * 1024);  // 50KB album art
decoder.set_max_metadata_size(FLAC_METADATA_TYPE_VORBIS_COMMENT, 4 * 1024);  // 4KB tags
// Or use convenience method:
decoder.set_max_album_art_size(50 * 1024);

// Streaming buffer for reading file in chunks
// 512KB is large enough for even high-resolution multi-channel FLAC frames
const size_t streaming_buffer_size = 512 * 1024;
std::vector<uint8_t> streaming_buffer(streaming_buffer_size);
size_t bytes_in_buffer = 0;

// Read initial chunk from file
size_t bytes_read = read_from_file(streaming_buffer.data(), streaming_buffer_size);
bytes_in_buffer = bytes_read;

// Read header with streaming (handles files with large metadata)
size_t header_bytes_consumed = 0;
while (true) {
    FLACDecoderResult result = decoder.read_header(streaming_buffer.data(), bytes_in_buffer);
    header_bytes_consumed = decoder.get_bytes_index();

    if (result == FLAC_DECODER_SUCCESS) {
        // Header complete
        break;
    } else if (result == FLAC_DECODER_HEADER_OUT_OF_DATA) {
        // Need more data - shift remaining data to start of buffer
        size_t remaining = bytes_in_buffer - header_bytes_consumed;
        if (remaining > 0) {
            std::memmove(streaming_buffer.data(),
                        streaming_buffer.data() + header_bytes_consumed,
                        remaining);
        }
        bytes_in_buffer = remaining;

        // Read more data from file
        size_t new_bytes = read_from_file(streaming_buffer.data() + bytes_in_buffer,
                                         streaming_buffer_size - bytes_in_buffer);
        if (new_bytes == 0) {
            // No more data but header incomplete
            return FLAC_DECODER_ERROR_UNEXPECTED_EOF;
        }
        bytes_in_buffer += new_bytes;
        header_bytes_consumed = 0;  // Reset since we moved data
    } else {
        // Handle error
        return result;
    }
}

// Move remaining data after header to beginning of buffer
size_t remaining_after_header = bytes_in_buffer - header_bytes_consumed;
if (remaining_after_header > 0 && header_bytes_consumed > 0) {
    std::memmove(streaming_buffer.data(),
                streaming_buffer.data() + header_bytes_consumed,
                remaining_after_header);
}
bytes_in_buffer = remaining_after_header;

// Get stream info
printf("Sample rate: %u Hz\n", decoder.get_sample_rate());
printf("Channels: %u\n", decoder.get_num_channels());
printf("Bit depth: %u bits\n", decoder.get_sample_depth());

// Get MD5 signature for validation (optional)
const uint8_t* md5 = decoder.get_md5_signature();
printf("MD5 signature: ");
for (int i = 0; i < 16; i++) {
    printf("%02x", md5[i]);
}
printf("\n");

// Allocate output buffer
size_t output_size = decoder.get_output_buffer_size_bytes();
std::vector<uint8_t> output_buffer(output_size);

// Refill buffer if needed before starting frame decode
if (bytes_in_buffer < streaming_buffer_size / 2) {
    size_t new_bytes = read_from_file(streaming_buffer.data() + bytes_in_buffer,
                                     streaming_buffer_size - bytes_in_buffer);
    bytes_in_buffer += new_bytes;
}

// Decode frames with streaming
while (bytes_in_buffer > 0) {
    uint32_t samples_decoded;
    FLACDecoderResult result = decoder.decode_frame(
        streaming_buffer.data(), bytes_in_buffer,
        output_buffer.data(), &samples_decoded
    );

    if (result == FLAC_DECODER_SUCCESS) {
        // Process decoded audio
        size_t bytes_to_write = samples_decoded * decoder.get_num_channels() *
                                decoder.get_output_bytes_per_sample();
        write_audio_output(output_buffer.data(), bytes_to_write);

        // Move remaining data to start of buffer (safe after CRC validation)
        size_t bytes_consumed = decoder.get_bytes_index();
        bytes_in_buffer -= bytes_consumed;

        if (bytes_in_buffer > 0) {
            std::memmove(streaming_buffer.data(),
                        streaming_buffer.data() + bytes_consumed,
                        bytes_in_buffer);
        }

        // Refill buffer if running low
        if (bytes_in_buffer < streaming_buffer_size / 2) {
            size_t new_bytes = read_from_file(streaming_buffer.data() + bytes_in_buffer,
                                             streaming_buffer_size - bytes_in_buffer);
            bytes_in_buffer += new_bytes;
        }

    } else if (result == FLAC_DECODER_NO_MORE_FRAMES) {
        break;  // End of stream
    } else if (result == FLAC_DECODER_ERROR_OUT_OF_DATA) {
        // Frame incomplete - try to read more data
        size_t new_bytes = read_from_file(streaming_buffer.data() + bytes_in_buffer,
                                         streaming_buffer_size - bytes_in_buffer);
        if (new_bytes == 0) {
            // No more data available - unexpected end of file
            return FLAC_DECODER_ERROR_UNEXPECTED_EOF;
        }
        bytes_in_buffer += new_bytes;
        continue;  // Try decoding again
    } else {
        // Handle other errors
        return result;
    }
}
```

## Performance Characteristics

### ESP32-S3 @ 240MHz
- **16-bit/48kHz stereo**: ~3-4% CPU Usage (depends on how frequently the 64-bit LPC restoration path is necessary)
- **24-bit/48kHz stereo**: ~5% CPU Usage (always uses the 64-bit LPC restoration path)

### Memory Usage
- **Decoder object**: ~200 bytes
- **Block samples buffer**: `max_block_size * channels * 4 bytes` (typically 16-64KB)
- **Metadata**: Variable based on configured limits

## Build Configuration

The decoder automatically adapts to the build environment:

- **ESP-IDF builds**: Includes assembly optimizations, PSRAM support
- **Standalone builds**: Uses ANSI C implementations only
- **Optimization**: Built with `-O3` for maximum performance

## Testing

The decoder is tested against the [FLAC decoder test bench](https://github.com/ietf-wg-cellar/flac-test-files) suite of files. It matches ffmpeg's output for all test files that follow the [FLAC specification's streaming subset](https://www.rfc-editor.org/rfc/rfc9639.html#streamable-subset).

The test suite (`host_examples/flac_to_wav/test_flac_decoder.py`) validates:
- Bit-perfect PCM output compared to ffmpeg
- MD5 signature verification against the header value (expected to mismatch for 12-bit and 20-bit files due to byte-alignment padding)
- Proper handling of various FLAC features

Passing Test Scenarios:
- Various bit depths (8, 12, 16, 20, 24, 32 bits)
  - Output is always byte-aligned with zero-padded LSBs
- Various sample rates
- Various number of channels (1 to 8 channels)
- Large files with embedded album art
- Streaming scenarios with small buffers

## Known Limitations

- No seeking support (would require SEEKTABLE parsing)
- MD5 signature is extracted but not automatically validated during decoding (can be accessed via `get_md5_signature()` for external validation)

## References

- [FLAC Format Specification](https://xiph.org/flac/format.html)
- [Nayuki's Simple FLAC Implementation](https://www.nayuki.io/res/simple-flac-implementation/)
- [Mike Hansen's C++ FLAC decoder port](https://github.com/synesthesiam/flac-decoder)
- [FLAC decoder test bench](https://github.com/ietf-wg-cellar/flac-test-files)
- [Xtensa ISA Reference](https://www.cadence.com/content/dam/cadence-www/global/en_US/documents/tools/ip/tensilica-ip/isa-summary.pdf)