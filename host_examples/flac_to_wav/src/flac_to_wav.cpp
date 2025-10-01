#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <cstdint>
#include <algorithm>
#include <cstddef>
#include "flac_decoder.h"
#include "md5.h"

// Pack samples for MD5 computation according to FLAC spec
// Samples should be sign-extended to the next whole byte in little-endian order
void pack_samples_for_md5(const uint8_t* padded_samples, uint8_t* packed_output,
                          uint32_t num_samples, uint32_t bits_per_sample) {
    uint32_t bytes_per_padded_sample = (bits_per_sample + 7) / 8;
    uint32_t shift_amount = (bits_per_sample % 8 != 0) ? (8 - (bits_per_sample % 8)) : 0;

    // For byte-aligned samples, no repacking needed
    if (shift_amount == 0) {
        std::memcpy(packed_output, padded_samples, num_samples * bytes_per_padded_sample);
        return;
    }

    // For non-byte-aligned samples, unshift and sign-extend
    for (uint32_t i = 0; i < num_samples; i++) {
        const uint8_t* sample_ptr = padded_samples + (i * bytes_per_padded_sample);
        uint8_t* output_ptr = packed_output + (i * bytes_per_padded_sample);

        // Read the padded sample (LSB-padded by decoder)
        int32_t sample = 0;
        for (uint32_t byte = 0; byte < bytes_per_padded_sample; byte++) {
            sample |= ((int32_t)sample_ptr[byte]) << (byte * 8);
        }

        // Right-shift to remove padding
        sample >>= shift_amount;

        // Sign-extend to fill the container
        // Create a mask for the sign bit
        int32_t sign_bit = 1 << (bits_per_sample - 1);
        if (sample & sign_bit) {
            // Negative number - extend sign bits
            int32_t extension_mask = ~((1 << bits_per_sample) - 1);
            sample |= extension_mask;
        }

        // Write back in little-endian
        for (uint32_t byte = 0; byte < bytes_per_padded_sample; byte++) {
            output_ptr[byte] = (sample >> (byte * 8)) & 0xFF;
        }
    }
}

struct WAVHeader {
    char riff[4];           // "RIFF"
    uint32_t file_size;     // File size - 8
    char wave[4];           // "WAVE"
    char fmt[4];            // "fmt "
    uint32_t fmt_size;      // Format chunk size (16 for PCM, 40 for EXTENSIBLE)
    uint16_t audio_format;  // Audio format (1 for PCM, 0xFFFE for EXTENSIBLE)
    uint16_t num_channels;  // Number of channels
    uint32_t sample_rate;   // Sample rate
    uint32_t byte_rate;     // Byte rate
    uint16_t block_align;   // Block align
    uint16_t bits_per_sample; // Bits per sample (container size)
};

struct WAVExtensibleHeader {
    uint16_t cb_size;           // Size of extension (22)
    uint16_t valid_bits_per_sample; // Valid bits in container
    uint32_t channel_mask;      // Channel position mask
    uint8_t sub_format[16];     // Format GUID (first 2 bytes are format code)
};

struct WAVDataChunk {
    char data[4];           // "data"
    uint32_t data_size;     // Data chunk size
};

void write_wav_header(std::ofstream& file, uint32_t sample_rate, uint16_t num_channels, 
                      uint16_t bits_per_sample, uint32_t num_samples) {
    WAVHeader header;
    WAVExtensibleHeader ext_header;
    WAVDataChunk data_chunk;
    
    // Calculate container size (round up to byte boundary)
    uint16_t container_bits = ((bits_per_sample + 7) / 8) * 8;
    uint32_t bytes_per_sample = container_bits / 8;
    
    // Determine if we need WAVE_FORMAT_EXTENSIBLE
    // Use extensible format for non-standard bit depths or high bit depths
    bool use_extensible = (bits_per_sample == 12) || (bits_per_sample == 20) || 
                          (bits_per_sample == 24) || (bits_per_sample == 32) ||
                          (num_channels > 2);
    
    // RIFF chunk
    std::memcpy(header.riff, "RIFF", 4);
    
    // WAVE format
    std::memcpy(header.wave, "WAVE", 4);
    
    // fmt chunk
    std::memcpy(header.fmt, "fmt ", 4);
    
    if (use_extensible) {
        // Use WAVE_FORMAT_EXTENSIBLE
        header.fmt_size = 40;  // Extended format
        header.audio_format = 0xFFFE;  // WAVE_FORMAT_EXTENSIBLE
        header.num_channels = num_channels;
        header.sample_rate = sample_rate;
        header.bits_per_sample = container_bits;
        header.byte_rate = sample_rate * num_channels * bytes_per_sample;
        header.block_align = num_channels * bytes_per_sample;
        
        // Extension
        ext_header.cb_size = 22;
        ext_header.valid_bits_per_sample = bits_per_sample;
        ext_header.channel_mask = (num_channels == 1) ? 0x4 : 0x3;  // Mono or Stereo
        
        // GUID for PCM: {00000001-0000-0010-8000-00aa00389b71}
        static const uint8_t pcm_guid[16] = {
            0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10, 0x00,
            0x80, 0x00, 0x00, 0xaa, 0x00, 0x38, 0x9b, 0x71
        };
        std::memcpy(ext_header.sub_format, pcm_guid, 16);
    } else {
        // Use standard PCM format
        header.fmt_size = 16;  // PCM format
        header.audio_format = 1;  // PCM
        header.num_channels = num_channels;
        header.sample_rate = sample_rate;
        header.bits_per_sample = bits_per_sample;
        header.byte_rate = sample_rate * num_channels * bytes_per_sample;
        header.block_align = num_channels * bytes_per_sample;
    }
    
    // data chunk
    std::memcpy(data_chunk.data, "data", 4);
    data_chunk.data_size = num_samples * num_channels * bytes_per_sample;
    
    // Calculate file size
    uint32_t fmt_chunk_size = 8 + header.fmt_size;  // "fmt " + size + data
    uint32_t data_chunk_size = 8 + data_chunk.data_size;  // "data" + size + data
    header.file_size = 4 + fmt_chunk_size + data_chunk_size;  // "WAVE" + chunks
    
    // Write header
    file.write(reinterpret_cast<char*>(&header), sizeof(header));
    if (use_extensible) {
        file.write(reinterpret_cast<char*>(&ext_header), sizeof(ext_header));
    }
    file.write(reinterpret_cast<char*>(&data_chunk), sizeof(data_chunk));
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input.flac> <output.wav>" << std::endl;
        return 1;
    }
    
    const char* input_file = argv[1];
    const char* output_file = argv[2];
    
    // Open input FLAC file
    std::ifstream flac_file(input_file, std::ios::binary);
    if (!flac_file) {
        std::cerr << "Error: Could not open input file: " << input_file << std::endl;
        return 1;
    }

    // Initialize FLAC decoder
    esp_audio_libs::flac::FLACDecoder decoder;

    // Streaming buffer for reading file in chunks
    // This demonstrates proper streaming for embedded systems with limited RAM
    // 512KB is large enough for even high-resolution multi-channel FLAC frames
    const size_t streaming_buffer_size = 512 * 1024;  // 512KB buffer
    std::vector<uint8_t> streaming_buffer(streaming_buffer_size);
    size_t bytes_in_buffer = 0;

    // Read initial chunk
    flac_file.read(reinterpret_cast<char*>(streaming_buffer.data()), streaming_buffer_size);
    bytes_in_buffer = flac_file.gcount();

    if (bytes_in_buffer == 0) {
        std::cerr << "Error: Could not read from input file" << std::endl;
        return 1;
    }

    // Read FLAC header with streaming
    std::cout << "Reading FLAC header..." << std::endl;
    size_t header_bytes_consumed = 0;

    while (true) {
        auto result = decoder.read_header(streaming_buffer.data(), bytes_in_buffer);
        header_bytes_consumed = decoder.get_bytes_index();

        if (result == esp_audio_libs::flac::FLAC_DECODER_SUCCESS) {
            // Header complete
            break;
        } else if (result == esp_audio_libs::flac::FLAC_DECODER_HEADER_OUT_OF_DATA) {
            // Need more data - shift remaining data and read more
            size_t remaining = bytes_in_buffer - header_bytes_consumed;
            if (remaining > 0) {
                std::memmove(streaming_buffer.data(),
                            streaming_buffer.data() + header_bytes_consumed,
                            remaining);
            }
            bytes_in_buffer = remaining;

            // Read more data
            if (!flac_file.eof()) {
                flac_file.read(reinterpret_cast<char*>(streaming_buffer.data() + bytes_in_buffer),
                              streaming_buffer_size - bytes_in_buffer);
                bytes_in_buffer += flac_file.gcount();
            }

            if (bytes_in_buffer == remaining) {
                // No new data read and still out of data
                std::cerr << "Error: Unexpected end of file while reading header" << std::endl;
                return 1;
            }
            header_bytes_consumed = 0;  // Reset since we moved data
        } else {
            std::cerr << "Error: Failed to read FLAC header. Error code: " << result << std::endl;
            return 1;
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
    
    // Get audio parameters
    uint32_t sample_rate = decoder.get_sample_rate();
    uint32_t num_channels = decoder.get_num_channels();
    uint32_t bits_per_sample = decoder.get_sample_depth();
    uint32_t total_samples = decoder.get_num_samples();
    uint32_t max_block_size = decoder.get_max_block_size();
    
    std::cout << "FLAC file info:" << std::endl;
    std::cout << "  Sample rate: " << sample_rate << " Hz" << std::endl;
    std::cout << "  Channels: " << num_channels << std::endl;
    std::cout << "  Bits per sample: " << bits_per_sample << std::endl;
    std::cout << "  Total samples: " << total_samples << std::endl;
    std::cout << "  Max block size: " << max_block_size << std::endl;

    // Print MD5 signature from header
    const uint8_t* md5_sig = decoder.get_md5_signature();
    std::cout << "  MD5 signature: ";
    for (int i = 0; i < 16; i++) {
        printf("%02x", md5_sig[i]);
    }
    std::cout << std::endl;
    
    // Check for metadata
    const auto& metadata = decoder.get_metadata_blocks();
    std::cout << "  Metadata blocks: " << metadata.size() << std::endl;
    for (const auto& block : metadata) {
        std::cout << "    - Type " << block.type << ", size: " << block.length << " bytes" << std::endl;
    }
    
    // Open output WAV file
    std::ofstream wav_file(output_file, std::ios::binary);
    if (!wav_file) {
        std::cerr << "Error: Could not create output file: " << output_file << std::endl;
        return 1;
    }
    
    // Write WAV header with actual bit depth from FLAC
    // If total_samples is 0 (unknown), we'll write a placeholder and update it later
    write_wav_header(wav_file, sample_rate, num_channels, bits_per_sample, total_samples);
    
    // Remember position after header for potential update
    std::streampos header_end_pos = wav_file.tellp();
    
    // Prepare output buffer (size in bytes)
    uint32_t bytes_per_sample_out = (bits_per_sample + 7) / 8;
    std::vector<uint8_t> output_buffer(max_block_size * num_channels * bytes_per_sample_out);

    // Initialize MD5 context for verification
    MD5 md5_ctx;
    bool md5_all_zero = true;
    for (int i = 0; i < 16; i++) {
        if (md5_sig[i] != 0) {
            md5_all_zero = false;
            break;
        }
    }

    // Decode FLAC frames with streaming
    std::cout << "Decoding FLAC frames..." << std::endl;
    uint32_t frames_decoded = 0;
    uint32_t samples_per_channel_decoded = 0;

    // Refill buffer if needed before starting frame decode
    if (bytes_in_buffer < streaming_buffer_size / 2 && !flac_file.eof()) {
        flac_file.read(reinterpret_cast<char*>(streaming_buffer.data() + bytes_in_buffer),
                      streaming_buffer_size - bytes_in_buffer);
        bytes_in_buffer += flac_file.gcount();
    }

    // Decode frames with streaming
    while (bytes_in_buffer > 0) {
        uint32_t num_samples = 0;
        auto result = decoder.decode_frame(
            streaming_buffer.data(),
            bytes_in_buffer,
            output_buffer.data(),
            &num_samples
        );

        if (result == esp_audio_libs::flac::FLAC_DECODER_SUCCESS) {
            // Update MD5 with decoded samples (in little-endian byte order as per FLAC spec)
            uint32_t bytes_to_write = num_samples * bytes_per_sample_out;
            if (!md5_all_zero) {
                std::vector<uint8_t> md5_buffer(bytes_to_write);

                // For 8-bit samples, the decoder outputs unsigned (adds 128)
                // but MD5 must be computed on signed samples per FLAC spec
                if (bits_per_sample == 8) {
                    for (size_t i = 0; i < bytes_to_write; i++) {
                        md5_buffer[i] = output_buffer[i] - 128;
                    }
                } else {
                    // For non-byte-aligned samples, unpack and repack with sign-extension
                    pack_samples_for_md5(output_buffer.data(), md5_buffer.data(),
                                        num_samples, bits_per_sample);
                }

                md5_ctx.update(md5_buffer.data(), bytes_to_write);
            }

            // Write decoded samples to WAV file
            wav_file.write(reinterpret_cast<char*>(output_buffer.data()), bytes_to_write);

            // Track samples per channel for progress
            uint32_t samples_this_frame = num_samples / num_channels;
            samples_per_channel_decoded += samples_this_frame;
            frames_decoded++;

            // Update progress (only if total_samples is known)
            if (total_samples > 0 && samples_per_channel_decoded % (sample_rate * 10) == 0) {
                std::cout << "  Decoded " << samples_per_channel_decoded << " / " << total_samples
                         << " samples per channel (" << (samples_per_channel_decoded * 100 / total_samples) << "%)" << std::endl;
            }

            // Get bytes consumed and move remaining data to start of buffer
            size_t bytes_consumed = decoder.get_bytes_index();
            bytes_in_buffer -= bytes_consumed;

            if (bytes_in_buffer > 0) {
                std::memmove(streaming_buffer.data(),
                            streaming_buffer.data() + bytes_consumed,
                            bytes_in_buffer);
            }

            // Refill buffer if running low
            if (bytes_in_buffer < streaming_buffer_size / 2 && !flac_file.eof()) {
                flac_file.read(reinterpret_cast<char*>(streaming_buffer.data() + bytes_in_buffer),
                              streaming_buffer_size - bytes_in_buffer);
                bytes_in_buffer += flac_file.gcount();
            }

        } else if (result == esp_audio_libs::flac::FLAC_DECODER_NO_MORE_FRAMES) {
            std::cout << "Reached end of FLAC file." << std::endl;
            break;
        } else if (result == esp_audio_libs::flac::FLAC_DECODER_ERROR_OUT_OF_DATA) {
            // Frame incomplete - need more data
            // If buffer is full, we can't read more - this means the frame is larger than our buffer
            if (bytes_in_buffer >= streaming_buffer_size && flac_file.eof()) {
                std::cerr << "Error: Frame larger than buffer size (" << streaming_buffer_size << " bytes)" << std::endl;
                std::cerr << "  Frames decoded: " << frames_decoded << std::endl;
                std::cerr << "  Samples per channel decoded: " << samples_per_channel_decoded << " / " << total_samples << std::endl;
                wav_file.close();
                return 1;
            }

            // Try to read more data
            if (!flac_file.eof() && bytes_in_buffer < streaming_buffer_size) {
                flac_file.read(reinterpret_cast<char*>(streaming_buffer.data() + bytes_in_buffer),
                              streaming_buffer_size - bytes_in_buffer);
                size_t new_bytes = flac_file.gcount();

                if (new_bytes > 0) {
                    bytes_in_buffer += new_bytes;
                    continue;  // Try decoding again
                }
            }

            // Truly out of data
            std::cerr << "Error: Unexpected end of file while decoding frame" << std::endl;
            std::cerr << "  Frames decoded: " << frames_decoded << std::endl;
            std::cerr << "  Samples per channel decoded: " << samples_per_channel_decoded << " / " << total_samples << std::endl;
            std::cerr << "  Buffer has " << bytes_in_buffer << " bytes, EOF: " << (flac_file.eof() ? "yes" : "no") << std::endl;
            wav_file.close();
            return 1;
        } else {
            std::cerr << "Error: Failed to decode frame. Error code: " << result << std::endl;
            std::cerr << "  Frames decoded: " << frames_decoded << std::endl;
            std::cerr << "  Samples per channel decoded: " << samples_per_channel_decoded << " / " << total_samples << std::endl;
            wav_file.close();
            return 1;
        }
    }

    flac_file.close();
    
    // If total_samples was unknown (0) or incorrect, update the WAV header with actual sample count
    if (samples_per_channel_decoded != total_samples && samples_per_channel_decoded > 0) {
        wav_file.close();
        
        // Reopen file in read-write mode to update header
        std::fstream update_file(output_file, std::ios::binary | std::ios::in | std::ios::out);
        if (update_file) {
            // Update file size in RIFF header (offset 4)
            uint32_t bytes_per_sample = (bits_per_sample + 7) / 8;
            uint32_t data_size = samples_per_channel_decoded * num_channels * bytes_per_sample;
            uint32_t file_size = static_cast<uint32_t>(header_end_pos) - 8 + data_size;  // Total size minus 8 bytes for "RIFF" and size field
            
            update_file.seekp(4);
            update_file.write(reinterpret_cast<char*>(&file_size), 4);
            
            // Update data chunk size - we know it's at a fixed position for our simple WAV format
            // fmt chunk is 8 + 16 bytes, then data chunk header starts at offset 36
            // data size field is at offset 40
            update_file.seekp(40);
            update_file.write(reinterpret_cast<char*>(&data_size), 4);
            
            update_file.close();
        }
    } else {
        wav_file.close();
    }
    
    std::cout << "Successfully converted FLAC to WAV!" << std::endl;
    std::cout << "Frames decoded: " << frames_decoded << std::endl;
    std::cout << "Samples per channel decoded: " << samples_per_channel_decoded << std::endl;
    std::cout << "Output file: " << output_file << std::endl;

    // Perform MD5 verification
    std::cout << "\n=== MD5 Verification ===" << std::endl;
    if (md5_all_zero) {
        std::cout << "Status: SKIPPED (no MD5 signature in file)" << std::endl;
    } else {
        // Finalize MD5 and compare
        auto computed_md5 = md5_ctx.finalize();

        bool md5_match = true;
        for (int i = 0; i < 16; i++) {
            if (computed_md5[i] != md5_sig[i]) {
                md5_match = false;
                break;
            }
        }

        std::cout << "Expected MD5: ";
        for (int i = 0; i < 16; i++) {
            printf("%02x", md5_sig[i]);
        }
        std::cout << std::endl;

        std::cout << "Computed MD5: ";
        for (int i = 0; i < 16; i++) {
            printf("%02x", computed_md5[i]);
        }
        std::cout << std::endl;

        if (md5_match) {
            std::cout << "Result: PASS - MD5 signatures match!" << std::endl;
        } else {
            std::cout << "Result: FAIL - MD5 signatures do NOT match!" << std::endl;
        }
    }

    return 0;
}