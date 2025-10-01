/// @file flac_crc.h
/// @brief CRC checksum functions for FLAC frame validation
///
/// This module provides CRC-8 and CRC-16 calculation functions using lookup tables
/// for efficient validation of FLAC frame headers and frame data integrity.

#pragma once

#include <cstdint>
#include <cstddef>

namespace esp_audio_libs {
namespace flac {

/// @brief Calculate CRC-8 checksum over a buffer
///
/// Used to validate FLAC frame headers. The CRC-8 in FLAC covers the frame header
/// from the sync code through the end of the header (before audio data).
///
/// @param data Pointer to data buffer
/// @param len Length of data in bytes
/// @return Calculated CRC-8 checksum
uint8_t calculate_crc8(const uint8_t *data, std::size_t len);

/// @brief Calculate CRC-16 checksum over a buffer
///
/// Used to validate FLAC frame data integrity. The CRC-16 in FLAC covers the entire
/// frame from the sync code through all audio data (excluding the CRC-16 itself).
///
/// @param data Pointer to data buffer
/// @param len Length of data in bytes
/// @return Calculated CRC-16 checksum
uint16_t calculate_crc16(const uint8_t *data, std::size_t len);

}  // namespace flac
}  // namespace esp_audio_libs