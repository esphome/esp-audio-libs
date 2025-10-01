/// @file flac_lpc.h
/// @brief FLAC Linear Predictive Coding (LPC) helper functions
///
/// This module provides functions for LPC restoration in FLAC decoding, including
/// overflow detection and both 32-bit and 64-bit arithmetic implementations.

#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace esp_audio_libs {
namespace flac {

/// @brief Check if 32-bit arithmetic is safe for LPC restoration
///
/// Analyzes the LPC coefficients, sample depth, and shift value to determine if
/// 32-bit arithmetic can be used without overflow. If this returns false, 64-bit
/// arithmetic must be used for safe decoding.
///
/// @param sample_depth Bit depth of the audio samples (e.g., 16, 24)
/// @param coefs Pointer to array of LPC coefficients
/// @param order Number of LPC coefficients (predictor order)
/// @param shift Right shift amount applied after prediction
/// @return true if 32-bit arithmetic is safe, false if 64-bit is required
bool can_use_32bit_lpc(uint32_t sample_depth, const int32_t *coefs, uint32_t order, int shift);

/// @brief Restore linear prediction using 32-bit arithmetic
///
/// Performs LPC restoration using 32-bit integer arithmetic. This is faster but can
/// overflow with high-resolution audio or large coefficients. Use can_use_32bit_lpc()
/// to verify safety before calling this function.
///
/// @param sub_frame_buffer Buffer containing warm-up samples followed by residuals (modified in-place)
/// @param num_of_samples Total number of samples in the buffer
/// @param coefs Vector of LPC coefficients
/// @param shift Right shift amount to apply after prediction
void restore_linear_prediction_32bit(int32_t *sub_frame_buffer, std::size_t num_of_samples,
                                     const std::vector<int32_t> &coefs, int32_t shift);

/// @brief Restore linear prediction using 64-bit arithmetic
///
/// Performs LPC restoration using 64-bit integer arithmetic. This is slower but safe
/// for all valid FLAC streams, including high-resolution audio.
///
/// @param sub_frame_buffer Buffer containing warm-up samples followed by residuals (modified in-place)
/// @param num_of_samples Total number of samples in the buffer
/// @param coefs Vector of LPC coefficients
/// @param shift Right shift amount to apply after prediction
void restore_linear_prediction_64bit(int32_t *sub_frame_buffer, std::size_t num_of_samples,
                                     const std::vector<int32_t> &coefs, int32_t shift);

}  // namespace flac
}  // namespace esp_audio_libs