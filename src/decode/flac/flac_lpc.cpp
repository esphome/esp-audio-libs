#include "flac_lpc.h"
#include "flac_lpc_asm.h"

#include <algorithm>
#include <climits>
#include <cstdint>

namespace esp_audio_libs {
namespace flac {

// ============================================================================
// LPC Overflow Detection Helpers
// These functions determine if 32-bit arithmetic is safe or if 64-bit is needed
// ============================================================================

// Helper function: Calculate the number of bits needed to represent a signed integer
static inline uint32_t bitmath_silog2(int64_t v) {
  if (v == 0)
    return 0;
  if (v == -1)
    return 2;

  // Handle INT64_MIN special case to avoid undefined behavior
  if (v == INT64_MIN)
    return 64;

  // Get absolute value (safe after INT64_MIN check)
  uint64_t abs_val = (v < 0) ? static_cast<uint64_t>(-v) : static_cast<uint64_t>(v);

#if defined(__GNUC__) || defined(__clang__)
  // Use builtin for counting leading zeros (more efficient)
  // __builtin_clzll returns the number of leading zero bits in a 64-bit value
  // We need to subtract from 64 to get the position of the highest set bit
  uint32_t bits = 64 - __builtin_clzll(abs_val);
#else
  // Fallback for compilers without __builtin_clzll
  uint32_t bits = 0;
  while (abs_val > 0) {
    bits++;
    abs_val >>= 1;
  }
#endif

  return bits + 1;  // Add one for sign bit
}

// Helper function: Calculate maximum prediction value before shift
static uint64_t lpc_max_prediction_value_before_shift(uint32_t subframe_bps, const int32_t *qlp_coeff, uint32_t order) {
  uint64_t max_abs_sample_value = static_cast<uint64_t>(1) << (subframe_bps - 1);
  uint32_t abs_sum_of_qlp_coeff = 0;

  for (uint32_t i = 0; i < order; i++) {
    abs_sum_of_qlp_coeff += std::abs(qlp_coeff[i]);
  }

  return max_abs_sample_value * abs_sum_of_qlp_coeff;
}

// Helper function: Calculate maximum bits needed for prediction before shift
static uint32_t lpc_max_prediction_before_shift_bps(uint32_t subframe_bps, const int32_t *qlp_coeff, uint32_t order) {
  return bitmath_silog2(lpc_max_prediction_value_before_shift(subframe_bps, qlp_coeff, order));
}

// Helper function: Calculate maximum bits needed for residual
static uint32_t lpc_max_residual_bps(uint32_t subframe_bps, const int32_t *qlp_coeff, uint32_t order,
                                     int lp_quantization) {
  uint64_t max_abs_sample_value = static_cast<uint64_t>(1) << (subframe_bps - 1);
  int64_t max_pred_before_shift =
      static_cast<int64_t>(lpc_max_prediction_value_before_shift(subframe_bps, qlp_coeff, order));
  uint64_t max_prediction_value_after_shift =
      static_cast<uint64_t>(-1 * ((-1 * max_pred_before_shift) >> lp_quantization));
  uint64_t max_residual_value = max_abs_sample_value + max_prediction_value_after_shift;
  return bitmath_silog2(max_residual_value);
}

// ============================================================================
// Public API
// ============================================================================

bool can_use_32bit_lpc(uint32_t sample_depth, const int32_t *coefs, uint32_t order, int shift) {
  return lpc_max_residual_bps(sample_depth, coefs, order, shift) <= 32 &&
         lpc_max_prediction_before_shift_bps(sample_depth, coefs, order) <= 32;
}

void restore_linear_prediction_32bit(int32_t *sub_frame_buffer, std::size_t num_of_samples,
                                     const std::vector<int32_t> &coefs, int32_t shift) {
#if (flac_lpc_asm_enabled == 1)
  // Use optimized assembly version for Xtensa
  restore_linear_prediction_32bit_asm(sub_frame_buffer, num_of_samples, coefs.data(), coefs.size(), shift);
#else
  // C++ implementation for other platforms
  const std::size_t number_of_coefficients = coefs.size();
  const std::size_t outer_loop_bound = num_of_samples - number_of_coefficients;

  for (std::size_t i = 0; i < outer_loop_bound; ++i) {
    int32_t sum = 0;
    for (std::size_t j = 0; j < number_of_coefficients; ++j) {
      sum += (sub_frame_buffer[i + j] * coefs[j]);
    }
    // Add prediction to existing residual (residual already in buffer)
    sub_frame_buffer[i + number_of_coefficients] += (sum >> shift);
  }
#endif
}

void restore_linear_prediction_64bit(int32_t *sub_frame_buffer, std::size_t num_of_samples,
                                     const std::vector<int32_t> &coefs, int32_t shift) {
#if (flac_lpc_asm_enabled == 1)
  // Use optimized 64-bit assembly version for Xtensa
  restore_linear_prediction_64bit_asm(sub_frame_buffer, num_of_samples, coefs.data(), coefs.size(), shift);
#else
  // C++ implementation for other platforms
  const std::size_t number_of_coefficients = coefs.size();
  const std::size_t outer_loop_bound = num_of_samples - number_of_coefficients;

  for (std::size_t i = 0; i < outer_loop_bound; ++i) {
    int64_t sum = 0;
    for (std::size_t j = 0; j < number_of_coefficients; ++j) {
      sum += static_cast<int64_t>(sub_frame_buffer[i + j]) * static_cast<int64_t>(coefs[j]);
    }
    // Add prediction to existing residual (residual already in buffer)
    sub_frame_buffer[i + number_of_coefficients] += static_cast<int32_t>(sum >> shift);
  }
#endif
}

}  // namespace flac
}  // namespace esp_audio_libs