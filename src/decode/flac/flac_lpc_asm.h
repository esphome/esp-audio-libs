#ifndef _FLAC_LPC_ASM_H_
#define _FLAC_LPC_ASM_H_

#include "flac_lpc_platform.h"
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#if (flac_lpc_asm_enabled == 1)

/// @brief Optimized 32-bit assembly implementation of linear prediction restoration for FLAC decoding
///
/// @param buffer Pointer to sub_frame_buffer (int32_t*) - contains warm-up samples followed by residuals
/// @param num_samples Total number of samples in the buffer
/// @param coefficients Pointer to LPC coefficients array (int32_t*)
/// @param order Number of coefficients (predictor order)
/// @param shift Shift amount for the prediction
///
/// @return 0 on success, -1 on error (null pointers or invalid parameters)
///
/// @note This function modifies the buffer in-place, replacing residuals with
///       restored sample values.
/// @note Optimized for orders 1-12 with fully unrolled loops; uses a generic
///       loop for higher orders.
int restore_linear_prediction_32bit_asm(
    int32_t* buffer,
    size_t num_samples,
    const int32_t* coefficients,
    size_t order,
    int32_t shift
);

/// @brief Optimized 64-bit assembly implementation of linear prediction restoration for FLAC decoding
///
/// This function uses MULL and MULSH instructions to perform 64-bit arithmetic
/// for high-resolution audio where 32-bit arithmetic could overflow.
///
/// @param buffer Pointer to the sub-frame buffer (int32_t*). The first 'order' samples
///               should contain the warm-up samples, and the remaining samples contain
///               the residuals that will be restored in-place.
/// @param num_samples Total number of samples in the buffer (including warm-up samples)
/// @param coefficients Pointer to the array of LPC coefficients
/// @param order Number of coefficients (predictor order)
/// @param shift Right shift amount to apply after accumulation (0-31)
///
/// @return 0 on success, -1 on error (null pointers)
///
/// @note This function modifies the buffer in-place, replacing residuals with
///       restored sample values.
/// @note Optimized for orders 1-12 with fully unrolled loops; uses a generic
///       loop for higher orders.
int restore_linear_prediction_64bit_asm(
    int32_t* buffer,
    size_t num_samples,
    const int32_t* coefficients,
    size_t order,
    int32_t shift
);

#endif // flac_lpc_asm_enabled

#ifdef __cplusplus
}
#endif

#endif // _FLAC_LPC_ASM_H_