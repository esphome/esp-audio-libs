/**
 * flac_lpc_32_asm.S
 *
 * Optimized 32-bit assembly implementation of linear prediction restoration for FLAC decoding
 * on ESP32/Xtensa architecture.
 *
 * This implementation uses 32-bit integer arithmetic for efficient processing
 * of standard resolution audio where overflow is not possible.
 *
 * Features:
 * - 32-bit multiplication and accumulation
 * - Jump table dispatch for orders 1-12
 * - Fully unrolled inner loops for common orders
 * - Register-cached coefficients for orders 1-5
 * - Stack-based coefficients for orders 6-12
 * - Interleaved memory operations: loads for next MAC are issued during current multiply,
 *   hiding memory latency and improving instruction throughput
 * - Generic fallback for orders > 12
 */

#include "flac_lpc_platform.h"
#if (flac_lpc_asm_enabled == 1)

.section .text
.global restore_linear_prediction_32bit_asm
.type restore_linear_prediction_32bit_asm, @function
.align 4

/**
 * @brief Optimized 32-bit linear prediction restoration for FLAC decoding
 *
 * @param a2 Pointer to sub_frame_buffer (int32_t*)
 * @param a3 Number of samples (num_of_samples)
 * @param a4 Pointer to coefficients array (int32_t*)
 * @param a5 Number of coefficients (order)
 * @param a6 Shift amount (0-31)
 *
 * @return 0 on success, -1 on error
 *
 * Uses 32-bit arithmetic: sum = Σ(buffer[i+j] * coef[j]) >> shift
 */
restore_linear_prediction_32bit_asm:
    entry a1, 64                    // Reserve stack space for coefficients

    // Validate inputs
    beqz a2, .Lerror               // Check for null buffer
    beqz a4, .Lerror               // Check for null coefficients
    beqz a5, .Ldone                // If order is 0, nothing to do

    // Calculate outer loop bound: num_samples - order
    sub a7, a3, a5                  // a7 = outer_loop_bound
    blti a7, 1, .Ldone              // If no iterations needed (a7 <= 0), exit

    // Jump table dispatch based on coefficient order
    movi a8, 12
    blt a8, a5, .Lgeneric           // If 12 < order (order > 12), use generic version

    // Calculate jump address: jump_table + (order-1)*4
    addi a8, a5, -1                 // a8 = order - 1
    movi a9, .Ljump_table
    addx4 a8, a8, a9                // a8 = jump_table + (order-1)*4
    l32i a8, a8, 0                  // Load jump address
    jx a8                           // Jump to specific handler

.align 4
.Ljump_table:
    .word .Lorder_1
    .word .Lorder_2
    .word .Lorder_3
    .word .Lorder_4
    .word .Lorder_5
    .word .Lorder_6
    .word .Lorder_7
    .word .Lorder_8
    .word .Lorder_9
    .word .Lorder_10
    .word .Lorder_11
    .word .Lorder_12

// ============================================================================
// Order 1: Single coefficient
// ============================================================================
.align 4
.Lorder_1:
    l32i a8, a4, 0                  // Load coefficient into register
    mov a9, a2                      // Working buffer pointer
    ssr a6                          // Set shift register outside loop

    loop a7, .Lloop_1_end
.Lloop_1:
    l32i a11, a9, 0                 // Load buffer[i]
    mull a11, a11, a8               // sum = buffer[i] * coef[0]
    sra a11, a11                    // sum >>= shift
    l32i a12, a9, 4                 // Load buffer[i+1] (residual)
    add a12, a12, a11               // Add prediction to residual
    s32i a12, a9, 4                 // Store result back

    addi a9, a9, 4                  // Advance buffer pointer
.Lloop_1_end:
    j .Ldone

// ============================================================================
// Order 2: Two coefficients
// ============================================================================
.align 4
.Lorder_2:
    l32i a8, a4, 0                  // coef[0]
    l32i a9, a4, 4                  // coef[1]
    mov a10, a2                     // Working buffer pointer
    ssr a6                          // Set shift register outside loop

    loop a7, .Lloop_2_end
.Lloop_2:
    l32i a12, a10, 0                // buffer[i]
    l32i a13, a10, 4                // buffer[i+1]
    mull a12, a12, a8               // buffer[i] * coef[0]
    mull a13, a13, a9               // buffer[i+1] * coef[1]
    add a12, a12, a13               // sum

    sra a12, a12                    // sum >>= shift

    l32i a13, a10, 8                // Load residual at buffer[i+2]
    add a13, a13, a12               // Add prediction
    s32i a13, a10, 8                // Store back

    addi a10, a10, 4
.Lloop_2_end:
    j .Ldone

// ============================================================================
// Order 3: Three coefficients
// ============================================================================
.align 4
.Lorder_3:
    l32i a8, a4, 0                  // coef[0]
    l32i a9, a4, 4                  // coef[1]
    l32i a13, a4, 8                 // coef[2]
    mov a10, a2                     // Working buffer pointer
    ssr a6                          // Set shift register outside loop

    loop a7, .Lloop_3_end
.Lloop_3:
    l32i a11, a10, 0                // buffer[i]
    l32i a12, a10, 4                // buffer[i+1]
    mull a11, a11, a8
    l32i a15, a10, 8                // buffer[i+2]
    mull a12, a12, a9
    add a11, a11, a12
    mull a15, a15, a13
    add a11, a11, a15               // sum

    sra a11, a11

    l32i a12, a10, 12               // Load residual
    add a12, a12, a11
    s32i a12, a10, 12

    addi a10, a10, 4
.Lloop_3_end:
    j .Ldone

// ============================================================================
// Order 4: Four coefficients
// ============================================================================
.align 4
.Lorder_4:
    l32i a8, a4, 0                  // coef[0]
    l32i a9, a4, 4                  // coef[1]
    l32i a13, a4, 8                 // coef[2]
    l32i a14, a4, 12                // coef[3]
    mov a10, a2                     // Working buffer pointer
    ssr a6                          // Set shift register outside loop

    loop a7, .Lloop_4_end
.Lloop_4:
    // Interleaved loads and multiplies for better pipelining
    l32i a11, a10, 0
    l32i a12, a10, 4
    mull a11, a11, a8
    mull a12, a12, a9
    add a11, a11, a12

    l32i a12, a10, 8
    mull a12, a12, a13
    add a11, a11, a12

    l32i a12, a10, 12
    mull a12, a12, a14
    add a11, a11, a12               // Final sum

    sra a11, a11

    l32i a12, a10, 16
    add a12, a12, a11
    s32i a12, a10, 16

    addi a10, a10, 4
.Lloop_4_end:
    j .Ldone

// ============================================================================
// Order 5: Five coefficients
// ============================================================================
.align 4
.Lorder_5:
    // Load all 5 coefficients into registers
    l32i a8, a4, 0
    l32i a9, a4, 4
    l32i a13, a4, 8
    l32i a14, a4, 12
    l32i a15, a4, 16
    mov a10, a2                     // Working buffer pointer
    ssr a6                          // Set shift register outside loop

    loop a7, .Lloop_5_end
.Lloop_5:
    l32i a11, a10, 0
    mull a11, a11, a8

    // a3 is now free to use as temp register since num_samples is no longer needed
    l32i a3, a10, 4                 // Temporarily use a3
    mull a3, a3, a9
    add a11, a11, a3

    l32i a3, a10, 8
    mull a3, a3, a13
    add a11, a11, a3

    l32i a3, a10, 12
    mull a3, a3, a14
    add a11, a11, a3

    l32i a3, a10, 16
    mull a3, a3, a15
    add a11, a11, a3

    sra a11, a11

    l32i a3, a10, 20
    add a3, a3, a11
    s32i a3, a10, 20

    addi a10, a10, 4
.Lloop_5_end:
    j .Ldone

// ============================================================================
// Order 6: Six coefficients (stack-based)
// ============================================================================
.align 4
.Lorder_6:
    // Copy all coefficients to stack
    l32i a8, a4, 0
    s32i a8, a1, 0
    l32i a8, a4, 4
    s32i a8, a1, 4
    l32i a8, a4, 8
    s32i a8, a1, 8
    l32i a8, a4, 12
    s32i a8, a1, 12
    l32i a8, a4, 16
    s32i a8, a1, 16
    l32i a8, a4, 20
    s32i a8, a1, 20

    mov a10, a2                     // Working buffer pointer
    ssr a6                          // Set shift register outside loop

    loop a7, .Lloop_6_end
.Lloop_6:
    movi a11, 0                     // sum = 0

    // Fully unrolled inner loop
    l32i a12, a10, 0
    l32i a13, a1, 0
    mull a12, a12, a13
    add a11, a11, a12

    l32i a12, a10, 4
    l32i a13, a1, 4
    mull a12, a12, a13
    add a11, a11, a12

    l32i a12, a10, 8
    l32i a13, a1, 8
    mull a12, a12, a13
    add a11, a11, a12

    l32i a12, a10, 12
    l32i a13, a1, 12
    mull a12, a12, a13
    add a11, a11, a12

    l32i a12, a10, 16
    l32i a13, a1, 16
    mull a12, a12, a13
    add a11, a11, a12

    l32i a12, a10, 20
    l32i a13, a1, 20
    mull a12, a12, a13
    add a11, a11, a12

    sra a11, a11

    l32i a12, a10, 24
    add a12, a12, a11
    s32i a12, a10, 24

    addi a10, a10, 4
.Lloop_6_end:
    j .Ldone

// ============================================================================
// Order 7: Seven coefficients (stack-based)
// ============================================================================
.align 4
.Lorder_7:
    // Copy all coefficients to stack
    movi a8, 0
    movi a11, 7                     // Load count for copy loop (7 coefficients)
    loop a11, .Lcopy_7_end
.Lcopy_7:
    addx4 a9, a8, a4
    l32i a9, a9, 0
    addx4 a12, a8, a1
    s32i a9, a12, 0
    addi a8, a8, 1
.Lcopy_7_end:
    mov a10, a2
    ssr a6                          // Set shift register outside loop
    loop a7, .Lloop_7_end
.Lloop_7:
    movi a11, 0                     // Initialize accumulator

    // Pre-load first operands
    l32i a12, a10, 0                // buffer[0]
    l32i a13, a1, 0                 // coefficient[0]

    // MAC 1: buffer[0] * coeff[0] with pre-loads
    l32i a14, a10, 4                // Pre-load buffer[1]
    mull a12, a12, a13              // Multiply
    l32i a13, a1, 4                 // Pre-load coefficient[1]
    add a11, a11, a12               // Accumulate

    // MAC 2: buffer[1] * coeff[1]
    l32i a12, a10, 8                // Pre-load buffer[2]
    mull a14, a14, a13              // Multiply
    l32i a13, a1, 8                 // Pre-load coefficient[2]
    add a11, a11, a14               // Accumulate

    // MAC 3: buffer[2] * coeff[2]
    l32i a14, a10, 12               // Pre-load buffer[3]
    mull a12, a12, a13              // Multiply
    l32i a13, a1, 12                // Pre-load coefficient[3]
    add a11, a11, a12               // Accumulate

    // MAC 4: buffer[3] * coeff[3]
    l32i a12, a10, 16               // Pre-load buffer[4]
    mull a14, a14, a13              // Multiply
    l32i a13, a1, 16                // Pre-load coefficient[4]
    add a11, a11, a14               // Accumulate

    // MAC 5: buffer[4] * coeff[4]
    l32i a14, a10, 20               // Pre-load buffer[5]
    mull a12, a12, a13              // Multiply
    l32i a13, a1, 20                // Pre-load coefficient[5]
    add a11, a11, a12               // Accumulate

    // MAC 6: buffer[5] * coeff[5]
    l32i a12, a10, 24               // Pre-load buffer[6]
    mull a14, a14, a13              // Multiply
    l32i a13, a1, 24                // Pre-load coefficient[6]
    add a11, a11, a14               // Accumulate

    // MAC 7: buffer[6] * coeff[6] - no more pre-loads
    mull a12, a12, a13              // Multiply
    add a11, a11, a12               // Accumulate

    // Shift and store result
    sra a11, a11                    // Arithmetic right shift using SAR
    l32i a12, a10, 28               // Load existing value at offset 28 (7*4)
    add a12, a12, a11               // Add filtered result
    s32i a12, a10, 28               // Store back at offset 28
    addi a10, a10, 4                // Move to next sample

.Lloop_7_end:
    j .Ldone

// ============================================================================
// Order 8: Eight coefficients (stack-based)
// ============================================================================
.align 4
.Lorder_8:
    // Copy all coefficients to stack
    movi a8, 0
    movi a11, 8                     // Load count for copy loop
    loop a11, .Lcopy_8_end
.Lcopy_8:
    addx4 a9, a8, a4
    l32i a9, a9, 0
    addx4 a12, a8, a1
    s32i a9, a12, 0
    addi a8, a8, 1
.Lcopy_8_end:
    mov a10, a2
    ssr a6                          // Set shift register outside loop
    loop a7, .Lloop_8_end
.Lloop_8:
    movi a11, 0                     // Initialize accumulator

    // Pre-load first operands
    l32i a12, a10, 0                // buffer[0]
    l32i a13, a1, 0                 // coefficient[0]

    // MAC 1: buffer[0] * coeff[0] with pre-loads
    l32i a14, a10, 4                // Pre-load buffer[1]
    mull a12, a12, a13              // Multiply
    l32i a13, a1, 4                 // Pre-load coefficient[1]
    add a11, a11, a12               // Accumulate

    // MAC 2: buffer[1] * coeff[1]
    l32i a12, a10, 8                // Pre-load buffer[2]
    mull a14, a14, a13              // Multiply
    l32i a13, a1, 8                 // Pre-load coefficient[2]
    add a11, a11, a14               // Accumulate

    // MAC 3: buffer[2] * coeff[2]
    l32i a14, a10, 12               // Pre-load buffer[3]
    mull a12, a12, a13              // Multiply
    l32i a13, a1, 12                // Pre-load coefficient[3]
    add a11, a11, a12               // Accumulate

    // MAC 4: buffer[3] * coeff[3]
    l32i a12, a10, 16               // Pre-load buffer[4]
    mull a14, a14, a13              // Multiply
    l32i a13, a1, 16                // Pre-load coefficient[4]
    add a11, a11, a14               // Accumulate

    // MAC 5: buffer[4] * coeff[4]
    l32i a14, a10, 20               // Pre-load buffer[5]
    mull a12, a12, a13              // Multiply
    l32i a13, a1, 20                // Pre-load coefficient[5]
    add a11, a11, a12               // Accumulate

    // MAC 6: buffer[5] * coeff[5]
    l32i a12, a10, 24               // Pre-load buffer[6]
    mull a14, a14, a13              // Multiply
    l32i a13, a1, 24                // Pre-load coefficient[6]
    add a11, a11, a14               // Accumulate

    // MAC 7: buffer[6] * coeff[6]
    l32i a14, a10, 28               // Pre-load buffer[7]
    mull a12, a12, a13              // Multiply
    l32i a13, a1, 28                // Pre-load coefficient[7]
    add a11, a11, a12               // Accumulate

    // MAC 8: buffer[7] * coeff[7] - no more pre-loads
    mull a14, a14, a13              // Multiply
    add a11, a11, a14               // Accumulate

    // Shift and store result
    sra a11, a11                    // Arithmetic right shift using SAR
    l32i a12, a10, 32               // Load existing value
    add a12, a12, a11               // Add filtered result
    s32i a12, a10, 32               // Store back
    addi a10, a10, 4                // Move to next sample

.Lloop_8_end:
    j .Ldone

// ============================================================================
// Order 9: Nine coefficients (stack-based)
// ============================================================================
.align 4
.Lorder_9:
    // Copy all coefficients to stack
    movi a8, 0
    movi a11, 9                     // Load count for copy loop
    loop a11, .Lcopy_9_end
.Lcopy_9:
    addx4 a9, a8, a4
    l32i a9, a9, 0
    addx4 a12, a8, a1
    s32i a9, a12, 0
    addi a8, a8, 1
.Lcopy_9_end:
    mov a10, a2
    ssr a6                          // Set shift register outside loop
    loop a7, .Lloop_9_end
.Lloop_9:
    movi a11, 0                     // Initialize accumulator

    // Pre-load first operands
    l32i a12, a10, 0                // buffer[0]
    l32i a13, a1, 0                 // coefficient[0]

    // MAC 1: buffer[0] * coeff[0] with pre-loads
    l32i a14, a10, 4                // Pre-load buffer[1]
    mull a12, a12, a13              // Multiply
    l32i a13, a1, 4                 // Pre-load coefficient[1]
    add a11, a11, a12               // Accumulate

    // MAC 2: buffer[1] * coeff[1]
    l32i a12, a10, 8                // Pre-load buffer[2]
    mull a14, a14, a13              // Multiply
    l32i a13, a1, 8                 // Pre-load coefficient[2]
    add a11, a11, a14               // Accumulate

    // MAC 3: buffer[2] * coeff[2]
    l32i a14, a10, 12               // Pre-load buffer[3]
    mull a12, a12, a13              // Multiply
    l32i a13, a1, 12                // Pre-load coefficient[3]
    add a11, a11, a12               // Accumulate

    // MAC 4: buffer[3] * coeff[3]
    l32i a12, a10, 16               // Pre-load buffer[4]
    mull a14, a14, a13              // Multiply
    l32i a13, a1, 16                // Pre-load coefficient[4]
    add a11, a11, a14               // Accumulate

    // MAC 5: buffer[4] * coeff[4]
    l32i a14, a10, 20               // Pre-load buffer[5]
    mull a12, a12, a13              // Multiply
    l32i a13, a1, 20                // Pre-load coefficient[5]
    add a11, a11, a12               // Accumulate

    // MAC 6: buffer[5] * coeff[5]
    l32i a12, a10, 24               // Pre-load buffer[6]
    mull a14, a14, a13              // Multiply
    l32i a13, a1, 24                // Pre-load coefficient[6]
    add a11, a11, a14               // Accumulate

    // MAC 7: buffer[6] * coeff[6]
    l32i a14, a10, 28               // Pre-load buffer[7]
    mull a12, a12, a13              // Multiply
    l32i a13, a1, 28                // Pre-load coefficient[7]
    add a11, a11, a12               // Accumulate

    // MAC 8: buffer[7] * coeff[7]
    l32i a12, a10, 32               // Pre-load buffer[8]
    mull a14, a14, a13              // Multiply
    l32i a13, a1, 32                // Pre-load coefficient[8]
    add a11, a11, a14               // Accumulate

    // MAC 9: buffer[8] * coeff[8] - no more pre-loads
    mull a12, a12, a13              // Multiply
    add a11, a11, a12               // Accumulate

    // Shift and store result
    sra a11, a11                    // Arithmetic right shift using SAR
    l32i a12, a10, 36               // Load existing value at offset 36 (9*4)
    add a12, a12, a11               // Add filtered result
    s32i a12, a10, 36               // Store back at offset 36
    addi a10, a10, 4                // Move to next sample

.Lloop_9_end:
    j .Ldone

// ============================================================================
// Order 10: Ten coefficients (stack-based)
// ============================================================================
.align 4
.Lorder_10:
    // Copy all coefficients to stack
    movi a8, 0
    movi a11, 10                    // Load count for copy loop
    loop a11, .Lcopy_10_end
.Lcopy_10:
    addx4 a9, a8, a4
    l32i a9, a9, 0
    addx4 a12, a8, a1
    s32i a9, a12, 0
    addi a8, a8, 1
.Lcopy_10_end:
    mov a10, a2
    ssr a6                          // Set shift register outside loop
    loop a7, .Lloop_10_end
.Lloop_10:
    movi a11, 0                     // Initialize accumulator

    // Pre-load first operands
    l32i a12, a10, 0                // buffer[0]
    l32i a13, a1, 0                 // coefficient[0]

    // MAC 1: buffer[0] * coeff[0] with pre-loads
    l32i a14, a10, 4                // Pre-load buffer[1]
    mull a12, a12, a13              // Multiply
    l32i a13, a1, 4                 // Pre-load coefficient[1]
    add a11, a11, a12               // Accumulate

    // MAC 2: buffer[1] * coeff[1]
    l32i a12, a10, 8                // Pre-load buffer[2]
    mull a14, a14, a13              // Multiply
    l32i a13, a1, 8                 // Pre-load coefficient[2]
    add a11, a11, a14               // Accumulate

    // MAC 3: buffer[2] * coeff[2]
    l32i a14, a10, 12               // Pre-load buffer[3]
    mull a12, a12, a13              // Multiply
    l32i a13, a1, 12                // Pre-load coefficient[3]
    add a11, a11, a12               // Accumulate

    // MAC 4: buffer[3] * coeff[3]
    l32i a12, a10, 16               // Pre-load buffer[4]
    mull a14, a14, a13              // Multiply
    l32i a13, a1, 16                // Pre-load coefficient[4]
    add a11, a11, a14               // Accumulate

    // MAC 5: buffer[4] * coeff[4]
    l32i a14, a10, 20               // Pre-load buffer[5]
    mull a12, a12, a13              // Multiply
    l32i a13, a1, 20                // Pre-load coefficient[5]
    add a11, a11, a12               // Accumulate

    // MAC 6: buffer[5] * coeff[5]
    l32i a12, a10, 24               // Pre-load buffer[6]
    mull a14, a14, a13              // Multiply
    l32i a13, a1, 24                // Pre-load coefficient[6]
    add a11, a11, a14               // Accumulate

    // MAC 7: buffer[6] * coeff[6]
    l32i a14, a10, 28               // Pre-load buffer[7]
    mull a12, a12, a13              // Multiply
    l32i a13, a1, 28                // Pre-load coefficient[7]
    add a11, a11, a12               // Accumulate

    // MAC 8: buffer[7] * coeff[7]
    l32i a12, a10, 32               // Pre-load buffer[8]
    mull a14, a14, a13              // Multiply
    l32i a13, a1, 32                // Pre-load coefficient[8]
    add a11, a11, a14               // Accumulate

    // MAC 9: buffer[8] * coeff[8]
    l32i a14, a10, 36               // Pre-load buffer[9]
    mull a12, a12, a13              // Multiply
    l32i a13, a1, 36                // Pre-load coefficient[9]
    add a11, a11, a12               // Accumulate

    // MAC 10: buffer[9] * coeff[9] - no more pre-loads
    mull a14, a14, a13              // Multiply
    add a11, a11, a14               // Accumulate

    // Shift and store result
    sra a11, a11                    // Arithmetic right shift using SAR
    l32i a12, a10, 40               // Load existing value at offset 40 (10*4)
    add a12, a12, a11               // Add filtered result
    s32i a12, a10, 40               // Store back at offset 40
    addi a10, a10, 4                // Move to next sample

.Lloop_10_end:
    j .Ldone

// ============================================================================
// Order 11: Eleven coefficients (stack-based)
// ============================================================================
.align 4
.Lorder_11:
    // Copy all coefficients to stack
    movi a8, 0
    movi a11, 11                    // Load count for copy loop
    loop a11, .Lcopy_11_end
.Lcopy_11:
    addx4 a9, a8, a4
    l32i a9, a9, 0
    addx4 a12, a8, a1
    s32i a9, a12, 0
    addi a8, a8, 1
.Lcopy_11_end:
    mov a10, a2
    ssr a6                          // Set shift register outside loop
    loop a7, .Lloop_11_end
.Lloop_11:
    movi a11, 0                     // Initialize accumulator

    // Pre-load first operands
    l32i a12, a10, 0                // buffer[0]
    l32i a13, a1, 0                 // coefficient[0]

    // MAC 1: buffer[0] * coeff[0] with pre-loads
    l32i a14, a10, 4                // Pre-load buffer[1]
    mull a12, a12, a13              // Multiply
    l32i a13, a1, 4                 // Pre-load coefficient[1]
    add a11, a11, a12               // Accumulate

    // MAC 2: buffer[1] * coeff[1]
    l32i a12, a10, 8                // Pre-load buffer[2]
    mull a14, a14, a13              // Multiply
    l32i a13, a1, 8                 // Pre-load coefficient[2]
    add a11, a11, a14               // Accumulate

    // MAC 3: buffer[2] * coeff[2]
    l32i a14, a10, 12               // Pre-load buffer[3]
    mull a12, a12, a13              // Multiply
    l32i a13, a1, 12                // Pre-load coefficient[3]
    add a11, a11, a12               // Accumulate

    // MAC 4: buffer[3] * coeff[3]
    l32i a12, a10, 16               // Pre-load buffer[4]
    mull a14, a14, a13              // Multiply
    l32i a13, a1, 16                // Pre-load coefficient[4]
    add a11, a11, a14               // Accumulate

    // MAC 5: buffer[4] * coeff[4]
    l32i a14, a10, 20               // Pre-load buffer[5]
    mull a12, a12, a13              // Multiply
    l32i a13, a1, 20                // Pre-load coefficient[5]
    add a11, a11, a12               // Accumulate

    // MAC 6: buffer[5] * coeff[5]
    l32i a12, a10, 24               // Pre-load buffer[6]
    mull a14, a14, a13              // Multiply
    l32i a13, a1, 24                // Pre-load coefficient[6]
    add a11, a11, a14               // Accumulate

    // MAC 7: buffer[6] * coeff[6]
    l32i a14, a10, 28               // Pre-load buffer[7]
    mull a12, a12, a13              // Multiply
    l32i a13, a1, 28                // Pre-load coefficient[7]
    add a11, a11, a12               // Accumulate

    // MAC 8: buffer[7] * coeff[7]
    l32i a12, a10, 32               // Pre-load buffer[8]
    mull a14, a14, a13              // Multiply
    l32i a13, a1, 32                // Pre-load coefficient[8]
    add a11, a11, a14               // Accumulate

    // MAC 9: buffer[8] * coeff[8]
    l32i a14, a10, 36               // Pre-load buffer[9]
    mull a12, a12, a13              // Multiply
    l32i a13, a1, 36                // Pre-load coefficient[9]
    add a11, a11, a12               // Accumulate

    // MAC 10: buffer[9] * coeff[9]
    l32i a12, a10, 40               // Pre-load buffer[10]
    mull a14, a14, a13              // Multiply
    l32i a13, a1, 40                // Pre-load coefficient[10]
    add a11, a11, a14               // Accumulate

    // MAC 11: buffer[10] * coeff[10] - no more pre-loads
    mull a12, a12, a13              // Multiply
    add a11, a11, a12               // Accumulate

    // Shift and store result
    sra a11, a11                    // Arithmetic right shift using SAR
    l32i a12, a10, 44               // Load existing value at offset 44 (11*4)
    add a12, a12, a11               // Add filtered result
    s32i a12, a10, 44               // Store back at offset 44
    addi a10, a10, 4                // Move to next sample

.Lloop_11_end:
    j .Ldone

// ============================================================================
// Order 12: Twelve coefficients (stack-based)
// ============================================================================
.align 4
.Lorder_12:
    // Copy all coefficients to stack
    movi a8, 0
    movi a11, 12                    // Load count for copy loop
    loop a11, .Lcopy_12_end
.Lcopy_12:
    addx4 a9, a8, a4
    l32i a9, a9, 0
    addx4 a12, a8, a1
    s32i a9, a12, 0
    addi a8, a8, 1
.Lcopy_12_end:
    mov a10, a2
    ssr a6                          // Set shift register outside loop
    loop a7, .Lloop_12_end
.Lloop_12:
    movi a11, 0                     // Initialize accumulator

    // Pre-load first operands
    l32i a12, a10, 0                // buffer[0]
    l32i a13, a1, 0                 // coefficient[0]

    // MAC 1: buffer[0] * coeff[0] with pre-loads
    l32i a14, a10, 4                // Pre-load buffer[1]
    mull a12, a12, a13              // Multiply
    l32i a13, a1, 4                 // Pre-load coefficient[1]
    add a11, a11, a12               // Accumulate

    // MAC 2: buffer[1] * coeff[1]
    l32i a12, a10, 8                // Pre-load buffer[2]
    mull a14, a14, a13              // Multiply
    l32i a13, a1, 8                 // Pre-load coefficient[2]
    add a11, a11, a14               // Accumulate

    // MAC 3: buffer[2] * coeff[2]
    l32i a14, a10, 12               // Pre-load buffer[3]
    mull a12, a12, a13              // Multiply
    l32i a13, a1, 12                // Pre-load coefficient[3]
    add a11, a11, a12               // Accumulate

    // MAC 4: buffer[3] * coeff[3]
    l32i a12, a10, 16               // Pre-load buffer[4]
    mull a14, a14, a13              // Multiply
    l32i a13, a1, 16                // Pre-load coefficient[4]
    add a11, a11, a14               // Accumulate

    // MAC 5: buffer[4] * coeff[4]
    l32i a14, a10, 20               // Pre-load buffer[5]
    mull a12, a12, a13              // Multiply
    l32i a13, a1, 20                // Pre-load coefficient[5]
    add a11, a11, a12               // Accumulate

    // MAC 6: buffer[5] * coeff[5]
    l32i a12, a10, 24               // Pre-load buffer[6]
    mull a14, a14, a13              // Multiply
    l32i a13, a1, 24                // Pre-load coefficient[6]
    add a11, a11, a14               // Accumulate

    // MAC 7: buffer[6] * coeff[6]
    l32i a14, a10, 28               // Pre-load buffer[7]
    mull a12, a12, a13              // Multiply
    l32i a13, a1, 28                // Pre-load coefficient[7]
    add a11, a11, a12               // Accumulate

    // MAC 8: buffer[7] * coeff[7]
    l32i a12, a10, 32               // Pre-load buffer[8]
    mull a14, a14, a13              // Multiply
    l32i a13, a1, 32                // Pre-load coefficient[8]
    add a11, a11, a14               // Accumulate

    // MAC 9: buffer[8] * coeff[8]
    l32i a14, a10, 36               // Pre-load buffer[9]
    mull a12, a12, a13              // Multiply
    l32i a13, a1, 36                // Pre-load coefficient[9]
    add a11, a11, a12               // Accumulate

    // MAC 10: buffer[9] * coeff[9]
    l32i a12, a10, 40               // Pre-load buffer[10]
    mull a14, a14, a13              // Multiply
    l32i a13, a1, 40                // Pre-load coefficient[10]
    add a11, a11, a14               // Accumulate

    // MAC 11: buffer[10] * coeff[10]
    l32i a14, a10, 44               // Pre-load buffer[11]
    mull a12, a12, a13              // Multiply
    l32i a13, a1, 44                // Pre-load coefficient[11]
    add a11, a11, a12               // Accumulate

    // MAC 12: buffer[11] * coeff[11] - no more pre-loads
    mull a14, a14, a13              // Multiply
    add a11, a11, a14               // Accumulate

    // Shift and store result
    sra a11, a11                    // Arithmetic right shift using SAR
    l32i a12, a10, 48               // Load existing value at offset 48 (12*4)
    add a12, a12, a11               // Add filtered result
    s32i a12, a10, 48               // Store back at offset 48
    addi a10, a10, 4                // Move to next sample

.Lloop_12_end:
    j .Ldone

// ============================================================================
// Generic version for any order > 12
// ============================================================================
.align 4
.Lgeneric:
    mov a10, a2                     // Working buffer pointer
    ssr a6                          // Set shift register outside loops

    loop a7, .Louter_generic_end
.Louter_generic:
    movi a11, 0                     // sum = 0
    movi a12, 0                     // j = 0 (inner loop counter)
    mov a13, a10                    // Current buffer position for inner loop

    loop a5, .Linner_generic_end
.Linner_generic:
    l32i a14, a13, 0                // buffer[i+j]
    addx4 a15, a12, a4              // &coef[j]
    l32i a15, a15, 0                // coef[j]
    mull a14, a14, a15              // buffer[i+j] * coef[j]
    add a11, a11, a14               // sum += ...

    addi a12, a12, 1                // j++
    addi a13, a13, 4                // Advance buffer pointer for inner loop
.Linner_generic_end:

    sra a11, a11                    // sum >>= shift

    addx4 a12, a5, a10              // &buffer[i + order]
    l32i a13, a12, 0                // Load residual
    add a13, a13, a11               // Add prediction
    s32i a13, a12, 0                // Store result

    addi a10, a10, 4                // Advance buffer pointer
.Louter_generic_end:
    j .Ldone

// ============================================================================
// Exit points
// ============================================================================
.align 4
.Lerror:
    movi a2, -1                     // Return error code
    retw.n

.align 4
.Ldone:
    movi a2, 0                      // Return success
    retw.n

#endif  // flac_lpc_asm_enabled
