/**
 * flac_lpc_64_asm.S
 *
 * Optimized 64-bit assembly implementation of linear prediction restoration for FLAC decoding
 * on ESP32/Xtensa architecture.
 *
 * This implementation uses MULL and MULSH instructions to perform 64-bit arithmetic
 * for high-resolution audio where 32-bit arithmetic could overflow.
 *
 * Features:
 * - 64-bit multiplication using MULL/MULSH instruction pair
 * - Jump table dispatch for orders 1-12
 * - Fully unrolled inner loops for common orders
 * - Register-cached coefficients for orders 1-4
 * - Stack-based coefficients for orders 5-12
 * - Interleaved memory operations: loads for next MAC are issued during current multiply,
 *   hiding memory latency and improving instruction throughput
 * - Generic fallback for orders > 12
 * - Proper carry propagation for 64-bit accumulation
 * - Efficient 64-bit shifting using SRC instruction
 */

#include "flac_lpc_platform.h"
#if (flac_lpc_asm_enabled == 1)

.section .text
.global restore_linear_prediction_64bit_asm
.type restore_linear_prediction_64bit_asm, @function
.align 4

/**
 * @brief Optimized 64-bit linear prediction restoration for FLAC decoding
 *
 * @param a2 Pointer to sub_frame_buffer (int32_t*)
 * @param a3 Number of samples (num_of_samples)
 * @param a4 Pointer to coefficients array (int32_t*)
 * @param a5 Number of coefficients (order)
 * @param a6 Shift amount (0-31)
 *
 * @return 0 on success, -1 on error
 *
 * Uses 64-bit arithmetic: sum = Σ(buffer[i+j] * coef[j]) >> shift
 */
restore_linear_prediction_64bit_asm:
    entry a1, 80                    // Reserve stack space for coefficients and temp storage

    // Validate inputs
    beqz a2, .L64_error            // Check for null buffer
    beqz a4, .L64_error            // Check for null coefficients
    beqz a5, .L64_done             // If order is 0, nothing to do

    // Calculate outer loop bound: num_samples - order
    sub a7, a3, a5                 // a7 = outer_loop_bound
    blti a7, 1, .L64_done          // If no iterations needed, exit

    // Jump table dispatch based on coefficient order
    movi a8, 12
    blt a8, a5, .L64_generic       // If order > 12, use generic version

    // Calculate jump address: jump_table + (order-1)*4
    addi a8, a5, -1                // a8 = order - 1
    movi a9, .L64_jump_table
    addx4 a8, a8, a9               // a8 = jump_table + (order-1)*4
    l32i a8, a8, 0                 // Load jump address
    jx a8                          // Jump to specific handler

.align 4
.L64_jump_table:
    .word .L64_order_1
    .word .L64_order_2
    .word .L64_order_3
    .word .L64_order_4
    .word .L64_order_5
    .word .L64_order_6
    .word .L64_order_7
    .word .L64_order_8
    .word .L64_order_9
    .word .L64_order_10
    .word .L64_order_11
    .word .L64_order_12

// ============================================================================
// Order 1: Single coefficient with 64-bit accumulation
// ============================================================================
.align 4
.L64_order_1:
    l32i a8, a4, 0                 // Load coefficient into register
    mov a9, a2                     // Working buffer pointer
    ssr a6                         // Set shift amount in SAR

    loop a7, .L64_loop_1_end
.L64_loop_1:
    l32i a11, a9, 0                // Load buffer[i]

    // 64-bit multiplication: buffer[i] * coef[0]
    mull a12, a11, a8              // Low 32 bits in a12
    mulsh a13, a11, a8             // High 32 bits in a13 (signed)

    // Shift 64-bit result right by 'shift' amount
    // Using SRC instruction to combine and shift
    src a14, a13, a12              // result = (high:low) >> shift

    // Add prediction to residual
    l32i a15, a9, 4                // Load buffer[i+1] (residual)
    add a15, a15, a14              // Add prediction to residual
    s32i a15, a9, 4                // Store result back

    addi a9, a9, 4                 // Advance buffer pointer
.L64_loop_1_end:
    j .L64_done

// ============================================================================
// Order 2: Two coefficients with 64-bit accumulation
// ============================================================================
.align 4
.L64_order_2:
    l32i a8, a4, 0                 // coef[0]
    l32i a9, a4, 4                 // coef[1]
    mov a10, a2                    // Working buffer pointer
    ssr a6                         // Set shift amount

    loop a7, .L64_loop_2_end
.L64_loop_2:
    // Initialize 64-bit accumulator
    movi a11, 0                    // sum_low = 0
    movi a12, 0                    // sum_high = 0

    // First: buffer[i] * coef[0]
    l32i a15, a10, 0               // buffer[i]
    mull a3, a15, a8               // Low bits (use a3 as temp)
    mulsh a14, a15, a8             // High bits (use a14 for high bits)
    add a11, a11, a3               // Add to sum_low directly in a11
    bgeu a11, a3, .L64_nc2_1      // Check carry from low addition
    addi a12, a12, 1               // Propagate carry
.L64_nc2_1:
    add a12, a12, a14              // Add high bits to sum_high

    // Second: buffer[i+1] * coef[1]
    l32i a15, a10, 4               // buffer[i+1]
    mull a3, a15, a9               // Low bits
    mulsh a14, a15, a9             // High bits
    add a11, a11, a3               // Add to sum_low directly
    bgeu a11, a3, .L64_nc2_2      // Check carry properly
    addi a12, a12, 1               // Propagate carry
.L64_nc2_2:
    add a12, a12, a14              // Add high bits

    // Shift right
    src a15, a12, a11              // result = (high:low) >> shift

    // Add to residual
    l32i a3, a10, 8                // Load residual at buffer[i+2]
    add a3, a3, a15                // Add prediction
    s32i a3, a10, 8                // Store back

    addi a10, a10, 4
.L64_loop_2_end:
    j .L64_done

// ============================================================================
// Order 3: Three coefficients with 64-bit accumulation
// ============================================================================
.align 4
.L64_order_3:
    l32i a8, a4, 0                 // coef[0]
    l32i a9, a4, 4                 // coef[1]
    l32i a13, a4, 8                // coef[2]
    mov a10, a2                    // Working buffer pointer
    ssr a6

    loop a7, .L64_loop_3_end
.L64_loop_3:
    // Initialize 64-bit accumulator
    movi a11, 0                    // sum_low = 0
    movi a12, 0                    // sum_high = 0

    // First: buffer[i] * coef[0]
    l32i a15, a10, 0
    mull a3, a15, a8               // Low bits (use a3 as temp)
    mulsh a14, a15, a8             // High bits (use a14 for high bits)
    add a11, a11, a3               // Add to sum_low directly in a11
    bgeu a11, a3, .L64_nc3_1      // Check carry from low addition
    addi a12, a12, 1               // Propagate carry
.L64_nc3_1:
    add a12, a12, a14              // Add high bits to sum_high

    // Second: buffer[i+1] * coef[1]
    l32i a15, a10, 4
    mull a3, a15, a9
    mulsh a14, a15, a9             // High bits
    add a11, a11, a3               // Add to sum_low directly
    bgeu a11, a3, .L64_nc3_2      // Check carry properly
    addi a12, a12, 1
.L64_nc3_2:
    add a12, a12, a14              // Add high bits

    // Third: buffer[i+2] * coef[2]
    l32i a15, a10, 8
    mull a3, a15, a13
    mulsh a14, a15, a13            // High bits
    add a11, a11, a3               // Add to sum_low directly
    bgeu a11, a3, .L64_nc3_3      // Check carry properly
    addi a12, a12, 1
.L64_nc3_3:
    add a12, a12, a14              // Add high bits

    // Shift right
    src a15, a12, a11

    // Add to residual
    l32i a3, a10, 12
    add a3, a3, a15
    s32i a3, a10, 12

    addi a10, a10, 4
.L64_loop_3_end:
    j .L64_done

// ============================================================================
// Order 4: Four coefficients with 64-bit accumulation
// ============================================================================
.align 4
.L64_order_4:
    l32i a8, a4, 0                 // coef[0]
    l32i a9, a4, 4                 // coef[1]
    l32i a13, a4, 8                // coef[2]
    l32i a14, a4, 12               // coef[3]
    mov a10, a2                    // Working buffer pointer
    ssr a6

    loop a7, .L64_loop_4_end
.L64_loop_4:
    // Initialize 64-bit accumulator
    movi a11, 0                    // sum_low = 0
    movi a12, 0                    // sum_high = 0

    // Accumulate all 4 products with carry propagation
    // buffer[i] * coef[0]
    l32i a15, a10, 0
    mull a3, a15, a8               // Use a3 as temp (like order 3)
    mulsh a15, a15, a8             // High bits in a15
    add a11, a11, a3
    bgeu a11, a3, .L64_nc4_1
    addi a12, a12, 1
.L64_nc4_1:
    add a12, a12, a15

    // buffer[i+1] * coef[1]
    l32i a15, a10, 4
    mull a3, a15, a9
    mulsh a15, a15, a9
    add a11, a11, a3
    bgeu a11, a3, .L64_nc4_2
    addi a12, a12, 1
.L64_nc4_2:
    add a12, a12, a15

    // buffer[i+2] * coef[2]
    l32i a15, a10, 8
    mull a3, a15, a13
    mulsh a15, a15, a13
    add a11, a11, a3
    bgeu a11, a3, .L64_nc4_3
    addi a12, a12, 1
.L64_nc4_3:
    add a12, a12, a15

    // buffer[i+3] * coef[3]
    l32i a15, a10, 12
    mull a3, a15, a14
    mulsh a15, a15, a14
    add a11, a11, a3
    bgeu a11, a3, .L64_nc4_4
    addi a12, a12, 1
.L64_nc4_4:
    add a12, a12, a15

    // Shift and store
    src a15, a12, a11

    l32i a3, a10, 16
    add a3, a3, a15
    s32i a3, a10, 16

    addi a10, a10, 4
.L64_loop_4_end:
    j .L64_done

// ============================================================================
// Order 5: Five coefficients (stack-based) with interleaved memory loads
// ============================================================================
.align 4
.L64_order_5:
    // Copy coefficients to stack
    movi a8, 0
    movi a11, 5                    // Loop count for copy loop (5 coefficients)
    loop a11, .L64_copy_5_end
.L64_copy_5:
    addx4 a9, a8, a4
    l32i a9, a9, 0
    addx4 a12, a8, a1
    s32i a9, a12, 0
    addi a8, a8, 1
.L64_copy_5_end:
    mov a10, a2
    ssr a6
    loop a7, .L64_loop_5_end
.L64_loop_5:
    movi a11, 0                    // sum_low (accumulator low)
    movi a12, 0                    // sum_high (accumulator high)

    // Pre-load first operands
    l32i a13, a10, 0               // buffer[0]
    l32i a14, a1, 0                // coefficient[0]

    // MAC 1: buffer[0] * coeff[0] - with pre-load of next
    l32i a15, a10, 4               // Pre-load buffer[1]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 4                // Pre-load coefficient[1]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc5_1      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc5_1:
    add a12, a12, a8               // Add high parts

    // MAC 2: buffer[1] * coeff[1] - using pre-loaded values
    l32i a13, a10, 8               // Pre-load buffer[2]
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    l32i a14, a1, 8                // Pre-load coefficient[2]
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc5_2      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc5_2:
    add a12, a12, a8               // Add high parts

    // MAC 3: buffer[2] * coeff[2]
    l32i a15, a10, 12              // Pre-load buffer[3]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 12               // Pre-load coefficient[3]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc5_3      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc5_3:
    add a12, a12, a8               // Add high parts

    // MAC 4: buffer[3] * coeff[3]
    l32i a13, a10, 16              // Pre-load buffer[4]
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    l32i a14, a1, 16               // Pre-load coefficient[4]
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc5_4      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc5_4:
    add a12, a12, a8               // Add high parts

    // MAC 5: buffer[4] * coeff[4] - no more pre-loads needed
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc5_5      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc5_5:
    add a12, a12, a8               // Add high parts

    // Now a11:a12 contains the 64-bit sum (low:high)
    // Perform 64-bit arithmetic right shift by amount in SAR (from a6)
    src a11, a12, a11              // Funnel shift using SAR - result in a11

    // Load existing value, add result, and store
    l32i a13, a10, 20              // Load existing value at output (5*4 = 20)
    add a13, a13, a11              // Add filtered result
    s32i a13, a10, 20              // Store back result
    addi a10, a10, 4               // Move to next sample

.L64_loop_5_end:
    j .L64_done

// ============================================================================
// Order 6: Six coefficients (stack-based)
// ============================================================================
.align 4
.L64_order_6:
    movi a8, 0
    movi a11, 6                    // Loop count for copy loop (6 coefficients)
    loop a11, .L64_copy_6_end
.L64_copy_6:
    addx4 a9, a8, a4
    l32i a9, a9, 0
    addx4 a12, a8, a1
    s32i a9, a12, 0
    addi a8, a8, 1
.L64_copy_6_end:
    mov a10, a2
    ssr a6
    loop a7, .L64_loop_6_end
.L64_loop_6:
    movi a11, 0                    // sum_low (accumulator low)
    movi a12, 0                    // sum_high (accumulator high)

    // Pre-load first operands
    l32i a13, a10, 0               // buffer[0]
    l32i a14, a1, 0                // coefficient[0]

    // MAC 1: buffer[0] * coeff[0] - with pre-load of next
    l32i a15, a10, 4               // Pre-load buffer[1]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 4                // Pre-load coefficient[1]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc6_1      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc6_1:
    add a12, a12, a8               // Add high parts

    // MAC 2: buffer[1] * coeff[1] - using pre-loaded values
    l32i a13, a10, 8               // Pre-load buffer[2]
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    l32i a14, a1, 8                // Pre-load coefficient[2]
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc6_2      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc6_2:
    add a12, a12, a8               // Add high parts

    // MAC 3: buffer[2] * coeff[2]
    l32i a15, a10, 12              // Pre-load buffer[3]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 12               // Pre-load coefficient[3]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc6_3      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc6_3:
    add a12, a12, a8               // Add high parts

    // MAC 4: buffer[3] * coeff[3]
    l32i a13, a10, 16              // Pre-load buffer[4]
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    l32i a14, a1, 16               // Pre-load coefficient[4]
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc6_4      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc6_4:
    add a12, a12, a8               // Add high parts

    // MAC 5: buffer[4] * coeff[4]
    l32i a15, a10, 20              // Pre-load buffer[5]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 20               // Pre-load coefficient[5]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc6_5      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc6_5:
    add a12, a12, a8               // Add high parts

    // MAC 6: buffer[5] * coeff[5] - no more pre-loads needed
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc6_6      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc6_6:
    add a12, a12, a8               // Add high parts

    // Now a11:a12 contains the 64-bit sum (low:high)
    // Perform 64-bit arithmetic right shift by amount in SAR (from a6)
    src a11, a12, a11              // Funnel shift using SAR - result in a11

    // Load existing value, add result, and store
    l32i a13, a10, 24              // Load existing value at output (6*4 = 24)
    add a13, a13, a11              // Add filtered result
    s32i a13, a10, 24              // Store back result
    addi a10, a10, 4               // Move to next sample

.L64_loop_6_end:
    j .L64_done

// ============================================================================
// Order 7: Seven coefficients (stack-based)
// ============================================================================
.align 4
.L64_order_7:
    movi a8, 0
    movi a11, 7                    // Loop count for copy loop (7 coefficients)
    loop a11, .L64_copy_7_end
.L64_copy_7:
    addx4 a9, a8, a4
    l32i a9, a9, 0
    addx4 a12, a8, a1
    s32i a9, a12, 0
    addi a8, a8, 1
.L64_copy_7_end:
    mov a10, a2
    ssr a6
    loop a7, .L64_loop_7_end
.L64_loop_7:
    movi a11, 0                    // sum_low (accumulator low)
    movi a12, 0                    // sum_high (accumulator high)

    // Pre-load first operands
    l32i a13, a10, 0               // buffer[0]
    l32i a14, a1, 0                // coefficient[0]

    // MAC 1: buffer[0] * coeff[0] - with pre-load of next
    l32i a15, a10, 4               // Pre-load buffer[1]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 4                // Pre-load coefficient[1]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc7_1      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc7_1:
    add a12, a12, a8               // Add high parts

    // MAC 2: buffer[1] * coeff[1] - using pre-loaded values
    l32i a13, a10, 8               // Pre-load buffer[2]
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    l32i a14, a1, 8                // Pre-load coefficient[2]
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc7_2      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc7_2:
    add a12, a12, a8               // Add high parts

    // MAC 3: buffer[2] * coeff[2]
    l32i a15, a10, 12              // Pre-load buffer[3]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 12               // Pre-load coefficient[3]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc7_3      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc7_3:
    add a12, a12, a8               // Add high parts

    // MAC 4: buffer[3] * coeff[3]
    l32i a13, a10, 16              // Pre-load buffer[4]
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    l32i a14, a1, 16               // Pre-load coefficient[4]
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc7_4      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc7_4:
    add a12, a12, a8               // Add high parts

    // MAC 5: buffer[4] * coeff[4]
    l32i a15, a10, 20              // Pre-load buffer[5]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 20               // Pre-load coefficient[5]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc7_5      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc7_5:
    add a12, a12, a8               // Add high parts

    // MAC 6: buffer[5] * coeff[5]
    l32i a13, a10, 24              // Pre-load buffer[6]
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    l32i a14, a1, 24               // Pre-load coefficient[6]
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc7_6      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc7_6:
    add a12, a12, a8               // Add high parts

    // MAC 7: buffer[6] * coeff[6] - no more pre-loads needed
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc7_7      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc7_7:
    add a12, a12, a8               // Add high parts

    // Now a11:a12 contains the 64-bit sum (low:high)
    // Perform 64-bit arithmetic right shift by amount in SAR (from a6)
    src a11, a12, a11              // Funnel shift using SAR - result in a11

    // Load existing value, add result, and store
    l32i a13, a10, 28              // Load existing value at output (7*4 = 28)
    add a13, a13, a11              // Add filtered result
    s32i a13, a10, 28              // Store back result
    addi a10, a10, 4               // Move to next sample

.L64_loop_7_end:
    j .L64_done

// ============================================================================
// Order 8: Eight coefficients (stack-based)
// ============================================================================
.align 4
.L64_order_8:
    movi a8, 0
    movi a11, 8                    // Loop count for copy loop
    loop a11, .L64_copy_8_end
.L64_copy_8:
    addx4 a9, a8, a4
    l32i a9, a9, 0
    addx4 a12, a8, a1
    s32i a9, a12, 0
    addi a8, a8, 1
.L64_copy_8_end:
    mov a10, a2
    ssr a6
    loop a7, .L64_loop_8_end
.L64_loop_8:
    movi a11, 0                    // sum_low (accumulator low)
    movi a12, 0                    // sum_high (accumulator high)

    // Pre-load first operands
    l32i a13, a10, 0               // buffer[0]
    l32i a14, a1, 0                // coefficient[0]

    // MAC 1: buffer[0] * coeff[0] - with pre-load of next
    l32i a15, a10, 4               // Pre-load buffer[1]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 4                // Pre-load coefficient[1]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc8_1      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc8_1:
    add a12, a12, a8               // Add high parts

    // MAC 2: buffer[1] * coeff[1] - using pre-loaded values
    l32i a13, a10, 8               // Pre-load buffer[2]
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    l32i a14, a1, 8                // Pre-load coefficient[2]
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc8_2      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc8_2:
    add a12, a12, a8               // Add high parts

    // MAC 3: buffer[2] * coeff[2]
    l32i a15, a10, 12              // Pre-load buffer[3]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 12               // Pre-load coefficient[3]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc8_3      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc8_3:
    add a12, a12, a8               // Add high parts

    // MAC 4: buffer[3] * coeff[3]
    l32i a13, a10, 16              // Pre-load buffer[4]
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    l32i a14, a1, 16               // Pre-load coefficient[4]
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc8_4      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc8_4:
    add a12, a12, a8               // Add high parts

    // MAC 5: buffer[4] * coeff[4]
    l32i a15, a10, 20              // Pre-load buffer[5]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 20               // Pre-load coefficient[5]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc8_5      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc8_5:
    add a12, a12, a8               // Add high parts

    // MAC 6: buffer[5] * coeff[5]
    l32i a13, a10, 24              // Pre-load buffer[6]
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    l32i a14, a1, 24               // Pre-load coefficient[6]
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc8_6      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc8_6:
    add a12, a12, a8               // Add high parts

    // MAC 7: buffer[6] * coeff[6]
    l32i a15, a10, 28              // Pre-load buffer[7]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 28               // Pre-load coefficient[7]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc8_7      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc8_7:
    add a12, a12, a8               // Add high parts

    // MAC 8: buffer[7] * coeff[7] - no more pre-loads needed
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc8_8      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc8_8:
    add a12, a12, a8               // Add high parts

    // Now a11:a12 contains the 64-bit sum (low:high)
    // Perform 64-bit arithmetic right shift by amount in SAR (from a6)
    src a11, a12, a11              // Funnel shift using SAR - result in a11

    // Load existing value, add result, and store
    l32i a13, a10, 32              // Load existing value at output
    add a13, a13, a11              // Add filtered result
    s32i a13, a10, 32              // Store back result
    addi a10, a10, 4               // Move to next sample

.L64_loop_8_end:
    j .L64_done

// ============================================================================
// Order 9: Nine coefficients (stack-based) with interleaved memory loads
// ============================================================================
.align 4
.L64_order_9:
    movi a8, 0
    movi a11, 9                    // Loop count for copy loop (9 coefficients)
    loop a11, .L64_copy_9_end
.L64_copy_9:
    addx4 a9, a8, a4
    l32i a9, a9, 0
    addx4 a12, a8, a1
    s32i a9, a12, 0
    addi a8, a8, 1
.L64_copy_9_end:
    mov a10, a2
    ssr a6
    loop a7, .L64_loop_9_end
.L64_loop_9:
    movi a11, 0                    // sum_low (accumulator low)
    movi a12, 0                    // sum_high (accumulator high)

    // Pre-load first operands
    l32i a13, a10, 0               // buffer[0]
    l32i a14, a1, 0                // coefficient[0]

    // MAC 1: buffer[0] * coeff[0] - with pre-load of next
    l32i a15, a10, 4               // Pre-load buffer[1]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 4                // Pre-load coefficient[1]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc9_1      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc9_1:
    add a12, a12, a8               // Add high parts

    // MAC 2: buffer[1] * coeff[1] - using pre-loaded values
    l32i a13, a10, 8               // Pre-load buffer[2]
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    l32i a14, a1, 8                // Pre-load coefficient[2]
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc9_2      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc9_2:
    add a12, a12, a8               // Add high parts

    // MAC 3: buffer[2] * coeff[2]
    l32i a15, a10, 12              // Pre-load buffer[3]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 12               // Pre-load coefficient[3]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc9_3      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc9_3:
    add a12, a12, a8               // Add high parts

    // MAC 4: buffer[3] * coeff[3]
    l32i a13, a10, 16              // Pre-load buffer[4]
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    l32i a14, a1, 16               // Pre-load coefficient[4]
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc9_4      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc9_4:
    add a12, a12, a8               // Add high parts

    // MAC 5: buffer[4] * coeff[4]
    l32i a15, a10, 20              // Pre-load buffer[5]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 20               // Pre-load coefficient[5]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc9_5      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc9_5:
    add a12, a12, a8               // Add high parts

    // MAC 6: buffer[5] * coeff[5]
    l32i a13, a10, 24              // Pre-load buffer[6]
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    l32i a14, a1, 24               // Pre-load coefficient[6]
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc9_6      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc9_6:
    add a12, a12, a8               // Add high parts

    // MAC 7: buffer[6] * coeff[6]
    l32i a15, a10, 28              // Pre-load buffer[7]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 28               // Pre-load coefficient[7]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc9_7      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc9_7:
    add a12, a12, a8               // Add high parts

    // MAC 8: buffer[7] * coeff[7]
    l32i a13, a10, 32              // Pre-load buffer[8]
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    l32i a14, a1, 32               // Pre-load coefficient[8]
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc9_8      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc9_8:
    add a12, a12, a8               // Add high parts

    // MAC 9: buffer[8] * coeff[8] - no more pre-loads needed
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc9_9      // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc9_9:
    add a12, a12, a8               // Add high parts

    // Now a11:a12 contains the 64-bit sum (low:high)
    // Perform 64-bit arithmetic right shift by amount in SAR (from a6)
    src a11, a12, a11              // Funnel shift using SAR - result in a11

    // Load existing value, add result, and store
    l32i a13, a10, 36              // Load existing value at output (9*4 = 36)
    add a13, a13, a11              // Add filtered result
    s32i a13, a10, 36              // Store back result
    addi a10, a10, 4               // Move to next sample

.L64_loop_9_end:
    j .L64_done

.align 4
.L64_order_10:
    movi a8, 0
    movi a11, 10                   // Loop count for copy loop (10 coefficients)
    loop a11, .L64_copy_10_end
.L64_copy_10:
    addx4 a9, a8, a4
    l32i a9, a9, 0
    addx4 a12, a8, a1
    s32i a9, a12, 0
    addi a8, a8, 1
.L64_copy_10_end:
    mov a10, a2
    ssr a6
    loop a7, .L64_loop_10_end
.L64_loop_10:
    movi a11, 0                    // sum_low (accumulator low)
    movi a12, 0                    // sum_high (accumulator high)

    // Pre-load first operands
    l32i a13, a10, 0               // buffer[0]
    l32i a14, a1, 0                // coefficient[0]

    // MAC 1: buffer[0] * coeff[0] - with pre-load of next
    l32i a15, a10, 4               // Pre-load buffer[1]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 4                // Pre-load coefficient[1]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc10_1     // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc10_1:
    add a12, a12, a8               // Add high parts

    // MAC 2: buffer[1] * coeff[1] - using pre-loaded values
    l32i a13, a10, 8               // Pre-load buffer[2]
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    l32i a14, a1, 8                // Pre-load coefficient[2]
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc10_2     // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc10_2:
    add a12, a12, a8               // Add high parts

    // MAC 3: buffer[2] * coeff[2]
    l32i a15, a10, 12              // Pre-load buffer[3]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 12               // Pre-load coefficient[3]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc10_3     // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc10_3:
    add a12, a12, a8               // Add high parts

    // MAC 4: buffer[3] * coeff[3]
    l32i a13, a10, 16              // Pre-load buffer[4]
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    l32i a14, a1, 16               // Pre-load coefficient[4]
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc10_4     // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc10_4:
    add a12, a12, a8               // Add high parts

    // MAC 5: buffer[4] * coeff[4]
    l32i a15, a10, 20              // Pre-load buffer[5]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 20               // Pre-load coefficient[5]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc10_5     // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc10_5:
    add a12, a12, a8               // Add high parts

    // MAC 6: buffer[5] * coeff[5]
    l32i a13, a10, 24              // Pre-load buffer[6]
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    l32i a14, a1, 24               // Pre-load coefficient[6]
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc10_6     // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc10_6:
    add a12, a12, a8               // Add high parts

    // MAC 7: buffer[6] * coeff[6]
    l32i a15, a10, 28              // Pre-load buffer[7]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 28               // Pre-load coefficient[7]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc10_7     // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc10_7:
    add a12, a12, a8               // Add high parts

    // MAC 8: buffer[7] * coeff[7]
    l32i a13, a10, 32              // Pre-load buffer[8]
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    l32i a14, a1, 32               // Pre-load coefficient[8]
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc10_8     // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc10_8:
    add a12, a12, a8               // Add high parts

    // MAC 9: buffer[8] * coeff[8]
    l32i a15, a10, 36              // Pre-load buffer[9]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 36               // Pre-load coefficient[9]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc10_9     // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc10_9:
    add a12, a12, a8               // Add high parts

    // MAC 10: buffer[9] * coeff[9] - no more pre-loads needed
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc10_10    // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc10_10:
    add a12, a12, a8               // Add high parts

    // Now a11:a12 contains the 64-bit sum (low:high)
    // Perform 64-bit arithmetic right shift by amount in SAR (from a6)
    src a11, a12, a11              // Funnel shift using SAR - result in a11

    // Load existing value, add result, and store
    l32i a13, a10, 40              // Load existing value at output (10*4 = 40)
    add a13, a13, a11              // Add filtered result
    s32i a13, a10, 40              // Store back result
    addi a10, a10, 4               // Move to next sample

.L64_loop_10_end:
    j .L64_done

// ============================================================================
// Order 11: Eleven coefficients (stack-based) with interleaved memory loads
// ============================================================================
.align 4
.L64_order_11:
    movi a8, 0
    movi a11, 11                   // Loop count for copy loop (11 coefficients)
    loop a11, .L64_copy_11_end
.L64_copy_11:
    addx4 a9, a8, a4
    l32i a9, a9, 0
    addx4 a12, a8, a1
    s32i a9, a12, 0
    addi a8, a8, 1
.L64_copy_11_end:
    mov a10, a2
    ssr a6
    loop a7, .L64_loop_11_end
.L64_loop_11:
    movi a11, 0                    // sum_low (accumulator low)
    movi a12, 0                    // sum_high (accumulator high)

    // Pre-load first operands
    l32i a13, a10, 0               // buffer[0]
    l32i a14, a1, 0                // coefficient[0]

    // MAC 1: buffer[0] * coeff[0] - with pre-load of next
    l32i a15, a10, 4               // Pre-load buffer[1]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 4                // Pre-load coefficient[1]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc11_1     // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc11_1:
    add a12, a12, a8               // Add high parts

    // MAC 2: buffer[1] * coeff[1] - using pre-loaded values
    l32i a13, a10, 8               // Pre-load buffer[2]
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    l32i a14, a1, 8                // Pre-load coefficient[2]
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc11_2     // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc11_2:
    add a12, a12, a8               // Add high parts

    // MAC 3: buffer[2] * coeff[2]
    l32i a15, a10, 12              // Pre-load buffer[3]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 12               // Pre-load coefficient[3]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc11_3     // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc11_3:
    add a12, a12, a8               // Add high parts

    // MAC 4: buffer[3] * coeff[3]
    l32i a13, a10, 16              // Pre-load buffer[4]
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    l32i a14, a1, 16               // Pre-load coefficient[4]
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc11_4     // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc11_4:
    add a12, a12, a8               // Add high parts

    // MAC 5: buffer[4] * coeff[4]
    l32i a15, a10, 20              // Pre-load buffer[5]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 20               // Pre-load coefficient[5]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc11_5     // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc11_5:
    add a12, a12, a8               // Add high parts

    // MAC 6: buffer[5] * coeff[5]
    l32i a13, a10, 24              // Pre-load buffer[6]
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    l32i a14, a1, 24               // Pre-load coefficient[6]
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc11_6     // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc11_6:
    add a12, a12, a8               // Add high parts

    // MAC 7: buffer[6] * coeff[6]
    l32i a15, a10, 28              // Pre-load buffer[7]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 28               // Pre-load coefficient[7]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc11_7     // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc11_7:
    add a12, a12, a8               // Add high parts

    // MAC 8: buffer[7] * coeff[7]
    l32i a13, a10, 32              // Pre-load buffer[8]
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    l32i a14, a1, 32               // Pre-load coefficient[8]
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc11_8     // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc11_8:
    add a12, a12, a8               // Add high parts

    // MAC 9: buffer[8] * coeff[8]
    l32i a15, a10, 36              // Pre-load buffer[9]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 36               // Pre-load coefficient[9]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc11_9     // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc11_9:
    add a12, a12, a8               // Add high parts

    // MAC 10: buffer[9] * coeff[9]
    l32i a13, a10, 40              // Pre-load buffer[10]
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    l32i a14, a1, 40               // Pre-load coefficient[10]
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc11_10    // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc11_10:
    add a12, a12, a8               // Add high parts

    // MAC 11: buffer[10] * coeff[10] - no more pre-loads needed
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc11_11    // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc11_11:
    add a12, a12, a8               // Add high parts

    // Now a11:a12 contains the 64-bit sum (low:high)
    // Perform 64-bit arithmetic right shift by amount in SAR (from a6)
    src a11, a12, a11              // Funnel shift using SAR - result in a11

    // Load existing value, add result, and store
    l32i a13, a10, 44              // Load existing value at output (11*4 = 44)
    add a13, a13, a11              // Add filtered result
    s32i a13, a10, 44              // Store back result
    addi a10, a10, 4               // Move to next sample

.L64_loop_11_end:
    j .L64_done

// ============================================================================
// Order 12: Twelve coefficients (stack-based) with interleaved memory loads
// ============================================================================
.align 4
.L64_order_12:
    movi a8, 0
    movi a11, 12                   // Loop count for copy loop (12 coefficients)
    loop a11, .L64_copy_12_end
.L64_copy_12:
    addx4 a9, a8, a4
    l32i a9, a9, 0
    addx4 a12, a8, a1
    s32i a9, a12, 0
    addi a8, a8, 1
.L64_copy_12_end:
    mov a10, a2
    ssr a6
    loop a7, .L64_loop_12_end
.L64_loop_12:
    movi a11, 0                    // sum_low (accumulator low)
    movi a12, 0                    // sum_high (accumulator high)

    // Pre-load first operands
    l32i a13, a10, 0               // buffer[0]
    l32i a14, a1, 0                // coefficient[0]

    // MAC 1: buffer[0] * coeff[0] - with pre-load of next
    l32i a15, a10, 4               // Pre-load buffer[1]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 4                // Pre-load coefficient[1]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc12_1     // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc12_1:
    add a12, a12, a8               // Add high parts

    // MAC 2: buffer[1] * coeff[1] - using pre-loaded values
    l32i a13, a10, 8               // Pre-load buffer[2]
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    l32i a14, a1, 8                // Pre-load coefficient[2]
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc12_2     // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc12_2:
    add a12, a12, a8               // Add high parts

    // MAC 3: buffer[2] * coeff[2]
    l32i a15, a10, 12              // Pre-load buffer[3]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 12               // Pre-load coefficient[3]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc12_3     // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc12_3:
    add a12, a12, a8               // Add high parts

    // MAC 4: buffer[3] * coeff[3]
    l32i a13, a10, 16              // Pre-load buffer[4]
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    l32i a14, a1, 16               // Pre-load coefficient[4]
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc12_4     // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc12_4:
    add a12, a12, a8               // Add high parts

    // MAC 5: buffer[4] * coeff[4]
    l32i a15, a10, 20              // Pre-load buffer[5]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 20               // Pre-load coefficient[5]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc12_5     // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc12_5:
    add a12, a12, a8               // Add high parts

    // MAC 6: buffer[5] * coeff[5]
    l32i a13, a10, 24              // Pre-load buffer[6]
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    l32i a14, a1, 24               // Pre-load coefficient[6]
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc12_6     // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc12_6:
    add a12, a12, a8               // Add high parts

    // MAC 7: buffer[6] * coeff[6]
    l32i a15, a10, 28              // Pre-load buffer[7]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 28               // Pre-load coefficient[7]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc12_7     // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc12_7:
    add a12, a12, a8               // Add high parts

    // MAC 8: buffer[7] * coeff[7]
    l32i a13, a10, 32              // Pre-load buffer[8]
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    l32i a14, a1, 32               // Pre-load coefficient[8]
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc12_8     // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc12_8:
    add a12, a12, a8               // Add high parts

    // MAC 9: buffer[8] * coeff[8]
    l32i a15, a10, 36              // Pre-load buffer[9]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 36               // Pre-load coefficient[9]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc12_9     // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc12_9:
    add a12, a12, a8               // Add high parts

    // MAC 10: buffer[9] * coeff[9]
    l32i a13, a10, 40              // Pre-load buffer[10]
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    l32i a14, a1, 40               // Pre-load coefficient[10]
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc12_10    // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc12_10:
    add a12, a12, a8               // Add high parts

    // MAC 11: buffer[10] * coeff[10]
    l32i a15, a10, 44              // Pre-load buffer[11]
    mulsh a8, a13, a14             // High 32 bits of product
    mull a13, a13, a14             // Low 32 bits of product
    l32i a14, a1, 44               // Pre-load coefficient[11]
    // 64-bit Addition: accumulator += product
    add a11, a13, a11              // Add low parts (result in a11)
    bgeu a11, a13, .L64_nc12_11    // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc12_11:
    add a12, a12, a8               // Add high parts

    // MAC 12: buffer[11] * coeff[11] - no more pre-loads needed
    mulsh a8, a15, a14             // High 32 bits of product
    mull a15, a15, a14             // Low 32 bits of product
    // 64-bit Addition: accumulator += product
    add a11, a15, a11              // Add low parts (result in a11)
    bgeu a11, a15, .L64_nc12_12    // Check for carry
    addi a12, a12, 1               // Add carry to high
.L64_nc12_12:
    add a12, a12, a8               // Add high parts

    // Now a11:a12 contains the 64-bit sum (low:high)
    // Perform 64-bit arithmetic right shift by amount in SAR (from a6)
    src a11, a12, a11              // Funnel shift using SAR - result in a11

    // Load existing value, add result, and store
    l32i a13, a10, 48              // Load existing value at output (12*4 = 48)
    add a13, a13, a11              // Add filtered result
    s32i a13, a10, 48              // Store back result
    addi a10, a10, 4               // Move to next sample

.L64_loop_12_end:
    j .L64_done

// ============================================================================
// Generic version for any order > 12 with 64-bit accumulation
// ============================================================================
.align 4
.L64_generic:
    mov a10, a2                    // Working buffer pointer
    ssr a6                         // Set shift register

    loop a7, .L64_outer_generic_end
.L64_outer_generic:
    movi a11, 0                    // sum_low = 0
    movi a12, 0                    // sum_high = 0
    movi a13, 0                    // j = 0 (inner loop counter)
    mov a14, a10                   // Current buffer position for inner loop

    loop a5, .L64_inner_generic_end
.L64_inner_generic:
    l32i a15, a14, 0               // buffer[i+j]
    addx4 a9, a13, a4              // &coef[j]
    l32i a3, a9, 0                 // coef[j] (use a3 for coefficient)

    // 64-bit multiplication and accumulation
    mull a9, a15, a3               // Low 32 bits (result in a9)
    add a9, a11, a9                // New sum_low
    bgeu a9, a11, .L64_nc_gen      // Check for carry
    addi a12, a12, 1               // Propagate carry
.L64_nc_gen:
    mov a11, a9                    // Update sum_low
    mulsh a9, a15, a3              // High 32 bits
    add a12, a12, a9               // Add to sum_high

    addi a13, a13, 1               // j++
    addi a14, a14, 4               // Advance buffer pointer for inner loop
.L64_inner_generic_end:

    // Shift 64-bit result
    src a13, a12, a11              // result = (high:low) >> shift

    // Add to residual
    addx4 a14, a5, a10             // &buffer[i + order]
    l32i a15, a14, 0               // Load residual
    add a15, a15, a13              // Add prediction
    s32i a15, a14, 0               // Store result

    addi a10, a10, 4               // Advance buffer pointer
.L64_outer_generic_end:
    j .L64_done

// ============================================================================
// Exit points
// ============================================================================
.align 4
.L64_error:
    movi a2, -1                    // Return error code
    retw.n

.align 4
.L64_done:
    movi a2, 0                     // Return success
    retw.n

#endif  // flac_lpc_asm_enabled
