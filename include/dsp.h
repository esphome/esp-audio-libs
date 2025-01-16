// Copyright 2018-2019 Espressif Systems (Shanghai) PTE LTD
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef _DSP_H_
#define _DSP_H_

#include "dsp_platform.h"

#include <esp_err.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief      dot product of two float vectors
 * Dot product calculation for two floating point arrays: *dest += (src1[i] * src2[i]); i= [0..N)
 * The extension (_ansi) use ANSI C and could be compiled and run on any platform.
 * The extension (_ae32) is optimized for ESP32 chip.
 *
 * @param[in] src1  source array 1
 * @param[in] src2  source array 2
 * @param dest  destination pointer
 * @param[in] len   length of input arrays
 */
esp_err_t dsps_dotprod_f32_ansi(const float *src1, const float *src2, float *dest, int len);
esp_err_t dsps_dotprod_f32_ae32(const float *src1, const float *src2, float *dest, int len);
esp_err_t dsps_dotprod_f32_aes3(const float *src1, const float *src2, float *dest, int len);

/**
 * @brief   multiply constant
 *
 * The function multiplies input array to the constant value
 * x[i*step_out] = y[i*step_in]*C; i=[0..len)
 * The implementation use ANSI C and could be compiled and run on any platform
 *
 * @param[in] input: input array
 * @param output: output array
 * @param len: amount of operations for arrays
 * @param C: constant value
 * @param step_in: step over input array (by default should be 1)
 * @param step_out: step over output array (by default should be 1)
 *
 * @return
 *      - ESP_OK on success
 *      - One of the error codes from DSP library
 */
esp_err_t dsps_mulc_s16_ae32(const int16_t *input, int16_t *output, int len, int16_t C, int step_in, int step_out);
esp_err_t dsps_mulc_s16_ansi(const int16_t *input, int16_t *output, int len, int16_t C, int step_in, int step_out);

/**
 * @brief   add two arrays
 *
 * The function add one input array to another
 * out[i*step_out] = input1[i*step1] + input2[i*step2]; i=[0..len)
 * The implementation use ANSI C and could be compiled and run on any platform
 *
 * @param[in] input1: input array 1
 * @param[in] input2: input array 2
 * @param output: output array
 * @param len: amount of operations for arrays
 * @param step1: step over input array 1 (by default should be 1)
 * @param step2: step over input array 2 (by default should be 1)
 * @param step_out: step over output array (by default should be 1)
 *
 * @return
 *      - ESP_OK on success
 *      - One of the error codes from DSP library
 */
esp_err_t dsps_add_s16_ansi(const int16_t *input1, const int16_t *input2, int16_t *output, int len, int step1,
                            int step2, int step_out, int shift);
esp_err_t dsps_add_s16_ae32(const int16_t *input1, const int16_t *input2, int16_t *output, int len, int step1,
                            int step2, int step_out, int shift);
esp_err_t dsps_add_s16_aes3(const int16_t *input1, const int16_t *input2, int16_t *output, int len, int step1,
                            int step2, int step_out, int shift);

#if (dsps_dotprod_f32_aes3_enabled == 1)
#define dsps_dotprod_f32 dsps_dotprod_f32_aes3
#elif (dotprod_f32_ae32_enabled == 1)
#define dsps_dotprod_f32 dsps_dotprod_f32_ae32
#else
#define dsps_dotprod_f32 dsps_dotprod_f32_ansi
#endif  // dsps_dotprod_f32_ae32_enabled

#if (dsps_mulc_s16_ae32_enabled == 1)
#define dsps_mulc_s16 dsps_mulc_s16_ae32
#else
#define dsps_mulc_s16 dsps_mulc_s16_ansi
#endif  // dsps_mulc_s16_ae32_enabled

#if (dsps_add_s16_aes3_enabled == 1)
#define dsps_add_s16 dsps_add_s16_aes3
#elif (dsps_add_s16_ae32_enabled == 1)
#define dsps_add_s16 dsps_add_s16_ae32
#else
#define dsps_add_s16 dsps_add_s16_ansi
#endif  // dsps_add_s16_aes3_enabled

#ifdef __cplusplus
}
#endif

#endif  // __DSP_H__