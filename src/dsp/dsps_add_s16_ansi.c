/*
 * SPDX-FileCopyrightText: 2023-2024 Espressif Systems (Shanghai) CO LTD
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#include "dsp.h"
#include <stddef.h>

esp_err_t dsps_add_s16_ansi(const int16_t *input1, const int16_t *input2, int16_t *output, int len, int step1,
                            int step2, int step_out, int shift) {
  if (NULL == input1) {
    return ESP_FAIL;
  }
  if (NULL == input2) {
    return ESP_FAIL;
  }
  if (NULL == output) {
    return ESP_FAIL;
  }

  for (int i = 0; i < len; i++) {
    int32_t acc = (int32_t) input1[i * step1] + (int32_t) input2[i * step2];
    output[i * step_out] = acc >> shift;
  }
  return ESP_OK;
}
