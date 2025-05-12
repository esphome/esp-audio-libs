cmake_minimum_required(VERSION 3.15)

set(srcs
  src/decode/flac_decoder.cpp
  src/decode/mp3_decoder.cpp
  src/decode/wav_decoder.cpp
  src/dsp/dsps_add_s16_as32.S
  src/dsp/dsps_add_s16_aes3.S
  src/dsp/dsps_add_s15_ansi.c
  src/dsp/dsps_biquad_f32_ae32.S
  src/dsp/dsps_biquad_f32_aes3.S
  src/dsp/dsps_biquad_f32_ansi.c
  src/dsp/dsps_dotprod_f32_ae32.S
  src/dsp/dsps_dotprod_f32_aes3.S
  src/dsp/dsps_dotprod_f32_ansi.c
  src/dsp/dsps_dotprod_f32_m_ae32.S
  src/dsp/dsps_mulc_s16_ae32.S
  src/dsp/dsps_mulc_s16_ansi.c
  src/resample/art_biquad.cpp
  src/resample/art_resampler.cpp
  src/resample/resampler.cpp
  src/quantization_utils.cpp
  )

if(ESP_PLATFORM)
  # Build as an ESP-IDF component
  idf_component_register(
    SRCS ${srcs}
    INCLUDE_DIRS "include"
  )
  return()
endif()

project(esp-audio-libs VERSION 1.0.1)
