// A C++ class wrapper for the ART resampler

#pragma once

#include "art_biquad.h"
#include "art_resampler.h"

#include <algorithm>
#include <esp_heap_caps.h>

namespace esp_audio_libs {
namespace resampler {

// Class that wraps the ART resampler for straightforward use

struct ResamplerResults {
  size_t frames_used;
  size_t frames_generated;
  size_t predicted_frames_used;
  uint32_t clipped_samples;
};

struct ResamplerConfiguration {
  float source_sample_rate;
  float target_sample_rate;
  uint8_t source_bits_per_sample;
  uint8_t target_bits_per_sample;
  uint8_t channels;
  bool use_pre_or_post_filter;
  bool subsample_interpolate;
  uint16_t number_of_taps;
  uint16_t number_of_filters;
};

class Resampler {
 public:
  Resampler(size_t input_buffer_samples, size_t output_buffer_samples)
      : input_buffer_samples_(input_buffer_samples), output_buffer_samples_(output_buffer_samples) {}
  ~Resampler();

  /// @brief Initializes the resampler
  /// @param config ResamplerConfiguration
  /// @return True if buffers were allocated succesfully, false otherwise.
  bool initialize(ResamplerConfiguration &config);

  /// @brief Resamples the input samples to the initalized sample rate
  /// @param input Pointer to source samples as a uint8_t buffer
  /// @param output Pointer to write resampled samples as a uint8_t buffer
  /// @param input_frames_available Frames available at the input source pointer
  /// @param output_frames_free Frames free at the output sink pointer
  /// @param gain_db Gain (in dB) to apply before resampling
  /// @return (ResamplerResults) Information about the number of frames used and  generated
  ResamplerResults resample(const uint8_t *input_buffer, uint8_t *output_buffer, size_t input_frames_available,
                            size_t output_frames_free, float gain_db);

 protected:
  float *float_input_buffer_{nullptr};
  size_t input_buffer_samples_;

  float *float_output_buffer_{nullptr};
  size_t output_buffer_samples_;

  art_resampler::Resample *resampler_{nullptr};

  art_resampler::Biquad lowpass_[2][2];
  art_resampler::BiquadCoefficients lowpass_coeff_;

  uint16_t number_of_taps_;
  uint16_t number_of_filters_;

  float sample_ratio_{1.0};
  float lowpass_ratio_{1.0};

  bool pre_filter_{false};
  bool post_filter_{false};
  bool requires_resampling_{false};

  uint8_t input_bits_;
  uint8_t output_bits_;
  uint8_t channels_;
};
}  // namespace resampler
}  // namespace esp_audio_libs
