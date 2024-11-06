// A C++ class wrapper for the ART resampler

#pragma once

#include "art_resampler.h"
#include "art_biquad.h"

#include <algorithm>
#include <esp_heap_caps.h>

namespace resampler {

// Class that wraps the ART resampler for straightforward use

class Resampler {
 public:
  Resampler(size_t input_buffer_samples, size_t output_buffer_samples)
      : input_buffer_samples_(input_buffer_samples), output_buffer_samples_(output_buffer_samples) {}
  ~Resampler();

  /// @brief Initializes the resampler
  /// @param target_sample_rate Output sample rate
  /// @param source_sample_rate Input sample rate
  /// @param channels Source audio's number of channels
  /// @param number_of_taps Number of taps per filter. Must be a multiple of 4.
  /// @param number_of_filters Number of windowed-sinc filters. Must be greater than 1.
  /// @param use_pre_post_filter Enable a cascading biquad filter before or after the FIR filter
  /// @param subsample_interpolate Linearly interpolate the output sample with the two nearest sinc filters
  /// @return True if buffers were allocated succesfully, false otherwise.
  bool initialize(float target_sample_rate, float source_sample_rate, uint8_t channels, uint16_t number_of_taps,
                  uint16_t number_of_filters, bool use_pre_post_filter, bool subsample_interpolate);

  /// @brief Resamples the input samples to the initalized sample rate
  /// @param input Pointer to int16_t source samples
  /// @param output Pointer to write resamples int16_t samples.
  /// @param input_frames_available Frames available at the input source pointer
  /// @param output_frames_free Frames free at the output sink pointer
  /// @param frames_used size_t passed-by-reference variable that will store the number of frames processed from the
  /// input source
  /// @param frames_generated size_t passed-by-reference variable that will store the number of output frames generated
  void resample(const int16_t *input, int16_t *output, size_t input_frames_available, size_t output_frames_free,
                size_t &frames_used, size_t &frames_generated);

 protected:
  float *float_input_buffer_{nullptr};
  size_t input_buffer_samples_;

  float *float_output_buffer_{nullptr};
  size_t output_buffer_samples_;

  Resample *resampler_{nullptr};

  Biquad lowpass_[2][2];
  BiquadCoefficients lowpass_coeff_;

  uint16_t number_of_taps_;
  uint16_t number_of_filters_;

  float sample_ratio_{1.0};
  float lowpass_ratio_{1.0};

  bool pre_filter_{false};
  bool post_filter_{false};

  uint8_t channels_;
};
}  // namespace resampler
