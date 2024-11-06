#include "resampler.h"

namespace resampler {

Resampler::~Resampler() {
  if (this->resampler_ != nullptr) {
    resampleFree(this->resampler_);
  }

  if (this->float_input_buffer_ != nullptr) {
    free(this->float_input_buffer_);
  }
  if (this->float_output_buffer_ != nullptr) {
    free(this->float_output_buffer_);
  }
};

bool Resampler::initialize(float target_sample_rate, float source_sample_rate, uint8_t channels,
                           uint16_t number_of_taps, uint16_t number_of_filters, bool use_pre_post_filter,
                           bool subsample_interpolate) {
  this->channels_ = channels;
  this->number_of_taps_ = number_of_taps;
  this->number_of_filters_ = number_of_filters;

  this->float_input_buffer_ =
      (float *) heap_caps_malloc(this->input_buffer_samples_ * sizeof(float), MALLOC_CAP_SPIRAM | MALLOC_CAP_8BIT);

  this->float_output_buffer_ =
      (float *) heap_caps_malloc(this->output_buffer_samples_ * sizeof(float), MALLOC_CAP_SPIRAM | MALLOC_CAP_8BIT);

  if ((this->float_input_buffer_ == nullptr) || (this->float_output_buffer_ == nullptr)) {
    return false;
  }

  int flags = 0;

  if (subsample_interpolate) {
    flags |= SUBSAMPLE_INTERPOLATE;
  }

  this->sample_ratio_ = target_sample_rate / source_sample_rate;

  if (this->sample_ratio_ < 1.0) {
    this->lowpass_ratio_ -= (10.24 / this->number_of_taps_);

    if (this->lowpass_ratio_ < 0.84) {
      this->lowpass_ratio_ = 0.84;
    }

    if (this->lowpass_ratio_ < this->sample_ratio_) {
      // avoid discontinuities near unity sample ratios
      this->lowpass_ratio_ = this->sample_ratio_;
    }
  }
  if (this->lowpass_ratio_ * this->sample_ratio_ < 0.98 && use_pre_post_filter) {
    float cutoff = this->lowpass_ratio_ * this->sample_ratio_ / 2.0;
    biquad_lowpass(&this->lowpass_coeff_, cutoff);
    this->pre_filter_ = true;
  }

  if (this->lowpass_ratio_ / this->sample_ratio_ < 0.98 && use_pre_post_filter && !this->pre_filter_) {
    float cutoff = this->lowpass_ratio_ / this->sample_ratio_ / 2.0;
    biquad_lowpass(&this->lowpass_coeff_, cutoff);
    this->post_filter_ = true;
  }

  if (this->pre_filter_ || this->post_filter_) {
    for (int i = 0; i < this->channels_; ++i) {
      biquad_init(&this->lowpass_[i][0], &this->lowpass_coeff_, 1.0);
      biquad_init(&this->lowpass_[i][1], &this->lowpass_coeff_, 1.0);
    }
  }

  if (this->sample_ratio_ < 1.0) {
    this->resampler_ = resampleInit(this->channels_, this->number_of_taps_, this->number_of_filters_,
                                    this->sample_ratio_ * this->lowpass_ratio_, flags | INCLUDE_LOWPASS);
  } else if (this->lowpass_ratio_ < 1.0) {
    this->resampler_ = resampleInit(this->channels_, this->number_of_taps_, this->number_of_filters_,
                                    this->lowpass_ratio_, flags | INCLUDE_LOWPASS);
  } else {
    this->resampler_ = resampleInit(this->channels_, this->number_of_taps_, this->number_of_filters_, 1.0, flags);
  }

  resampleAdvancePosition(this->resampler_, this->number_of_taps_ / 2.0);

  return true;
}

void Resampler::resample(const int16_t *input, int16_t *output, size_t input_frames_available,
                         size_t output_frames_free, size_t &frames_used, size_t &frames_generated) {
  unsigned int necessary_frames = resampleGetRequiredSamples(this->resampler_, output_frames_free, this->sample_ratio_);

  unsigned int frames_to_process = std::min(input_frames_available, necessary_frames);

  for (unsigned int i = 0; i < frames_to_process * this->channels_; ++i) {
    this->float_input_buffer_[i] = static_cast<float>(input[i]) / 32768.0f;
  }

  if (this->pre_filter_) {
    for (int i = 0; i < this->channels_; ++i) {
      biquad_apply_buffer(&this->lowpass_[i][0], this->float_input_buffer_ + i, frames_to_process, this->channels_);
      biquad_apply_buffer(&this->lowpass_[i][1], this->float_input_buffer_ + i, frames_to_process, this->channels_);
    }
  }

  ResampleResult res = resampleProcessInterleaved(this->resampler_, this->float_input_buffer_, frames_to_process,
                                                  this->float_output_buffer_, output_frames_free, this->sample_ratio_);

  frames_generated = res.output_generated;

  if (this->post_filter_) {
    for (int i = 0; i < this->channels_; ++i) {
      biquad_apply_buffer(&this->lowpass_[i][0], this->float_output_buffer_ + i, frames_generated, this->channels_);
      biquad_apply_buffer(&this->lowpass_[i][1], this->float_output_buffer_ + i, frames_generated, this->channels_);
    }
  }

  frames_used = res.input_used;

  const size_t samples_generated = frames_generated * this->channels_;
  for (size_t i = 0; i < samples_generated; ++i) {
    const float sample = this->float_output_buffer_[i];
    if (sample >= 1.0f) {
      output[i] = 32767;
    } else if (sample < -1.0f) {
      output[i] = -32768;
    } else {
      output[i] = static_cast<int16_t>(sample * 32768);
    }
  }
}

}  // namespace resampler