# esp-audio-libs

A collection of libraries and functions that are useful for playing audio on ESP32 devices. It includes code based on the following:
- [esp-dsp](https://github.com/espressif/esp-dsp) assembly functions for floating point dot products and Q15 fixed point addition and constant multipliction.
    - Author: Espressif
    - License: Apache v2.0
- [ART-resampler](https://github.com/dbry/audio-resampler) for resampling audio, optimized with assembly dot product functions.
    - Author: David Bryant
    - License: BSD-3-Clause
- [wav-decoder](https://github.com/synesthesiam/wav-decoder) for parsing WAV file headers.
    - Author: Michael Hansen
    - License: MIT
- [flac-decoder](https://github.com/synesthesiam/flac-decoder) for decoding FLAC files.
    - Author: Michael Hansen
    - License: ESPHome License
    - Port of [Simple FLAC Decoder (Python)](https://www.nayuki.io/res/simple-flac-implementation/simple-decode-flac-to-wav.py) by Project Nayuki with MIT License
- [libHelix MP3 Decoder](https://en.wikipedia.org/wiki/Helix_Universal_Server) for decoding MP3 files
    - Author: Real Networks
    - License: RCSL