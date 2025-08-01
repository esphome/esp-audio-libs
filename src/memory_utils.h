#pragma once

#include <stddef.h>

namespace esp_audio_libs {
namespace internal {

/**
 * @brief Allocate memory preferring PSRAM, falling back to internal memory
 * 
 * This function first attempts to allocate memory from PSRAM. If that fails,
 * it falls back to allocating from internal memory.
 * 
 * @param size Size in bytes to allocate
 * @return Pointer to allocated memory, or nullptr if allocation fails
 */
void* alloc_psram_fallback(size_t size);

/**
 * @brief Free memory allocated by alloc_psram_fallback
 * 
 * @param ptr Pointer to memory to free
 */
void free_psram_fallback(void* ptr);

}  // namespace internal
}  // namespace esp_audio_libs