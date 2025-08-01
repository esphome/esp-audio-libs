#include "memory_utils.h"
#include <stdlib.h>

#ifdef ESP_PLATFORM
#include <esp_heap_caps.h>
#endif

namespace esp_audio_libs {
namespace internal {

void* alloc_psram_fallback(size_t size) {
#ifdef ESP_PLATFORM
    // Try to allocate from PSRAM first
    void* ptr = heap_caps_malloc(size, MALLOC_CAP_SPIRAM | MALLOC_CAP_8BIT);
    if (ptr == nullptr) {
        // Fall back to internal memory
        ptr = heap_caps_malloc(size, MALLOC_CAP_8BIT);
    }
    return ptr;
#else
    // For non-ESP platforms, use standard malloc
    return malloc(size);
#endif
}

void free_psram_fallback(void* ptr) {
#ifdef ESP_PLATFORM
    heap_caps_free(ptr);
#else
    free(ptr);
#endif
}

}  // namespace internal
}  // namespace esp_audio_libs