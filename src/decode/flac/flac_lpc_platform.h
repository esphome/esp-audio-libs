#ifndef _FLAC_LPC_PLATFORM_H_
#define _FLAC_LPC_PLATFORM_H_

#ifdef ESP_PLATFORM
#include "sdkconfig.h"
#endif

#ifdef __XTENSA__
#include <xtensa/config/core-isa.h>
#include <xtensa/config/core-matmap.h>

// Enable assembly optimization if we have the required Xtensa features
// The assembly uses MULL (multiply) and SSR (shift amount register) instructions
#if ((XCHAL_HAVE_LOOPS == 1) && (XCHAL_HAVE_MUL32 == 1))
#define flac_lpc_asm_enabled 1
#endif

#endif // __XTENSA__

#endif // _FLAC_LPC_PLATFORM_H_