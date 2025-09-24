#pragma once

#if defined(__GNUC__) || defined(__clang__)
#define MACRODR_DEPRECATED(msg) __attribute__((deprecated(msg)))
#elif defined(_MSC_VER)
#define MACRODR_DEPRECATED(msg) __declspec(deprecated(msg))
#else
#define MACRODR_DEPRECATED(msg)
#endif

