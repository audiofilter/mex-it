// We just need the mex_function header file
// This could be put into mex-it.h but this way allows naming of the mex function in matlab to this file name (without the .cxx)
#include <cstdint>
#include <iostream>
void mex_function(double d,
                  float f,
                  int8_t i8,
                  int16_t i16,
                  int32_t i32,
                  int64_t i64,
                  uint8_t u8,
                  uint16_t u16,
                  uint32_t u32,
                  uint64_t u64,
                  bool b,
                  double& sum) {
  if (b) {
    sum = d + f + double(i8 + i16 + i32 + i64) + double(u8 + u16 + u32 + u64);
  } else {
    sum = 0;
  }

  std::cout << "inputs validated ok and are = "
            << "(double) " << d << ","
            << "(float) " << f << ","
            << "(int8_t) " << (int)i8 << ","
            << "(int16_t) " << i16 << ","
            << "(int32_t) " << i32 << ","
            << "(int64_t) " << i64 << ","
            << "(uint8_t) " << (int)u8 << ","
            << "(uint16_t) " << u16 << ","
            << "(uint32_t) " << u32 << ","
            << "(uint64_t) " << u64 << ","
            << "(bool) " << b << "\n"
            << "result = " << sum << "\n";
}

#include "mex-it.h"
