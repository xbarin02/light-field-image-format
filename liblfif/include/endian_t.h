/**
* @file endian_t.h
* @author Drahomír Dlabaja (xdlaba02)
* @date 12. 5. 2019
* @copyright 2019 Drahomír Dlabaja
* @brief Functions which performs FDCT and IDCT.
*/

#include <algorithm>

#ifndef ENDIAN_T_H
#define ENDIAN_T_H

/**
* @brief  Function which performs endian swap.
* @param  data The data to be swapped.
* @return Swapped data.
*/
template <typename T>
inline T endianSwap(T data) {
  uint8_t *ptr = reinterpret_cast<uint8_t *>(&data);
  std::reverse(ptr, ptr + sizeof(T));
  return data;
}

/**
* @brief  Function which performs swap only if local endian != big endian.
* @param  data The data to be swapped.
* @return Swapped data.
*/
template <typename T>
inline T bigEndianSwap(T data) {
  if (htobe16(1) == 1) {
    return data;
  }
  else {
    return endianSwap(data);
  }
}

/**
* @brief Function which writes value to stream in big endian.
* @param data The data to be written.
* @param stream The stream to which the data will be written.
*/
template<typename T>
inline void writeValueToStream(T data, std::ostream &stream) {
  T dataBE = bigEndianSwap(data);
  stream.write(reinterpret_cast<const char *>(&dataBE), sizeof(dataBE));
}

/**
* @brief  Function which reads big endian value from stream.
* @param  stream The stream from which the data will be read.
* @return The data to be read.
*/
template<typename T>
inline T readValueFromStream(std::istream &stream) {
  T dataBE {};
  stream.read(reinterpret_cast<char *>(&dataBE), sizeof(dataBE));
  return bigEndianSwap(dataBE);
}

#endif
