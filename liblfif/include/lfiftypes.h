/**
* @file lfiftypes.h
* @author Drahomír Dlabaja (xdlaba02)
* @date 12. 5. 2019
* @copyright 2019 Drahomír Dlabaja
* @brief Some basic types used by the library.
*/

#ifndef LFIFTYPES_H
#define LFIFTYPES_H

#include "constpow.h"
#include "stack_allocator.h"

#include <cstdint>

#include <array>
#include <fstream>

using RGBUNIT   = uint16_t; /**< @brief Unit which is intended to containt RGB data.*/
using QDATAUNIT = int64_t;  /**< @brief Unit which is intended to containt quantized DCT coefficients.*/
using INPUTUNIT = float;    /**< @brief Unit which is intended to containt input/output YCbCr/YCoCg/other data.*/

using INPUTTRIPLET = std::array<INPUTUNIT, 3>; /**< @brief Unit which is intended to containt one triplet of input data.*/

template<typename T, size_t D>
class Block {
  Block(size_t block_size): m_block_size { block_size }, m_data {} {
    m_data = StackAllocator::malloc((pow(block_size, D) * sizeof(T));
  }

  Block(const Block &block): Block(block.size()) {
    operator =(block);
  }

  Block(const Block &&block) m_block_size { block.size() }, m_data { block.m_data } {
    block.m_block_size = 0;
    block.m_data = nullptr;
  }

  ~Block() {
    StackAllocator::free((pow(m_block_size, D) * sizeof(T));
  }

  T &operator [](size_t index) {
    return m_data[index];
  }

  T operator [](size_t index) const {
    return m_data[index];
  }

  Block &operator =(const Block &rhs) {
    assert(m_block_size == rhs.size());

    for (size_t i = 0; i < pow(m_block_size, D); i++) {
      m_data[i] = rhs[i];
    }
  }

  size_t size() {
    return m_block_size;
  }

private:
  size_t m_block_size;
  T     *m_data;
}

#endif
