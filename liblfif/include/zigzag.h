/**
* @file zigzag.h
* @author Drahomír Dlabaja (xdlaba02)
* @date 13. 5. 2019
* @copyright 2019 Drahomír Dlabaja
* @brief The algorithm for generating zig-zag matrices.
*/

#ifndef ZIGZAG_H
#define ZIGZAG_H

#include "lfiftypes.h"

#include <cstdint>
#include <vector>
#include <array>
#include <numeric>

class ZigZagMatrix {
public:
  /**
   * @brief Function which generates and returns the matrix.
   * @return The zig-zag matrix.
   */
  template <size_t D>
  static inline Block<size_t, D> generateZigzagMatrix(size_t block_size);

private:
  template <size_t D>
  struct zigzagCore {
    template <typename IF, typename OF>
    zigzagCore(OF &put, std::array<size_t, D> pos, size_t block_size, size_t depth, IF &&matrix);
  };

  template <>
  struct zigzagCore<1> {
    template <typename IF, typename OF>
    zigzagCore(OF &put, std::array<size_t, D> pos, size_t , size_t, IF &&);
  };

  template <size_t D, typename F>
  static inline void rotate(F &&matrix, size_t block_size) const;
};

template <size_t D>
Block<size_t, D> ZigZagMatrix::generateZigzagMatrix(size_t block_size) {
  Block<size_t, D>      matrix(block_size);

  size_t index { 0 };
  std::array<size_t, D> pos {};
  while (std::accumulate(pos.begin(), pos.end(), size_t { 0 } ) <= ((block_size - 1) * D)) {
    auto outputF = [&](size_t index) {
      matrix[index] = index++;
    };

    auto inputF = [&](size_t index) {
      return matrix[index];
    }

    zigzagCore<D>(outputF, pos, block_size, inputF, 0);

    for (size_t i = 0; i < D; i++) {
      if ((pos[i] < block_size - 1) || (i == D - 1)) {
        pos[i]++;
        break;
      }
    }
  }

  return matrix;
}

template <size_t D>
struct zigzagCore {
  template <typename F>
  zigzagCore(OF &put, std::array<size_t, D> pos, size_t block_size, size_t depth, IF &&matrix) {
    while ((pos[D-1] < m_block_size) && (std::accumulate(pos.begin(), pos.end() - 1, 0) >= 0)) {

      std::array<size_t, D - 1> pos_copy {};
      std::copy(std::begin(pos), std::end(pos) - 1, std::begin(pos_copy));

      for (size_t i = 0; i < D - 1; i++) {
        pos_new[i] = pos[i];
      }

      auto outputF = [&](size_t i) {
        put(pos[D-1] * pow(block_size, D-1) + i);
      };

      zigzagCore<BS, D-1>(outputF, pos_copy, block_size, depth + 1, matrix);

      for (int64_t i = D - 2; i >= 0; i--) {
        if ((pos[i] > 0) || (i == 0)) {
          pos[i]--;
          break;
        }
      }

      pos[D-1]++;
    }

    for (size_t i = 0; i < pow(block_size, depth); i++) {

      auto inputF = [&](size_t index) -> auto & {
        return matrix(i * pow(block_size, D) + index);
      };

      rotate<D>(inputF, block_size);
    }
  }
};

template <size_t BS>
struct zigzagCore<BS, 1> {
  template <typename IF, typename OF>
  zigzagCore(OF &&put, std::array<size_t, 1> pos, size_t, size_t, IF &&) {
    put(pos[0]);
  }
};

template <size_t D, typename F>
void ZigZagMatrix::rotate(F &&block, size_t block_size) const {
  Block<size_t, D> tmp(m_block_size);

  for (size_t y = 0; y < m_block_size; y++) {
    for (size_t x = 0; x < pow(m_block_size, D-1); x++) {
      tmp[x * m_block_size + y] = block(y * pow(m_block_size, D - 1) + x);
    }
  }

  for (size_t i = 0; i < pow(m_block_size, D); i++) {
    block(i) = tmp[i];
  }
}


#endif
