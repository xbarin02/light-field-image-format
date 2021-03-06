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

/**
 * @brief Function which rotates the input block clockwise.
 * @param table A block which shall be rotated.
 */
template <size_t BS, size_t D>
void rotate(size_t table[constpow(BS, D)]) {
  Block<size_t, BS, D> tmp {};

  for (size_t y = 0; y < BS; y++) {
    for (size_t x = 0; x < constpow(BS, D-1); x++) {
      tmp[x * BS + y] = table[y * constpow(BS, D - 1) + x];
    }
  }

  for (size_t i = 0; i < constpow(BS, D); i++) {
    table[i] = tmp[i];
  }
}

/**
 * @brief Structure which wraps the function for easy partial specialization.
 */
template <size_t BS, size_t D>
struct zigzagCore {

  /**
   * @brief Core function for zig-zag matrix generation. The function performs set of less dimensional zig-zag traversals in combination with rotations.
   * @param put   A callback function with internal counter which takes index in a block and puts traversal index in that place.
   * @param dims  Current position of a traversal pointer. Every call has its own copy.
   * @param table Pointer to (part of) traversal table for rotations.
   * @param depth The call depth which is important for rotating whole block and not only part of it.
   */
  template <typename F>
  zigzagCore(F &put, std::array<size_t, D> dims, size_t *table, size_t depth) {
    while ((dims[D-1] < BS) && (std::accumulate(dims.begin(), dims.end() - 1, 0) >= 0)) {
      std::array<size_t, D - 1> dims_new {};
      for (size_t i = 0; i < D - 1; i++) {
        dims_new[i] = dims[i];
      }

      auto put_new = [&](size_t i) {
        put(dims[D-1] * constpow(BS, D-1) + i);
      };

      zigzagCore<BS, D-1>(put_new, dims_new, table, depth + 1);

      for (int64_t i = D - 2; i >= 0; i--) {
        if ((dims[i] > 0) || (i == 0)) {
          dims[i]--;
          break;
        }
      }

      dims[D-1]++;
    }

    for (size_t i = 0; i < constpow(BS, depth); i++) {
      rotate<BS, D>(&table[i * constpow(BS, D)]);
    }
  }
};

/**
 * @brief The partial specialization for one-dimensional zig-zag.
 * @see ZigzagCore<BS,D>
 */
template <size_t BS>
struct zigzagCore<BS, 1> {

  /**
   * @brief The partial specialization for one-dimensional zig-zag.
   * @see ZigzagCore<BS,D>::ZigzagCore
   */
  template <typename F>
  zigzagCore(F &put, std::array<size_t, 1> dims, size_t *, size_t) {
    put(dims[0]);
  }
};

/**
 * @brief Function which generates and returns the matrix.
 * @return The zig-zag matrix.
 */
template <size_t BS, size_t D>
Block<size_t, BS, D> generateZigzagTable() {
  Block<size_t, BS, D> table {};

  size_t index = 0;

  std::array<size_t, D> dims {};

  while (std::accumulate(dims.begin(), dims.end(), size_t{} ) <= ((BS - 1) * D)) {
    auto put = [&](size_t i) {
      table[i] = index++;
    };

    zigzagCore<BS, D>(put, dims, table.data(), 0);

    for (size_t i = 0; i < D; i++) {
      if ((dims[i] < BS - 1) || (i == D - 1)) {
        dims[i]++;
        break;
      }
    }
  }

  return table;
}


#endif
