/**
* @file block.h
* @author Drahomír Dlabaja (xdlaba02)
* @date 12. 5. 2019
* @copyright 2019 Drahomír Dlabaja
* @brief Functions for block extraction and insertion.
*/

#ifndef BLOCK_H
#define BLOCK_H

#include "lfiftypes.h"

#include <cmath>

/**
* @brief Struct for block extraction which wraps static parameters for partial specialization.
*/
template<size_t D>
struct getBlock {

  /**
   * @brief Function which extracts block.
   * @param input  The input callback function with signature T input(size_t index), where T is type of extracted sample.
   * @param block  Index of a wanted block in memory.
   * @param dims   The dimensions of an extracted volume.
   * @param output The output callback function with signature void output(size_t index, T value), where T is type of extracted sample.
   */
  template <typename IF, typename OF>
  getBlock(IF &&input, size_t index, const size_t dims[D], OF &&output, size_t block_size) {
    size_t blocks_x = 1;
    size_t size_x   = 1;

    for (size_t i = 0; i < D - 1; i++) {
      blocks_x *= ceil(dims[i] / static_cast<double>(block_size));
      size_x *= dims[i];
    }

    for (size_t pixel = 0; pixel < block_size; pixel++) {
      size_t image_y = (index / blocks_x) * block_size + pixel;

      if (image_y >= dims[D-1]) {
        image_y = dims[D-1] - 1;
      }

      auto inputF = [&](size_t image_index) {
        return input(image_y * size_x + image_index);
      };

      auto outputF = [&](size_t pixel_index, const auto &value) {
        output(pixel * pow(block_size, D-1) + pixel_index, value);
      };

      getBlock<D-1>(inputF, index % blocks_x, dims, outputF, block_size);
    }
  }
};

/**
 * @brief The parital specialization for getting one dimensional block.
 * @see getBlock<block_size, D>
 */
template<>
struct getBlock<1> {

  /**
   * @brief The parital specialization for getting one dimensional block.
   * @see getBlock<block_size, D>::getBlock
   */
  template <typename IF, typename OF>
  getBlock(IF &&input, const size_t index, const size_t dims[1], OF &&output, size_t block_size) {
    for (size_t pixel = 0; pixel < block_size; pixel++) {
      size_t image = index * block_size + pixel;

      if (image >= dims[0]) {
        image = dims[0] - 1;
      }

      output(pixel, input(image));
    }
  }
};

/**
* @brief Struct for block insertion which wraps static parameters for partial specialization.
*/
template<size_t D>
struct putBlock {

  /**
   * @brief Function which inserts block.
   * @param input  The input callback function with signature T input(size_t index), where T is type of inserted sample.
   * @param block  Index of a put block block in memory.
   * @param dims   The dimensions of an inserted volume to which block will be inserted.
   * @param output The output callback function with signature void output(size_t index, T value), where T is type of inserted sample.
   */
  template <typename IF, typename OF>
  putBlock(IF &&input, size_t index, const size_t dims[D], OF &&output, size_t block_size) {
    size_t blocks_x = 1;
    size_t size_x   = 1;

    for (size_t i = 0; i < D - 1; i++) {
      blocks_x *= ceil(dims[i] / static_cast<double>(block_size));
      size_x *= dims[i];
    }

    for (size_t pixel = 0; pixel < block_size; pixel++) {
      size_t image = (index / blocks_x) * block_size + pixel;

      if (image >= dims[D-1]) {
        break;
      }

      auto inputF = [&](size_t pixel_index) {
        return input(pixel * pow(block_size, D-1) + pixel_index);
      };

      auto outputF = [&](size_t image_index, const auto &value) {
        return output(image * size_x + image_index, value);
      };

      putBlock<D-1>(inputF, index % blocks_x, dims, outputF, block_size);
    }
  }
};

/**
 * @brief The parital specialization for putting one dimensional block.
 */
template<>
struct putBlock<1> {

  /**
   * @brief The parital specialization for putting one dimensional block.
   */
  template <typename IF, typename OF>
  putBlock(IF &&input, size_t index, const size_t dims[1], OF &&output, size_t block_size) {
    for (size_t pixel = 0; pixel < block_size; pixel++) {
      size_t image = index * block_size + pixel;

      if (image >= dims[0]) {
        break;
      }

      output(image, input(pixel));
    }
  }
};

#endif
