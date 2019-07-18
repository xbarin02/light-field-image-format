/**
* @file traversal_table.h
* @author Drahomír Dlabaja (xdlaba02)
* @date 13. 5. 2019
* @copyright 2019 Drahomír Dlabaja
* @brief Module for generating linearization matrices.
*/

#ifndef TRAVERSAL_TABLE_H
#define TRAVERSAL_TABLE_H

#include "quant_table.h"
#include "zigzag.h"
#include "endian_t.h"

#include <iosfwd>

using REFBLOCKUNIT = double; /**< @brief Type intended to be used in reference block.*/
using TTABLEUNIT = size_t;   /**< @brief Type intended to be used in traversal matrix.*/

/**
 * @brief Reference block type.
 */
template<size_t D>
using ReferenceBlock = Block<REFBLOCKUNIT, D>;

/**
 * @brief Traversal table class.
 */
template <size_t D>
class TraversalTable {
public:

  /**
   * @brief Method which constructs traversal matrix by reference block.
   * @param reference_block The block from which the matrix shall be constructed.
   */
  void constructByReference(const ReferenceBlock<D> &reference_block);

  /**
   * @brief Method which constructs traversal matrix from quantization matrix.
   * @param quant_table The quantization matrix from which the matrix shall be constructed.
   */
  void constructByQuantTable(const QuantTable<D> &quant_table);

  /**
   * @brief Method which constructs traversal matrix by Eucleidian distance from the DC coefficient.
   */
  void constructByRadius();

  /**
   * @brief Method which constructs traversal matrix by manhattan distance from the DC coefficient.
   */
  void constructByDiagonals();

  /**
   * @brief Method which constructs traversal matrix by Eucleidian distance from the AC coefficient with highest frequency.
   */
  void constructByHyperboloid();

  /**
   * @brief Method which constructs zig-zag matrix.
   */
  void constructZigzag();

  /**
   * @brief Method which writes the matrix to stream.
   * @param stream The stream to which the table shall be written.
   */
  void writeToStream(std::ostream &stream) const;

  /**
   * @brief Method which read the matrix from stream.
   * @param stream The stream from which the table shall be read.
   */
  void readFromStream(std::istream &stream);

  /**
   * @brief Method used to index the internal matrix.
   * @param index The index of a linearization index.
   * @return Linearization index.
   */
  TTABLEUNIT operator [](size_t index) const {
    return m_block[index];
  }

private:
  Block<TTABLEUNIT, D> m_block;
};

template <size_t D>
void TraversalTable<D>::constructByReference(const ReferenceBlock<D> &reference_block) {
  Block<std::pair<REFBLOCKUNIT, size_t>, D> srt(reference_block.size());

  for (size_t i = 0; i < pow(reference_block.size(), D); i++) {
    srt[i].first += reference_block[i];
    srt[i].second = i;
  }

  //do NOT sort DC coefficient, thus +1 at the begining.
  stable_sort(srt.begin() + 1, srt.end(), [](auto &left, auto &right) {
    return left.first > right.first;
  });

  m_block.resize(reference_block.size());

  for (size_t i = 0; i < pow(reference_block.size(), D); i++) {
    m_block[srt[i].second] = i;
  }
}

template <size_t D>
void TraversalTable<D>::constructByQuantTable(const QuantTable<BS, D> &quant_table) {
  ReferenceBlock<D> dummy(quant_table.size());

  for (size_t i = 0; i < pow(reference_block.size(), D); i++) {
    dummy[i] = -quant_table[i];
  }

  constructByReference(dummy);
}

template <size_t D>
void TraversalTable<D>::constructByRadius(size_t block_size) {
  ReferenceBlock<D> dummy(block_size);

  for (size_t i = 0; i < pow(BS, D); i++) {
    for (size_t j = 1; j <= D; j++) {
      size_t coord = (i % pow(BS, j)) / pow(BS, j - 1);
      dummy[i] -= coord * coord;
    }
  }

  constructByReference(dummy);
}

template <size_t D>
void TraversalTable<D>::constructByDiagonals(size_t block_size) {
  ReferenceBlock<D> dummy(block_size);

  for (size_t i = 0; i < pow(block_size, D); i++) {
    for (size_t j = 1; j <= D; j++) {
      dummy[i] -= (i % pow(block_size, j)) / pow(block_size, j - 1);
    }
  }

  constructByReference(dummy);
}

template <size_t D>
void TraversalTable<D>::constructByHyperboloid(size_t block_size) {
  ReferenceBlock<D> dummy(block_size);

  dummy.fill(1);

  for (size_t i = 0; i < pow(block_size, D); i++) {
    for (size_t j = 1; j <= D; j++) {
      size_t coord = (i % pow(block_size, j)) / pow(block_size, j - 1);
      dummy[i] *= (coord + 1);
    }
    dummy[i] *= -1;
  }

  constructByReference(dummy);
}

template <size_t D>
void TraversalTable<D>::constructZigzag(size_t block_size) {
  m_block = ZigZagMatrix::generateZigzagMatrix<D>(block_size);
}

template <size_t D>
void TraversalTable<D>::writeToStream(std::ostream &stream) const {
  size_t max_bits = ceil(log2(pow(m_block.size(), D)));

  if (max_bits <= 8) {
    for (size_t i = 0; i < pow(m_block.size(), D); i++) {
      writeValueToStream<uint8_t>(m_block[i], stream);
    }
  }
  else if (max_bits <= 16) {
    for (size_t i = 0; i < pow(m_block.size(), D); i++) {
      writeValueToStream<uint16_t>(m_block[i], stream);
    }
  }
  else if (max_bits <= 32) {
    for (size_t i = 0; i < pow(m_block.size(), D); i++) {
      writeValueToStream<uint32_t>(m_block[i], stream);
    }
  }
  else if (max_bits <= 64) {
    for (size_t i = 0; i < pow(m_block.size(), D); i++) {
      writeValueToStream<uint64_t>(m_block[i], stream);
    }
  }
}


template <size_t BS, size_t D>
void TraversalTable<BS, D>::readFromStream(std::istream &stream) {
  size_t max_bits = ceil(log2(pow(BS, D)));

  if (max_bits <= 8) {
    for (size_t i = 0; i < pow(BS, D); i++) {
      m_block[i] = readValueFromStream<uint8_t>(stream);
    }
  }
  else if (max_bits <= 16) {
    for (size_t i = 0; i < pow(BS, D); i++) {
      m_block[i] = readValueFromStream<uint16_t>(stream);
    }
  }
  else if (max_bits <= 32) {
    for (size_t i = 0; i < pow(BS, D); i++) {
      m_block[i] = readValueFromStream<uint32_t>(stream);
    }
  }
  else if (max_bits <= 64) {
    for (size_t i = 0; i < pow(BS, D); i++) {
      m_block[i] = readValueFromStream<uint64_t>(stream);
    }
  }
}

#endif
