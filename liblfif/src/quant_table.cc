#include "quant_table.h"

#include <algorithm>

template class QuantTable<2>;
template class QuantTable<3>;
template class QuantTable<4>;

using namespace std;

template <size_t D>
QuantTable<D> &
QuantTable<D>::baseDiagonalTable(const std::array<QTABLEUNIT, 64> &qtable, uint8_t quality) {
  for (size_t i = 0; i < constpow(BLOCK_SIZE, D); i++) {
    size_t diagonal = 0;
    for (size_t j = 0; j < D; j++) {
      diagonal += (i % constpow(BLOCK_SIZE, j + 1)) / constpow(BLOCK_SIZE, j);
    }
    m_block[i] = averageForDiagonal(qtable, diagonal);
  }

  scaleByQuality(quality);

  return *this;
}

template <size_t D>
QuantTable<D> &
QuantTable<D>::baseCopyTable(const std::array<QTABLEUNIT, 64> &qtable, uint8_t quality) {
  for (size_t slice = 0; slice < constpow(BLOCK_SIZE, D - 2); slice++) {
    for (size_t y = 0; y < BLOCK_SIZE; y++) {
      size_t orig_y = y < 8 ? y : 7;

      for (size_t x = 0; x < BLOCK_SIZE; x++) {
        size_t orig_x = x < 8 ? x : 7;

        m_block[(slice * BLOCK_SIZE + y) * BLOCK_SIZE + x] = qtable[orig_y * 8 + orig_x];
      }
    }
  }

  scaleByQuality(quality);

  return *this;
}

template <size_t D>
QuantTable<D> &
QuantTable<D>::baseUniform(uint8_t quality) {
  for (size_t i = 0; i < constpow(BLOCK_SIZE, D); i++) {
    m_block[i] = 101 - quality;
  }

  return *this;
}

template <size_t D>
QuantTable<D> &
QuantTable<D>::writeToStream(std::ofstream &stream) {
  for (size_t i = 0; i < constpow(BLOCK_SIZE, D); i++) {
    writeValueToStream<uint8_t>(stream, m_block[i]);
  }

  return *this;
}

template <size_t D>
QuantTable<D> &
QuantTable<D>::readFromStream(std::ifstream &stream) {
  for (size_t i = 0; i < constpow(BLOCK_SIZE, D); i++) {
    m_block[i] = readValueFromStream<uint8_t>(stream);
  }
  return *this;
}

template <size_t D>
uint64_t QuantTable<D>::operator [](size_t index) const {
  return m_block[index];
}

template <size_t D>
QTABLEUNIT
QuantTable<D>::averageForDiagonal(const std::array<QTABLEUNIT, 64> &qtable, size_t diagonal) {
  double sum = 0;
  size_t cnt = 0;

  for (size_t y = 0; y < 8; y++) {
    for (size_t x = 0; x < 8; x++) {
      if (x + y == diagonal) {
        sum += qtable[y*8+x];
        cnt++;
      }
    }
  }

  if (cnt) {
    return (sum / cnt) + 0.5;
  }

  return qtable[63];
}

template <size_t D>
QuantTable<D> &
QuantTable<D>::scaleByQuality(uint8_t quality) {
  double scale_coef = quality < 50 ? (50.0 / quality) : (200.0 - 2 * quality) / 100;

  for (size_t i = 0; i < constpow(BLOCK_SIZE, D); i++) {
    m_block[i] = clamp<double>(m_block[i] * scale_coef, 1, 255);
  }

  return *this;
}
