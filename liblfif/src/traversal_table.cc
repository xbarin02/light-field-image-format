#include "traversal_table.h"

#include <algorithm>

template class TraversalTable<2>;
template class TraversalTable<3>;
template class TraversalTable<4>;

template <size_t D>
TraversalTable<D> &TraversalTable<D>::constructByReference(const ReferenceBlock<D> &reference_block) {
  Block<std::pair<double, size_t>, D> srt {};

  for (size_t i = 0; i < constpow(8, D); i++) {
    srt[i].first += abs(reference_block[i]);
    srt[i].second = i;
  }

  stable_sort(srt.rbegin(), srt.rend());

  for (size_t i = 0; i < constpow(8, D); i++) {
    m_block[srt[i].second] = i;
  }

  return *this;
}

template <size_t D>
TraversalTable<D> &TraversalTable<D>::constructByRadius() {
  Block<std::pair<double, size_t>, D> srt {};

  for (size_t i = 0; i < constpow(8, D); i++) {
    for (size_t j = 1; j <= D; j++) {
      size_t coord = (i % constpow(8, j)) / constpow(8, j-1);
      srt[i].first += coord * coord;
    }
    srt[i].second = i;
  }

  stable_sort(srt.begin(), srt.end());

  for (size_t i = 0; i < constpow(8, D); i++) {
    m_block[srt[i].second] = i;
  }

  return *this;
}

template <size_t D>
TraversalTable<D> &TraversalTable<D>::constructByDiagonals() {
  Block<std::pair<double, size_t>, D> srt {};

  for (size_t i = 0; i < constpow(8, D); i++) {
    for (size_t j = 1; j <= D; j++) {
      srt[i].first += (i % constpow(8, j)) / constpow(8, j-1);
    }
    srt[i].second = i;
  }

  stable_sort(srt.begin(), srt.end());

  for (size_t i = 0; i < constpow(8, D); i++) {
    m_block[srt[i].second] = i;
  }

  return *this;
}

template <size_t D>
TraversalTable<D> &TraversalTable<D>::writeToStream(std::ofstream &stream) {
  writeToStreamHelper<TraversalTableUnit>(stream, m_block.data(), constpow(8, D));
  return *this;
}

template <size_t D>
TraversalTable<D> &TraversalTable<D>::readFromStream(std::ifstream &stream) {
  readFromStreamHelper<TraversalTableUnit>(stream, m_block.data(), constpow(8, D));
  return *this;
}

template <size_t D>
TraversalTableUnit TraversalTable<D>::operator [](size_t index) const {
  return m_block[index];
}
