/**
* @file predict.h
* @author Drahomír Dlabaja (xdlaba02)
* @date 9. 8. 2019
* @copyright 2019 Drahomír Dlabaja
* @brief Prediction stuff.
*/

#ifndef PREDICT_H
#define PREDICT_H

#include "lfiftypes.h"

#include <cstdint>
#include <cassert>
#include <cmath>


template<size_t BS, size_t D, typename T>
void getSlice(const Block<T, BS, D> &block, Block<T, BS, D - 1> &slice, size_t direction, size_t index) {
  assert(direction < D);
  assert(index < BS);

  for (size_t i { 0 }; i < constpow(BS, D - 1); i++) {
    size_t directional_index = i % constpow(BS, direction) + i / constpow(BS, direction) * constpow(BS, direction + 1);
    slice[i] = block[directional_index + index * constpow(BS, direction)];
  }
}

template<size_t BS, size_t D, typename T>
void putSlice(Block<T, BS, D> &block, const Block<T, BS, D - 1> &slice, size_t direction, size_t index) {
  assert(direction < D);
  assert(index < BS);

  for (size_t i { 0 }; i < constpow(BS, D - 1); i++) {
    size_t directional_index = i % constpow(BS, direction) + i / constpow(BS, direction) * constpow(BS, direction + 1);
    block[directional_index + index * constpow(BS, direction)] = slice[i];
  }
}

template<size_t BS, size_t D>
void predict_DC(Block<INPUTUNIT, BS, D> &predicted, const size_t block_dims[D], const INPUTUNIT *decoded) {
  INPUTUNIT avg { 0 };

  auto multDims = [&](size_t n) {
    size_t product { 1 };
    for (size_t i { 0 }; i < n; i++) {
      product *= block_dims[i] * BS;
    }
    return product;
  };

  for (size_t dir { 0 }; dir < D; dir++) {
    for (size_t i { 0 }; i < constpow(BS, D - 1); i++) {
      const INPUTUNIT *ptr { decoded };

      size_t pow_val { 0 };
      for (size_t d {0}; d < D; d++) {
        if (d == dir) {
          ptr -= multDims(d);
        }
        else {
          ptr += multDims(d) * ((i % constpow(BS, pow_val + 1)) / constpow(BS, pow_val));
          pow_val++;
        }
      }

      avg += *ptr;
    }
  }

  avg /= constpow(BS, D - 1) * D;

  for (size_t i { 0 }; i < constpow(BS, D); i++) {
    predicted[i] = avg;
  }
}

template <size_t BS, size_t D>
void project_neighbours_to_main_ref(Block<INPUTUNIT, BS * 2 + 1, D - 1> &main_ref, const int8_t direction[D], const INPUTUNIT *src, const size_t input_stride[D + 1]) {
  size_t main_ref_idx {};
  int64_t ref_offset {};

  for (size_t d = 0; d < D; d++) {
    if (direction[d] >= direction[main_ref_idx]) {
      main_ref_idx = d;
    }
  }

  for (size_t d { 0 }; d < D; d++) {
    if (d != main_ref_idx) {
      if (direction[d] > 0) {
        for (size_t i { 0 }; i < constpow(BS * 2 + 1, D - 1); i++) {
          int64_t src_idx {};

          for (size_t dd = 0; dd < D - 1; dd++) {
            int64_t idx = dd < d ? dd : dd + 1;
            src_idx += (i % constpow(BS * 2 + 1, dd + 1)) / constpow(BS * 2 + 1, dd) * input_stride[idx];
          }

          int64_t dst_idx {};
          std::array<int64_t, D> dst_pos {};

          //rotate index to direction d
          dst_idx = i % constpow(BS * 2 + 1, d) + i / constpow(BS * 2 + 1, d) * constpow(BS * 2 + 1, d + 1);

          bool overflow {};
          for (size_t dd { 0 }; dd < D; dd++) {
            dst_pos[dd] = dst_idx % constpow(BS * 2 + 1, dd + 1) / constpow(BS * 2 + 1, dd);

            if (dd != main_ref_idx) {
              if (direction[dd] >= 0) {
                dst_pos[dd] += BS;
              }
              if (direction[dd] <= 0){
                dst_pos[dd] += 1;
              }
            }

            if (dst_pos[dd] >= static_cast<int64_t>(BS * 2 + 1)) {
              overflow = true;
            }
          }

          if (dst_pos[main_ref_idx] > static_cast<int64_t>(BS)) {
            overflow = true;
          }

          if (overflow) {
            continue;
          }

          int64_t distance { dst_pos[main_ref_idx] };

          for (size_t dd { 0 }; dd < D; dd++) {
            dst_pos[dd] *= direction[main_ref_idx];
            dst_pos[dd] -= direction[dd] * distance;
            dst_pos[dd] /= direction[main_ref_idx];
          }

          for (size_t dd { 0 }; dd < D; dd++) {
            if (dst_pos[dd] < 0) {
              overflow = true;
            }
            else if (dst_pos[dd] >= static_cast<int64_t>(BS * 2 + 1)) {
              overflow = true;
            }
          }

          if (overflow) {
            continue;
          }

          dst_idx = 0;
          size_t pow {};
          for (size_t dd { 0 }; dd < D; dd++) {
            if (dd != main_ref_idx) {
              dst_idx += dst_pos[dd] * constpow(BS * 2 + 1, pow);
              pow++;
            }
          }

          // Tady by se taky mozna hodilo kontrolovat, jestli tam uz neni nejaka hodnota s mensi chybou
          main_ref[dst_idx] = src[src_idx];
        }
      }
    }
  }

  size_t idx {};
  for (size_t d { 0 }; d < D; d++) {
    if (d != main_ref_idx) {
      if (direction[d] >= 0) {
        ref_offset += constpow(BS * 2 + 1, idx) * BS;
      }
      if (direction[d] <= 0){
        ref_offset += constpow(BS * 2 + 1, idx);
      }
      idx++;
    }
  }

  // copy main samples to main reference
  for (size_t i { 0 }; i < constpow(BS * 2 + 1, D - 1); i++) {
    int64_t rot_idx {};
    int64_t src_idx {};
    size_t dst_idx {};

    rot_idx = i % constpow(BS * 2 + 1, main_ref_idx) + i / constpow(BS * 2 + 1, main_ref_idx) * constpow(BS * 2 + 1, main_ref_idx + 1);


    for (size_t dd = 0; dd < D; dd++) {
      src_idx += (rot_idx % constpow(BS * 2 + 1, dd + 1)) / constpow(BS * 2 + 1, dd) * input_stride[dd];
    }

    dst_idx = i + ref_offset;

    bool overflow {};
    for (size_t dd = 0; dd < D; dd++) {
      if (((dst_idx % constpow(BS * 2 + 1, dd + 1)) / constpow(BS * 2 + 1, dd)) < ((ref_offset % constpow(BS * 2 + 1, dd + 1)) / constpow(BS * 2 + 1, dd))) {
        overflow = true;
      }
    }

    if (overflow) {
      continue;
    }


    main_ref[dst_idx] = src[src_idx];
  }
}

template <size_t BS, size_t D>
struct interpolate {

  template <typename F>
  interpolate(F &&main_ref, const int64_t main_ref_pos[D], int8_t multiplier, INPUTUNIT &output) {
    int64_t pos = std::floor(main_ref_pos[D - 1] / static_cast<double>(multiplier));
    int64_t frac = main_ref_pos[D - 1] % multiplier;

    auto inputF1 = [&](size_t index) {
      return main_ref(pos * constpow(BS, D - 1) + index);
    };

    if (frac == 0) {
      interpolate<BS, D - 1>(inputF1, main_ref_pos, multiplier, output);
    }
    else {
      INPUTUNIT val1 {};
      INPUTUNIT val2 {};

      auto inputF2 = [&](size_t index) {
        return main_ref((pos + 1) * constpow(BS, D - 1) + index);
      };

      interpolate<BS, D - 1>(inputF1, main_ref_pos, multiplier, val1);
      interpolate<BS, D - 1>(inputF2, main_ref_pos, multiplier, val2);

      output = (val1 * (multiplier - frac) + val2 * frac) / multiplier;
    }
  }
};

template <size_t BS>
struct interpolate<BS, 0> {
  template <typename F>
  interpolate(F &&main_ref, const int64_t *, int8_t, INPUTUNIT &output) {
    output = main_ref(0);
  }
};

template <size_t BS, size_t D>
void predict_from_main_ref(Block<INPUTUNIT, BS, D> &output, const int8_t direction[D], const Block<INPUTUNIT, BS * 2 + 1, D - 1> &main_ref) {
  size_t main_ref_idx {};

  for (size_t d = 0; d < D; d++) {
    if (direction[d] >= direction[main_ref_idx]) {
      main_ref_idx = d;
    }
  }

  for (size_t i { 0 }; i < constpow(BS, D); i++) {
    std::array<int64_t, D> pos {};
    int64_t main_ref_pos[D - 1] {};

    for (size_t d { 0 }; d < D; d++) {
      pos[d] = ((i % constpow(BS, d + 1) / constpow(BS, d)) + 1);
    }

    for (size_t d { 0 }; d < D - 1; d++) {
      size_t idx = d < main_ref_idx ? d : d + 1;
      main_ref_pos[idx] = (pos[idx] * direction[main_ref_idx]) - (direction[idx] * pos[main_ref_idx]);

      // find offset for main neighbour to make space for projected samples
      if (direction[idx] >= 0) {
        main_ref_pos[idx] += BS * direction[main_ref_idx];
      }
    }

    auto inputF = [&](size_t index) {
      return main_ref[index];
    };

    interpolate<BS * 2 + 1, D - 1>(inputF, main_ref_pos, direction[main_ref_idx], output[i]);
  }
}

template <size_t BS, size_t D>
void predict_direction(Block<INPUTUNIT, BS, D> &output, const int8_t direction[D], const INPUTUNIT *src, const size_t input_stride[D + 1]) {
  Block<INPUTUNIT, BS * 2 + 1, D - 1> ref {};

  int64_t ptr_offset { 0 };
  size_t  main_ref   { 0 };

  // find which neighbouring block will be main
  for (size_t d = 0; d < D; d++) {
    if (direction[d] >= direction[main_ref]) {
      main_ref = d;
    }
  }

  if (direction[main_ref] <= 0) {
    return;
  }

  for (size_t d { 0 }; d < D; d++) {
    if (direction[d] > 0) {
      ptr_offset -= input_stride[d];
    }
  }

  project_neighbours_to_main_ref<BS, D>(ref, direction, &src[ptr_offset], input_stride);

  predict_from_main_ref<BS, D>(output, direction, ref);
}

#endif
