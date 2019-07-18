/**
* @file dct.h
* @author Drahomír Dlabaja (xdlaba02)
* @date 12. 5. 2019
* @copyright 2019 Drahomír Dlabaja
* @brief Functions which performs FDCT and IDCT.
*/

#ifndef DCT_H
#define DCT_H

#include "pow.h"

#include <cmath>
#include <array>

/**
* @brief Type used for the DCT computations.
*/
using DCTDATAUNIT = float;

class DCT {
public:
  DCT(size_t block_size): m_coefs(block_size) {
    for (size_t u = 0; u < block_size; u++) {
      for (size_t x = 0; x < block_size; x++) {
        m_coefs[u * block_size + x] = cos(((2 * x + 1) * u * M_PI ) / (2 * block_size)) * sqrt(2.0 / block_size) * (u == 0 ? (1 / sqrt(2)) : 1);
      }
    }
  }

  size_t blockSize() const {
    return m_coefs.size();
  }

  /**
  * @brief Struct for the FDCT which wraps static parameters for partial specialization.
  */
  template <size_t D>
  struct fdct {

    /**
    * @brief Function which performs the FDCT.
    * @param input  callback function returning samples from a block. Signature is DCTDATAUNIT input(size_t index).
    * @param output callback function for writing DCT coefficients. Signature is DCTDATAUNIT &output(size_t index).
    */
    template <typename IF, typename OF>
    fdct(IF &&input, OF &&output) {
      Block<DCTDATAUNIT, D> tmp(m_coefs.size());

      for (size_t slice = 0; slice < m_coefs.size(); slice++) {
        auto inputF = [&](size_t index) {
          return input(slice * pow(m_coefs.size(), D - 1) + index);
        };

        auto outputF = [&](size_t index) -> auto & {
          return tmp[slice * pow(m_coefs.size(), D - 1) + index];
        };

        fdct<D - 1>(inputF, outputF);
      }

      for (size_t noodle = 0; noodle < pow(m_coefs.size(), D - 1); noodle++) {
        auto inputF = [&](size_t index) {
          return tmp[index * pow(m_coefs.size(), D - 1) + noodle];
        };

        auto outputF = [&](size_t index) -> auto & {
          return output(index * pow(m_coefs.size(), D - 1) + noodle);
        };

        fdct<1>(inputF, outputF);
      }
    }
  };

  /**
   * @brief The parital specialization for performing one dimensional fdct.
   * @see fdct<D>
   */
  template<>
  struct fdct<1> {

    /**
     * @brief The parital specialization for performing one dimensional fdct.
     * @see fdct<D>::fdct
     */
    template <typename IF, typename OF>
    fdct(IF &&input, OF &&output) {
      for (size_t u = 0; u < m_coefs.size(); u++) {
        for (size_t x = 0; x < m_coefs.size(); x++) {
          output(u) += input(x) * m_coefs[u * m_coefs.size() + x];
        }
      }
    }
  };

  /**
  * @brief Struct for the IDCT which wraps static parameters for partial specialization.
  */
  template <size_t D>
  struct idct {

    /**
    * @brief Function which performs the IDCT.
    * @param input  callback function returning coefficients from decoded block. Signature is DCTDATAUNIT input(size_t index).
    * @param output callback function for writing output samples. Signature is DCTDATAUNIT &output(size_t index).
    */
    template <typename IF, typename OF>
    idct(IF &&input, OF &&output) {
      Block<DCTDATAUNIT, D> tmp(m_coefs.size());

      for (size_t slice = 0; slice < m_coefs.size(); slice++) {
        auto inputF = [&](size_t index) {
          return input(slice * pow(m_coefs.size(), D - 1) + index);
        };

        auto outputF = [&](size_t index) -> auto & {
          return tmp[slice * pow(m_coefs.size(), D - 1) + index];
        };

        idct<D - 1>(inputF, outputF);
      }

      for (size_t noodle = 0; noodle < pow(m_coefs.size(), D - 1); noodle++) {
        auto inputF = [&](size_t index) {
          return tmp[index * pow(m_coefs.size(), D - 1) + noodle];
        };

        auto outputF = [&](size_t index) -> auto & {
          return output(index * pow(m_coefs.size(), D - 1) + noodle);
        };

        idct<1>(inputF, outputF);
      }
    }
  };

  /**
   * @brief The parital specialization for putting one dimensional idct.
   * @see idct<BS, D>
   */
  template <>
  struct idct<1> {

    /**
     * @brief The parital specialization for putting one dimensional idct.
     * @see idct<BS, D>::idct
     */
    template <typename IF, typename OF>
    idct(IF &&input, OF &&output) {
      for (size_t x = 0; x < m_coefs.size(); x++) {
        for (size_t u = 0; u < m_coefs.size(); u++) {
          output(x) += input(u) * m_coefs[u * m_coefs.size() + x];
        }
      }
    }
  };

private:
  Block<DCTDATAUNIT, 2> m_coefs;
};

#endif
