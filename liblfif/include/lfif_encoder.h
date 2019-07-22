/**
* @file lfif_encoder.h
* @author Drahomír Dlabaja (xdlaba02)
* @date 12. 5. 2019
* @copyright 2019 Drahomír Dlabaja
* @brief Functions for encoding an image.
*/

#ifndef LFIF_ENCODER_H
#define LFIF_ENCODER_H

#include "block_compress_chain.h"
#include "cabac_contexts.h"
#include "bitstream.h"

#include <cstdint>
#include <ostream>
#include <sstream>

template<size_t D>
struct LfifEncoderParams {
  enum QuantTableType {
    DEFAULT = UNIFORM,
    DCTDIAG,
    DCTCOPY,
    FILLDIAG,
    FILLCOPY,
    UNIFORM,
    CUSTOM
  };

  enum TraversalTableType {
    DEFAULT = ZIGZAG,
    REFERENCE,
    ZIGZAG,
    QTABLE,
    RADIUS,
    DIAGONALS,
    HYPERBOLOID,
    CUSTOM
  };

  uint64_t img_dims[D+1] {}; /**< @brief Dimensions of a decoded image + image count. The multiple of all values should be equal to number of pixels in encoded image.*/
  uint8_t  color_depth   {}; /**< @brief Number of bits per sample used by each decoded channel.*/
  bool     use_huffman   {}; /**< @brief Huffman encoding will be used when this is true.*/
  uint8_t  block_size    {};

  std::vector<QuantTableType>     quant_table_types     {};
  std::vector<double>             qualities             {};
  std::vector<TraversalTableType> traversal_table_types {};

  uint8_t num_entropy_encoder_ctxs {};

  uint8_t num_channels {};

  std::vector<uint8_t> quant_table_idxs     {};
  std::vector<uint8_t> traversal_table_idxs {};
  std::vector<uint8_t> entropy_encoder_idxs {};
}

/**
* @brief Base structure for encoding an image.
*/
template<size_t D>
class LfifEncoder {
public:
  inline LfifEncoder(const LfifEncoderParams &params);

  const LfifEncoderParams &params;

  std::vector<QuantTable<D>>     quant_table;      /**< @brief Quantization matrices for luma and chroma.*/
  std::vector<TraversalTable<D>> traversal_table;  /**< @brief Traversal matrices for luma and chroma.*/

  //TODO ty dve hod do tridy kontext, encoder udelej univerzalni
  std::vector<ReferenceBlock<D>>             reference_block; /**< @brief Reference blocks for luma and chroma to be used for traversal optimization.*/
  std::vector<std::array<HuffmanWeights, 2>> huffman_weight; /**< @brief Weigth maps for luma and chroma and also for the DC and AC coefficients to be used for optimal Huffman encoding.*/
  std::vector<std::array<HuffmanEncoder, 2>> huffman_encoder; /**< @brief Huffman encoders for luma and chroma and also for the DC and AC coefficients.*/

  size_t blocks_cnt; /**< @brief Number of blocks in the encoded image.*/
  size_t pixels_cnt; /**< @brief Number of pixels in the encoded image.*/

  size_t amp_bits;    /**< @brief Number of bits sufficient to contain maximum DCT coefficient.*/
  size_t zeroes_bits; /**< @brief Number of bits sufficient to contain run-length of zero coefficients.*/
  size_t class_bits;  /**< @brief Number of bits sufficient to contain number of bits of maximum DCT coefficient.*/
  size_t max_zeroes;  /**< @brief Maximum length of zero coefficient run-length.*/

  /**
  * @brief Function which initializes quantization matrices.
  * @param table_type Type of the tables to be initialized.
  * @param quality Encoding quality from 1 to 100.
  * @return Nonzero if requested unknown table type.
  */
  inline int constructQuantizationTables(const std::string &table_type, float quality)

  /**
  * @brief Function which performs arbitrary scan of an image. This function prepares a block into the encoding structure buffers.
  * @param input Callback function specified by client with signature T input(size_t index), where T is pixel/sample type.
  * @param func Callback function which performs action on every block with signature void func(size_t channel), where channel is channel from which the extracted block is.
  */
  template<typename IF, typename F>
  inline void performScan(IF &&input, F &&func);

  /**
  * @brief Function which performs scan for linearization optimization.
  * @param input Callback function specified by client with signature T input(size_t index), where T is pixel/sample type.
  */
  template<typename F>
  void referenceScan(F &&input);

  /**
  * @brief Function which initializes linearization matrices.
  * @param table_type Type of the tables to be initialized.
  * @return Nonzero if requested unknown table type.
  */
  int constructTraversalTables(const std::string &table_type);

  /**
  * @brief Function which performs scan for huffman weights.
  * @param input Callback function specified by client with signature T input(size_t index), where T is pixel/sample type.
  */
  template<typename F>
  void huffmanScan(F &&input);

  /**
  * @brief Function which initializes Huffman encoder from weights.
  */
  void constructHuffmanTables();

  /**
  * @brief Function which writes tables to stream.
  * @param output Stream to which data will be written.
  */
  void writeHeader(std::ostream &output);

  /**
  * @brief Function which encodes the image to the stream.
  * @param input Callback function specified by client with signature T input(size_t index), where T is pixel/sample type.
  * @param output Output stream to which the image shall be encoded.
  */
  template<typename F>
  void outputScanHuffman(F &&input, std::ostream &output);

  /**
  * @brief Function which encodes the image to the stream with CABAC.
  * @param input Callback function specified by client with signature T input(size_t index), where T is pixel/sample type.
  * @param output Output stream to which the image shall be encoded.
  */
  template<typename F>
  void outputScanCABAC(F &&input, std::ostream &output)
};

template<size_t D>
LfifEncoder::LfifEncoder(const LfifEncoderParams &params): params { params }, dct(params.block_size) {
  quant_table = new QuantTable[params.num_quant_tables];
  traversal_table = new TraversalTable[params.num_traversal_tables];
}

template<size_t D>
void LfifEncoder<D>::initEncoder() {
  quant_tables[0]     = &quant_table[0];
  quant_tables[1]     = &quant_table[1];
  quant_tables[2]     = &quant_table[1];

  reference_blocks[0] = &reference_block[0];
  reference_blocks[1] = &reference_block[1];
  reference_blocks[2] = &reference_block[1];

  traversal_tables[0] = &traversal_table[0];
  traversal_tables[1] = &traversal_table[1];
  traversal_tables[2] = &traversal_table[1];

  huffman_weights[0]  =  huffman_weight[0];
  huffman_weights[1]  =  huffman_weight[1];
  huffman_weights[2]  =  huffman_weight[1];

  huffman_encoders[0] =  huffman_encoder[0];
  huffman_encoders[1] =  huffman_encoder[1];
  huffman_encoders[2] =  huffman_encoder[1];

  blocks_cnt = 1;
  pixels_cnt = 1;

  for (size_t i = 0; i < D; i++) {
    blocks_cnt *= ceil(img_dims[i] / static_cast<double>(block_size));
    pixels_cnt *= img_dims[i];
  }

  amp_bits = ceil(log2(pow(block_size, D))) + color_depth - D - (D/2);
  class_bits = RunLengthPair::classBits(amp_bits);
  zeroes_bits = RunLengthPair::zeroesBits(class_bits);
  max_zeroes = pow(2, zeroes_bits) - 1;

  for (size_t y = 0; y < 2; y++) {
    reference_block[y].fill(0);
    for (size_t x = 0; x < 2; x++) {
      huffman_weight[y][x].clear();
    }
  }
}

template<size_t D>
int LfifEncoder<D>::constructQuantizationTables(QuantTableType table_type, float quality) {
  auto scale_table = [&](const auto &table) {
    if (table_type == DCTDIAG)  return    scaleByDCT<2>(table, block_size);
    if (table_type == DCTCOPY)  return    scaleByDCT<2>(table, block_size);
    if (table_type == FILLDIAG) return scaleFillNear<2>(table, block_size);
    if (table_type == FILLCOPY) return scaleFillNear<2>(table, block_size);
    return QuantTable<2>();
  };

  auto extend_table = [&](const auto &table) {
    if (table_type == DCTDIAG)  return averageDiagonalTable<2, D>(table, block_size);
    if (table_type == DCTCOPY)  return            copyTable<2, D>(table, block_size);
    if (table_type == FILLDIAG) return averageDiagonalTable<2, D>(table, block_size);
    if (table_type == FILLCOPY) return            copyTable<2, D>(table, block_size);
    return QuantTable<D>();
  };

  if (table_type == DCTDIAG || table_type == FILLDIAG || table_type == DCTCOPY || table_type == FILLCOPY) {
    quant_table[0] = extend_table(scale_table(base_luma));
    quant_table[1] = extend_table(scale_table(base_chroma));
  }
  else if (table_type == DEFAULT || table_type == UNIFORM) {
    for (size_t i = 0; i < 2; i++) {
      quant_table[i] = uniformTable<D>(50);
    }
  }
  else {
    return -1;
  }

  for (size_t i = 0; i < 2; i++) {
    quant_table[i] = clampTable<D>(applyQuality<D>(quant_table[i], quality), 1, 255);
  }

  return 0;
}

template<size_t D>
template<typename IF, typename F>
void LfifEncoder<D>::performScan(IF &&input, F &&func) {
  Block<INPUTTRIPLET, D> current_block(m_block_size);
  Block<INPUTUNIT,    D>   input_block(m_block_size);

  auto outputF = [&](size_t index, const auto &value) {
    current_block[index] = value;
  };

  for (size_t img = 0; img < img_dims[D]; img++) {
    auto inputF = [&](size_t index) {
      return input(img * pixels_cnt + index);
    };

    for (size_t block = 0; block < blocks_cnt; block++) {
      getBlock<D>(inputF, block, img_dims, outputF);

      for (size_t channel = 0; channel < params.num_channels; channel++) {

        for (size_t i = 0; i < pow(m_block_size, D); i++) {
          input_block[i] = current_block[i][channel];
        }

        func(input_block, channel);
      }
    }
  }
}

template <size_t D>
template<typename F>
void LfifEncoder<D>::referenceScan(F &&input) {
  DCT                               dct(m_block_size);
  Block<DCTDATAUNIT, D>       dct_block(m_block_size);
  Block<QDATAUNIT,   D> quantized_block(m_block_size);

  auto compress_chain = [&](Block<INPUTUNIT, D> &input_block, size_t channel) {
    forwardDiscreteCosineTransform<D>(input_block,      dct_block,        dct);
                          quantize<D>(dct_block,        quantized_block,  quant_table[params.quant_table_idxs[channel]]);
               addToReferenceBlock<D>(quantized_block, reference_block[params.traversal_table_idxs[channel]]);
  };

  performScan(input, compress_chain);
}

template<size_t D>
int LfifEncoder<D>::constructTraversalTables(TraversalTableType table_type) {
  if (table_type == REFERENCE) {
    for (size_t i = 0; i < 2; i++) {
      traversal_table[i]
      .constructByReference(reference_block[i]);
    }
  }
  else if (table_type == DEFAULT || table_type == ZIGZAG) {
    for (size_t i = 0; i < 2; i++) {
      traversal_table[i]
      .constructZigzag();
    }
  }
  else if (table_type == QTABLE) {
    for (size_t i = 0; i < 2; i++) {
      traversal_table[i]
      .constructByQuantTable(quant_table[i]);
    }
  }
  else if (table_type == RADIUS) {
    for (size_t i = 0; i < 2; i++) {
      traversal_table[i]
      .constructByRadius();
    }
  }
  else if (table_type == DIAGONALS) {
    for (size_t i = 0; i < 2; i++) {
      traversal_table[i]
      .constructByDiagonals();
    }
  }
  else if (table_type == HYPERBOLOID) {
    for (size_t i = 0; i < 2; i++) {
      traversal_table[i]
      .constructByHyperboloid();
    }
  }
  else {
    return -1;
  }

  return 0;
}

template<size_t D>
template<typename F>
void LfifEncoder<D>::huffmanScan(F &&input) {
  DCT                                 dct(m_block_size);
  Block<DCTDATAUNIT,   D>       dct_block(m_block_size);
  Block<QDATAUNIT,     D> quantized_block(m_block_size);
  Block<RunLengthPair, D>       runlength(m_block_size);

  QDATAUNIT previous_DC [3] {};

  auto compress_chain = [&](Block<INPUTUNIT, D> &input_block, size_t channel) {
    forwardDiscreteCosineTransform<D>(input_block,     dct_block);
                          quantize<D>(dct_block,       quantized_block,         quant_table[params.quant_table_idxs[channel]]);
                      diffEncodeDC<D>(quantized_block, previous_DC[channel]);
                          traverse<D>(quantized_block, traversal_tables[params.traversal_table_idxs[channel]]);
                   runLengthEncode<D>(quantized_block, runlength,                max_zeroes);
                 huffmanAddWeights<D>(runlength,       huffman_weights[params.entropy_encoder_idxs[channel]], class_bits);
  };

  performScan(input, compress_chain);
}

template<size_t D>
void LfifEncoder<D>::constructHuffmanTables() {
  for (size_t y = 0; y < 2; y++) {
    for (size_t x = 0; x < 2; x++) {
      huffman_encoder[y][x]
      . generateFromWeights(huffman_weight[y][x]);
    }
  }
}

template<size_t D>
void LfifEncoder<D>::writeHeader(std::ostream &output) {
  output << "LFIF-" << D << "D\n";

  writeValueToStream<uint8_t>(block_size, output);
  writeValueToStream<uint8_t>(color_depth, output);

  for (size_t i = 0; i < D + 1; i++) {
    writeValueToStream<uint64_t>(img_dims[i], output);
  }

  for (size_t i = 0; i < params.num_quant_tables; i++) {
    writeToStream<D>(quant_table[i], output);
  }

  for (size_t i = 0; i < 2; i++) {
    traversal_table[i]
    . writeToStream(output);
  }

  writeValueToStream<uint8_t>(use_huffman, output);

  if (use_huffman) {
    for (size_t y = 0; y < 2; y++) {
      for (size_t x = 0; x < 2; x++) {
        huffman_encoder[y][x]
        . writeToStream(output);
      }
    }
  }
}

template<size_t D>
template<typename F>
void LfifEncoder<D>::outputScanHuffman(F &&input, std::ostream &output) {
  DCT                                 dct(m_block_size);
  Block<DCTDATAUNIT,   D>       dct_block(m_block_size);
  Block<QDATAUNIT,     D> quantized_block(m_block_size);
  Block<RunLengthPair, D>       runlength(m_block_size);

  QDATAUNIT previous_DC [3] {};
  OBitstream bitstream      {};

  bitstream.open(&output);

  auto compress_chain = [&](Block<INPUTUNIT, D> &input_block, size_t channel) {
    forwardDiscreteCosineTransform<D>(input_block,      dct_block);
                          quantize<D>(dct_block,        quantized_block,          *quant_tables[channel]);
                      diffEncodeDC<D>(quantized_block,  previous_DC[channel]);
                          traverse<D>(quantized_block, *traversal_tables[channel]);
                   runLengthEncode<D>(quantized_block,  runlength,                 max_zeroes);
             encodeToStreamHuffman<D>(runlength,        huffman_encoders[channel], bitstream, class_bits);
  };

  performScan(input, compress_chain);

  bitstream.flush();
}

template<size_t D>
template<typename F>
void LfifEncoder<D>::outputScanCABAC(F &&input, std::ostream &output) {
  Block<DCTDATAUNIT, D>       dct_block(m_block_size);
  Block<QDATAUNIT,   D> quantized_block(m_block_size);
  Block<RunLengthPair, D>     runlength(m_block_size);

  QDATAUNIT            previous_DC[3] {};
  CABACContexts<D>     contexts   [3] {};
  OBitstream           bitstream      {};
  CABACEncoder         cabac          {};

  bitstream.open(&output);
  cabac.init(bitstream);

  auto perform = [&](size_t channel) {
    forwardDiscreteCosineTransform<D>(input_block,      dct_block);
                          quantize<D>(dct_block,        quantized_block,          *quant_tables[channel]);
                      diffEncodeDC<D>(quantized_block,  previous_DC[channel]);
                          traverse<D>(quantized_block, *traversal_tables[channel]);
              encodeTraversedCABAC<D>(quantized_block,  cabac,                      contexts[channel]);
  };

  performScan(input, perform);

  cabac.terminate();
}

#endif
