/**
* @file cabac_contexts.h
* @author Drahomír Dlabaja (xdlaba02)
* @date 16. 7. 2019
* @copyright 2019 Drahomír Dlabaja
* @brief Context-Adaptive Binary Arithmetic Coding.
*/

#ifndef CABAC_CONTEXTS_H
#define CABAC_CONTEXTS_H

template <size_t D>
class CABACContexts {
public:
  CABACContexts(size_t block_size):
  significant_coef_flag_ctx(block_size),
  last_significant_coef_flag_ctx(block_size),
  m_block_size { block_size } {}

  static const size_t NUM_GREATER_ONE_CTXS { 5 };
  static const size_t NUM_ABS_LEVEL_CTXS   { 5 };

  CABAC::ContextModel           coded_block_flag_ctx;
  Block<CABAC::ContextModel, D> significant_coef_flag_ctx;
  Block<CABAC::ContextModel, D> last_significant_coef_flag_ctx;
  CABAC::ContextModel           coef_greater_one_ctx[NUM_GREATER_ONE_CTXS];
  CABAC::ContextModel           coef_abs_level_ctx[NUM_ABS_LEVEL_CTXS];

  size_t size() { return m_block_size; }

private:
  size_t m_block_size;
};

#endif
