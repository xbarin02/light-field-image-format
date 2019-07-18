/**
* @file stack_allocator.h
* @author Drahomír Dlabaja (xdlaba02)
* @date 18. 5. 2019
* @copyright 2019 Drahomír Dlabaja
* @brief Module for allocating memory on stack dynamically.
*/

#ifndef STACK_ALLOCATOR_H
#define STACK_ALLOCATOR_H

class StackAllocator {
public:
  static void resize(size_t size) {
    m_memory.resize(size);
  }

  char *malloc(size_t size) {
    char *mem_ptr = &m_memory[m_stack_ptr];
    m_stack_ptr += size;
    return mem_ptr;
  }

  void free(size_t size) {
    m_stack_ptr -= size;
  }

private:
  static size_t m_stack_ptr;
  static std::vector<char> m_memory;
}

#endif
