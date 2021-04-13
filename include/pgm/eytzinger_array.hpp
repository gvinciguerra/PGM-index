//
// Created by Роман Агеев on 4/12/21.
//

#ifndef PIECEWISEGEOMETRICMODELINDEX_EYTZINGER_ARRAY_HPP
#define PIECEWISEGEOMETRICMODELINDEX_EYTZINGER_ARRAY_HPP
#include <iostream>
#include <algorithm>

template<typename T>
class EytzingerArray {
  static const size_t multiplier = 64 / sizeof(T);
  static const size_t offset = multiplier + (multiplier >> 1);

  T* array = nullptr;
  size_t _size = 0;
  template<typename ForwardIterator>
  ForwardIterator reorder_data(ForwardIterator iter, size_t pos);

 public:
  EytzingerArray() = default;

  EytzingerArray(EytzingerArray const &other) : EytzingerArray() {
    auto copy = EytzingerArray(other.array, _size);
    swap(copy);
  }

  EytzingerArray &operator=(EytzingerArray const &other) {
    auto copy = EytzingerArray(other.array, _size);
    swap(copy);
    return *this;
  }

  EytzingerArray(EytzingerArray &&other) noexcept : EytzingerArray()  {
    swap(other);
  }

  EytzingerArray &operator=(EytzingerArray &&other) noexcept {
    swap(other);
    return *this;
  }
  ~EytzingerArray();

  size_t size() const {
    return _size;
  }

  template<typename ForwardIterator>
  EytzingerArray(ForwardIterator iter, size_t size);

  T & operator[](size_t idx) {
    return array[idx];
  }

  T const & operator[](size_t idx) const {
    return array[idx];
  }

  template<bool prefetch = true>
  size_t search(T x) const;

  void swap(EytzingerArray & other) noexcept {
    std::swap(array, other.array);
    std::swap(_size, other._size);
  }
};

template<typename T>
template<typename ForwardIterator>
ForwardIterator EytzingerArray<T>::reorder_data(ForwardIterator iter, size_t pos) {
  if (pos < _size) {
    iter = reorder_data(iter, (pos << 1u) + 1u);
    array[pos] = (*iter);
    ++iter;
    iter = reorder_data(iter, (pos << 1u) + 2u);
  }
  return iter;
}

template<typename T>
template<typename ForwardIterator>
EytzingerArray<T>::EytzingerArray(ForwardIterator iter, size_t size) : _size(size) {
  array = new T[size];
  reorder_data(iter, 0);
}

template<typename T>
template<bool builtin>
size_t EytzingerArray<T>::search(T x) const {
  size_t i = 0;
  while (i < _size) {
    if constexpr (builtin) __builtin_prefetch(array + (multiplier * i + offset));
    i = (x < array[i]) ? ((i << 1u) + 1u) : ((i << 1u) + 2u);
  }
  size_t j = (i + 1u) >> __builtin_ffs(~(i + 1u));
  return (j == 0) ? _size : j - 1;
}

template<typename T>
EytzingerArray<T>::~EytzingerArray() {
  delete[] array;
}

#endif //PIECEWISEGEOMETRICMODELINDEX_EYTZINGER_ARRAY_HPP
