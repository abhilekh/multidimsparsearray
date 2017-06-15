/**
 * This file is part of the MultiDimSparseArray library
 *
 * @license  BSD-3
 * @author   Abhilekh Agarwal
 */


#pragma once

#include <array>
#include <iostream>
#include "exception.h"

template<unsigned N>
class mdrange {
private:
    std::array<int, N> last;
    std::array<int, N> first, iter;
    bool overflow;
    int step;

    void init() {
        for (auto i = 0U; i < N; i++) {
            iter[i] = first[i];
            last[i] -= 1;
            if (last[i] < first[i]) {
                throw InvalidDimensionsException("Wrong Dimension");
            }
        }
        overflow = false;
    }

public:
    mdrange(std::array<int, N> end) :
            last(end), step(1) {
        first.fill(0);
        init();
    }

    mdrange(std::array<int, N> start, std::array<int, N> end) :
            last(end), first(start), step(1) {
        init();
    }

    const mdrange& begin() const {
        return *this;
    }
    const mdrange& end() const {
        return *this;
    }

    // Iterator functions
    bool operator!=(const mdrange&) const {
        return !overflow;
    }

    void operator++() {
        for (int i = N - 1; i > -1; i--) {
            if (last[i] == iter[i]) {
                continue;
            }
            iter[i] += step;
            if (i == N - 1)
                return;
            for (unsigned j = i + 1; j < N; j++) {
                iter[j] = first[j];
            }
            return;
        }
        overflow = true;
    }
    const std::array<int, N>& operator*() const {
        return iter;
    }
};
