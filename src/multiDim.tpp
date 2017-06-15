/**
 * This file is part of the MultiDimSparseArray library
 *
 * @license  BSD-3
 * @author   Abhilekh Agarwal
 */


#pragma once

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>
#include "patch/patch.hpp"
#include "SparseDim.tpp"
#include "exception.h"
#include "MultiDimIter.tpp"

template<typename T, unsigned N = 2>
class MultiDim: public SparseDim<T, N> {

public:
    MultiDim(const std::array<int, N> &dims_sz, T defval = T()) :
            SparseDim<T, N>(dims_sz, defval) {  // Constructor 1
    }

    MultiDim(const SparseDim<T, N> & matrix) :
            SparseDim<T, N>(matrix) {  // Constructor 2
    }

    MultiDim<T, N> & operator =(const MultiDim<T, N> & matrix) {
        this->printallowed = 0;
        if (this->printallowed)
            std::cout << __func__ << __LINE__ << std::endl;
        if (&matrix != this) {
            this->deepCopy(matrix);
        }

        return *this;
    }

    MultiDim<T, N> & operator =(const SparseDim<T, N> & matrix) {
        this->printallowed = 0;
        if (this->printallowed)
            std::cout << __func__ << __LINE__ << std::endl;
        if (&matrix != this) {
            this->deepCopy(matrix);
        }

        return *this;
    }

// === OPERATORS ==============================================

    MultiDim<T, N> operator +(const T &t) {
        return _oper([&t](T& i) {i+=t;}, false);
    }

    MultiDim<T, N> operator +(const MultiDim<T, N> &sd) {
        if (this->printallowed)
            std::cout << __func__ << __LINE__ << std::endl;
        if (this->_dims_sz != sd._dims_sz) {
            throw InvalidDimensionsException(
                    "Cannot add: matrices dimensions don't match.");
        }
        auto defval = this->_defval + sd._defval;
        MultiDim<T, N> result(*this);
        result._defval = defval;

        if (sd._defval == T()) {
            // Now result equlas part 1, and can only change when there is
            // a val in part b. Another assumption is
            // b's default value is 0 which makes add/subs as 1st val
            for (int _v : patch::xrange(this->_virdimsz)) {
                auto& rowlst = sd.rowlst_lst[_v];
                auto& collst = sd.collst_lst[_v];
                auto& vallst = sd.vallst_lst[_v];
                for (int row : patch::xrange(this->_dims_sz[N - 2])) {
                    for (int pos = rowlst[row]; pos < rowlst[row + 1]; pos++) {
                        auto chval = vallst.at(pos);
                        result._setop(_v, row, collst[pos],
                                [&chval](T i) {return (i + chval);});
                    }
                }
            }
        } else {
            // General case
            int vgrp = -1;
            for (auto arr : mdrange<N>(this->_dims_sz)) {
                if (!arr[N - 2] && !arr[N - 1]) {
                    vgrp++;
                }
                auto oval = this->_get(arr, vgrp);
                auto nval = oval + sd._get(arr, vgrp);
                if (nval != oval)
                    result._set(nval, arr, vgrp);
            }
        }
        return result;
    }

    MultiDim<T, N> operator -(const T &t) const {
        return _oper([&t](T& i) {i-=t;}, false);
    }

    MultiDim<T, N> operator -(const MultiDim<T, N> &sd) {
        if (this->printallowed)
            std::cout << __func__ << __LINE__ << std::endl;
        if (this->_dims_sz != sd._dims_sz) {
            throw InvalidDimensionsException(
                    "Cannot sub: matrices dimensions don't match.");
        }
        auto defval = this->_defval - sd._defval;
        MultiDim<T, N> result(*this);
        result._defval = defval;
        if (sd._defval == T()) {
            // Now result equals part 1, and can only change when there is
            // a val in part b. Another assumption is
            // b's default value is 0 which makes add/subs as 1st val
            for (int _v : patch::xrange(this->_virdimsz)) {
                auto& rowlst = sd.rowlst_lst[_v];
                auto& collst = sd.collst_lst[_v];
                auto& vallst = sd.vallst_lst[_v];
                for (int row : patch::xrange(this->_dims_sz[N - 2])) {
                    for (int pos = rowlst[row]; pos < rowlst[row + 1]; pos++) {
                        auto chval = vallst.at(pos);
                        result._setop(_v, row, collst[pos],
                                [&chval](T i) {return (i - chval);});
                    }
                }
            }
        } else {
            // General case
            int vgrp = -1;
            for (auto arr : mdrange<N>(this->_dims_sz)) {
                if (!arr[N - 2] && !arr[N - 1]) {
                    vgrp++;
                }
                auto oval = this->_get(arr, vgrp);
                auto nval = oval - sd._get(arr, vgrp);
                if (nval != oval)
                    result._set(nval, arr, vgrp);
            }
        }
        return result;
    }

    MultiDim<T, N> operator *(const T &t) const {
        if (t == T() && (this->vallst_lst[0].size() > 0)
                && (t * this->vallst_lst[0][0] == t)) {
            // This means either t is some thing like zero or
            // 1st item is identity element in my domain.
            // Making idempotent does not make sense.
            MultiDim<T, N> result(*this);
            result._defval = this->_defval * t;

            for (auto& vallst : result.vallst_lst)
                vallst.clear();
            for (auto& collst : result.collst_lst)
                collst.clear();
            for (auto& rowlst : result.rowlst_lst)
                std::fill(rowlst.begin(), rowlst.end(), 0);
            return result;
        } else {
            // Simple case
            return this->_oper([&t](T& i) {i*=t;}, false);
        }
    }

    MultiDim<T, N> operator /(const T &t) const {
        return this->_oper([&t](T& i) {i/=t;}, false);
    }

    MultiDim<T, N> operator %(const T &t) const {
        return this->_oper([&t](T& i) {i%=t;}, false);
    }

    template<typename UniOpType>
    MultiDim<T, N> oper(UniOpType op) const {
        return this->_oper(op, true);
    }

    // === FRIEND FUNCTIONS =========================================

    friend std::ostream & operator <<(std::ostream & os,
            const MultiDim<T, N> & sd) {
        int vgrp = -1;
        for (auto arr : mdrange<N>(sd._dims_sz)) {
            if (!arr[N - 1]) {
                os << "\n";
            }
            if (!arr[N - 2] && !arr[N - 1]) {
                os << "\nCase:";
                for (int a : patch::xrange(N - 1)) {
                    os << arr[a] << ",";
                }
                os << arr[N - 1] << "\n";
                vgrp++;
            }
            if (arr[N - 1]) {
                os << " ";
            }
            os << sd._get(arr, vgrp);
        }
        return os;
    }

protected:
    void deepCopy(const MultiDim<T, N> & md) {
        if (this->printallowed)
            std::cout << __func__ << __LINE__ << std::endl;
        this->_virdimsz = md._virdimsz;
        this->_defval = md._defval;
        std::copy(std::begin(md._dims_sz), std::end(md._dims_sz),
                std::begin(this->_dims_sz));
        this->vallst_lst = std::vector<std::vector<T> >(this->_virdimsz);
        this->collst_lst = std::vector<std::vector<int> >(this->_virdimsz);
        this->rowlst_lst = std::vector<std::vector<int> >(this->_virdimsz);
        for (int i : patch::xrange(this->_virdimsz)) {
            this->vallst_lst[i] = std::vector<T>(md.vallst_lst[i]);
            this->vallst_lst[i].reserve(10);
            this->collst_lst[i] = std::vector<int>(md.collst_lst[i]);
            this->collst_lst[i].reserve(10);
            // We make row one size larger in CSR
            this->rowlst_lst[i] = std::vector<int>(md.rowlst_lst[i]);
        }
    }

private:
    template<typename OpType>
    inline void _setop(int _virtualidx, int row, int col, OpType op) {
        // No bound check. It will make it slower
        if (this->printallowed > 1)
            std::cout << __func__ << __LINE__ << std::endl;
        int currCol = -1;
        auto& rowlst = this->rowlst_lst[_virtualidx];
        auto& collst = this->collst_lst[_virtualidx];
        int pos;
        for (pos = rowlst[row]; pos < rowlst[row + 1]; pos++) {
            currCol = collst[pos];

            if (currCol >= col) {
                break;
            }
        }
        if (currCol != col) {
            T val = op(T());
            if (val != this->_defval) {
                this->insert(_virtualidx, pos, row, col, val);
            }

        } else {
            T val = op(this->vallst_lst[_virtualidx][pos]);
            if (val == this->_defval) {
                this->remove(_virtualidx, pos, row);

            } else {
                this->vallst_lst[_virtualidx][pos] = val;
            }
        }
    }
};
