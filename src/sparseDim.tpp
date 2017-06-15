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
#include "exception.h"
#include "multiDimIter.tpp"

#define UNUSED(expr) do { (void)(expr); } while (0)

/*
 This is implementation of N dimensional array(N >=2). I could have done some
 thing fancy but CSR is alreay a good format. So I took last two indexes as
 row and col and prev to get a virdim. So if I have a array of dimension
 a[7][3][5][4]. So virdimsz = 7*3 =21 ( 0 to N-1-2 ) and I will have multiple row
 and col vector this will also keep insertion and deletion from vector fast.
 */

template<typename T, unsigned N = 2>
class SparseDim {
public:

    // === CREATION ==============================================
    SparseDim(const std::array<int, N> &dims_sz, T defval = T()) :
            _dims_sz(dims_sz) {
        static_assert(N>1, "Size cannot be less than 2");
        printallowed = 0;
        if (printallowed)
            std::cout << __func__ << __LINE__ << std::endl;
        _virdimsz = 1;
        for (int i : patch::xrange(N - 2)) {
            _virdimsz *= dims_sz[i];
        }
        for (int i : patch::xrange(N)) {
            if (dims_sz[i] < 1) {
                throw InvalidCoordinatesException(
                        "Coordinates cannot be 0 or -ve");
            }
        }
        _defval = defval;
        vallst_lst = std::vector<std::vector<T> >(_virdimsz);
        collst_lst = std::vector<std::vector<int> >(_virdimsz);
        rowlst_lst = std::vector<std::vector<int> >(_virdimsz);
        for (int i : patch::xrange(_virdimsz)) {
            vallst_lst[i].reserve(10);
            collst_lst[i].reserve(10);
            // We make row one size larger in CSR
            rowlst_lst[i] = std::vector<int>(dims_sz[N - 2] + 1, 0);
        }
    }

    SparseDim(const SparseDim<T, N> & matrix) { // copy constructor
        printallowed = 0;
        if (printallowed)
            std::cout << __func__ << __LINE__ << std::endl;
        deepCopy(matrix);
    }

    SparseDim<T, N> & operator =(const SparseDim<T, N> & matrix) {
        printallowed = 0;
        if (printallowed)
            std::cout << __func__ << __LINE__ << std::endl;
        if (&matrix != this) {
            this->deepCopy(matrix);
        }

        return *this;
    }

    void loglvl(int lvl) {
        printallowed = lvl;
    }

    // === VALUES ==============================================

    /*
     It does not get reference, as any change can make value become def value
     and vice versa. So I am only returning const value. This is obviously at
     performance cost.
     */
    T get(const std::array<int, N> &atloc) const {
        if (printallowed)
            std::cout << __func__ << __LINE__ << std::endl;
        return _get(atloc, -1);
    }

    void set(const T &val, const std::array<int, N> &atloc) {
        if (printallowed)
            std::cout << __func__ << __LINE__ << "val" << val << std::endl;
        _set(val, atloc, -1);
    }

    inline int getdimensions() const {
        return N;
    }

    inline unsigned getdimensionsSize(unsigned d) const {
        return (d < N) ? _dims_sz[d] : 0;
    }

    // === Disk save/retrieve ======================================
    int dump(const std::string &path) const {
        std::ofstream fout(path.c_str(), std::ios::out | std::ios::binary);
        if (!fout) {
            std::cout << "error dump" << 1 <<std::endl;;
            return 1;
        }

        if (!fout.good()) {
            std::cout << "error dump" << 2 <<std::endl;;
            return 2;
        }
        fout.write((char*) &_defval, sizeof(T));

        int sz = _dims_sz.size();
        fout.write((char*) &sz, sizeof(int));
        fout.write((char*) &_dims_sz[0], sz * sizeof(_dims_sz[0]));
        dumpVOfV(fout, rowlst_lst);
        dumpVOfV(fout, collst_lst);
        dumpVOfV(fout, vallst_lst);

        fout.close();
        return 0;
    }

    static SparseDim<T, N> load(const std::string &path, int& errcode) {
        errcode = 0;
        std::ifstream fin(path.c_str(), std::ios::in | std::ios::binary);
        if (!fin) {
            errcode = 1;
        }

        if (!fin.good()) {
            errcode = 2;
        }
        if(errcode > 0){
            std::cout << "error load" << errcode <<std::endl;
            std::array<T, N> dim;
            dim.fill(1);
            return SparseDim<T, N>(dim, T());
        }

        T defval;
        fin.read(reinterpret_cast<char*>(&defval), sizeof(T));

        std::array<T, N> dim;
        int sz;
        fin.read((char*) &sz, sizeof(int));
        assert(sz == N);
        fin.read(reinterpret_cast<char*>(&dim[0]), sz * sizeof(dim[0]));

        SparseDim<T, N> sm(dim, defval);

        loadVOfV(fin, sm.rowlst_lst);
        loadVOfV(fin, sm.collst_lst);
        loadVOfV(fin, sm.vallst_lst);
        fin.close();
        return sm;
    }

    // === FRIEND FUNCTIONS =========================================

    friend bool operator ==(const SparseDim<T, N> & a,
            const SparseDim<T, N> & b) {
        if (a.printallowed)
            std::cout << __func__ << __LINE__ << std::endl;
        if (a._defval == b._defval) {
            return ((a.vallst_lst == b.vallst_lst)
                    && (a.collst_lst == b.collst_lst) // I don't need
                    && (a.rowlst_lst == b.rowlst_lst) // to check _virdimsz
                    && (a._dims_sz == b._dims_sz));   // : Derived value
        }
        if (a._dims_sz != b._dims_sz) {
            return false;
        }
        bool match = true;
        int vgrp = -1;
        for (auto arr : mdrange<N>(a._dims_sz)) {
            if (!arr[N - 2] && !arr[N - 1]) {
                vgrp++;
            }
            if (a._get(arr, vgrp) != b._get(arr, vgrp)) {
                match = false;
                break;
            }
        }
        return match;
    }

    friend bool operator !=(const SparseDim<T, N> & a,
            const SparseDim<T, N> & b) {
        if (a.printallowed)
            std::cout << __func__ << __LINE__ << std::endl;
        return !(a == b);
    }

    friend std::ostream & operator <<(std::ostream & os,
            const SparseDim<T, N> & sd) {
        for (int _virtualidx : patch::xrange(sd._virdimsz)) {
            std::cout << "Struct{" << _virtualidx << "::";
            std::cout << "\nR::";
            for (auto i : sd.rowlst_lst[_virtualidx]) {
                os << i << " ";
            }
            std::cout << "\nC::";
            for (auto i : sd.collst_lst[_virtualidx]) {
                os << i << " ";
            }
            std::cout << "\nV::";
            for (auto i : sd.vallst_lst[_virtualidx]) {
                os << i << " ";
            }
            os << "}\n";
        }
        return os;
    }

protected:
    std::array<int, N> _dims_sz;
    size_t _virdimsz;
    std::vector<std::vector<T> > vallst_lst;
    std::vector<std::vector<int> > rowlst_lst, collst_lst;
    T _defval;
    int printallowed;

    template<typename UniOpType>
    inline SparseDim<T, N> _oper(UniOpType op, bool cacheval) const {
        SparseDim<T, N> result(*this);
        op(result._defval);
        if (cacheval) {
            std::unordered_map<T, T> values;
            values.insert( { this->_defval, result._defval });
            for (auto& vallst : result.vallst_lst) {
                for (T &i : vallst) {
                    if (!values.count(i)) {
                        T j = i;
                        op(j);
                        values.insert( { i, j });
                    }
                    i = values.at(i);
                }
            }
        } else {
            for (auto& vallst : result.vallst_lst) {
                for (T &i : vallst) {
                    op(i);
                }
            }
        }
        return result;
    }

    // === HELPERS / VALIDATORS ==============================================

    void deepCopy(const SparseDim<T, N> & sd) {
        if (printallowed)
            std::cout << __func__ << __LINE__ << std::endl;
        this->_virdimsz = sd._virdimsz;
        this->_defval = sd._defval;
        std::copy(std::begin(sd._dims_sz), std::end(sd._dims_sz),
                std::begin(_dims_sz));
        this->vallst_lst = std::vector<std::vector<T> >(_virdimsz);
        this->collst_lst = std::vector<std::vector<int> >(_virdimsz);
        this->rowlst_lst = std::vector<std::vector<int> >(_virdimsz);
        for (int i : patch::xrange(_virdimsz)) {
            this->vallst_lst[i] = std::vector<T>(sd.vallst_lst[i]);
            this->vallst_lst[i].reserve(10);
            this->collst_lst[i] = std::vector<int>(sd.collst_lst[i]);
            this->collst_lst[i].reserve(10);
            // We make row one size larger in CSR
            this->rowlst_lst[i] = std::vector<int>(sd.rowlst_lst[i]);
        }
    }

    inline size_t validateGetVirtualDim(
            const std::array<int, N> & coord) const {
        if (printallowed > 2)
            std::cout << __func__ << __LINE__ << std::endl;
        size_t pv = 0;
        unsigned mul = 1;
        if ((_virdimsz != vallst_lst.size()) || (_virdimsz != collst_lst.size())
                || (_virdimsz != rowlst_lst.size())) {
            throw InvalidStateException(
                    "Somehow someone changed size of array.");
        }

        size_t sz = _dims_sz.size();
        for (auto av : patch::binditer(coord)) {
            size_t actidx = sz - av.first - 1;
            if (av.second < 0 || av.second >= _dims_sz[actidx]) {
                throw InvalidCoordinatesException("Coordinates out of range.");
            }
            if (av.first > 1) {
                pv = (pv + mul * av.second);
                mul *= _dims_sz[actidx];
            }
        }
        return pv;
    }

    inline T _get(const std::array<int, N> &atloc, int _virtualidx = -1) const {
        // No bound check. It will make it slower
        if (printallowed > 1)
            std::cout << __func__ << __LINE__ << std::endl;
        int currCol;
        if (_virtualidx == -1) {
            _virtualidx = this->validateGetVirtualDim(atloc);
        }
        int row = atloc[N - 2];
        int col = atloc[N - 1];
        auto& rowlst = this->rowlst_lst[_virtualidx];
        auto& collst = this->collst_lst[_virtualidx];
        int pos;
        for (pos = rowlst[row]; pos < rowlst[row + 1]; pos++) {
            currCol = collst[pos];

            if (currCol == col) {
                return this->vallst_lst[_virtualidx][pos];
            } else if (currCol > col) {
                break;
            }
        }
        return _defval;
    }

    inline void _set(const T &val, const std::array<int, N> &atloc,
            int _virtualidx = -1) {
        // No bound check. It will make it slower
        if (printallowed > 1)
            std::cout << __func__ << __LINE__ << std::endl;
        int currCol = -1;
        if (_virtualidx == -1) {
            _virtualidx = this->validateGetVirtualDim(atloc);
        }

        int row = atloc[N - 2];
        int col = atloc[N - 1];
        auto& rowlst = this->rowlst_lst[_virtualidx];
        auto& collst = this->collst_lst[_virtualidx];
        auto& vallst = this->vallst_lst[_virtualidx];
        int pos;
        for (pos = rowlst[row]; pos < rowlst[row + 1]; pos++) {
            currCol = collst[pos];

            if (currCol >= col) {
                break;
            }
        }
        if (currCol != col) {
            if (val == _defval) {
                //its not there and we will not insert defval, so its noop
                return;
            }
            this->insert(_virtualidx, pos, row, col, val);
        } else {
            if (val == _defval) {
                // We are going to add def val, better remove it.
                this->remove(_virtualidx, pos, row);
            } else {
                vallst[pos] = val;
            }
        }
    }

    inline void insert(int _virtualidx, int dataindex, int row, int col,
            T val) {
        if (printallowed > 1)
            std::cout << __func__ << __LINE__ << ":" << val << std::endl;
        auto& vallst = this->vallst_lst[_virtualidx];
        auto& collst = this->collst_lst[_virtualidx];
        if (vallst.size() == 0) {
            vallst.push_back(val);
            collst.push_back(col);

        } else {
            vallst.insert(vallst.begin() + dataindex, val);
            collst.insert(collst.begin() + dataindex, col);
        }

        for (int i : patch::xrange(row + 1, this->_dims_sz[N - 2] + 1)) {
            // We make row one size larger in CSR
            this->rowlst_lst[_virtualidx][i] += 1;
        }
    }

    inline void remove(int _virtualidx, int dataindex, int row) {
        if (printallowed > 1)
            std::cout << __func__ << __LINE__ << std::endl;
        auto& vallst = this->vallst_lst[_virtualidx];
        auto& collst = this->collst_lst[_virtualidx];
        vallst.erase(vallst.begin() + dataindex);
        collst.erase(collst.begin() + dataindex);

        for (int i : patch::xrange(row + 1, this->_dims_sz[N - 2] + 1)) {
            // We make row one size larger in CSR
            this->rowlst_lst[_virtualidx][i] -= 1;
        }
    }

private:
    template<typename T1>
    void dumpVOfV(std::ofstream &fout,
            const std::vector<std::vector<T1>>& vOfv) const {
        auto sz = vOfv.size();
        fout.write(reinterpret_cast<const char*>(&sz), sizeof(int));
        for (auto _v : vOfv) {
            sz = _v.size();
            fout.write(reinterpret_cast<const char*>(&sz), sizeof(int));
            fout.write(reinterpret_cast<const char*>(&_v[0]), sz * sizeof(T1));
        }
    }

    template<typename T1>
    static void loadVOfV(std::ifstream &fin,
            std::vector<std::vector<T1>>& vOfv) {
        int sz;
        fin.read(reinterpret_cast<char*>(&sz), sizeof(int));
        vOfv.resize(sz);
        for (auto& _v : vOfv) {
            fin.read(reinterpret_cast<char*>(&sz), sizeof(int));
            _v.resize(sz);
            fin.read(reinterpret_cast<char*>(&_v[0]), sz * sizeof(T1));
        }
    }

};
