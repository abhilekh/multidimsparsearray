/**
 * This file is part of the MultiDimSparseArray library
 *
 * @license  BSD-3
 * @author   Abhilekh Agarwal
 */

#ifndef __SPARSEMATRIX_H__
#define __SPARSEMATRIX_H__

#include <algorithm>
#include <map>
#include <vector>
#include "MultiDim.tpp"

template<typename T>
class SparseMatrix: public MultiDim<T, 2U> {
public:
    // === CREATION ==============================================
    SparseMatrix(int rows, int columns, T defval = T()) :
        MultiDim<T, 2U>( { { rows, columns } }, defval) { // Constructor 1
        // this->printallowed = 1;
    }

    SparseMatrix(const MultiDim<T, 2U> & matrix) :
        MultiDim<T, 2U>(matrix) {  // Constructor 2
        //this->printallowed = 1;
    }

    SparseMatrix(const SparseDim<T, 2U> & matrix) :
        MultiDim<T, 2U>(matrix) {  // Constructor 3
        //this->printallowed = 1;
    }

    SparseMatrix<T> & operator =(const SparseMatrix<T> & matrix) {
        if (&matrix != this) {
            this->deepCopy(matrix);

        }

        return *this;
    }

    SparseMatrix<T> & operator =(const MultiDim<T, 2U> & matrix) {
        if (&matrix != this) {
            this->deepCopy(matrix);
        }
        return *this;
    }

    SparseMatrix<T> & operator =(const SparseDim<T, 2U> & matrix) {
        if (&matrix != this) {
            this->deepCopy(matrix);
        }
        return *this;
    }

    // === GETTERS / SETTERS ==============================================
    int getRowCount(void) const {
        return this->_dims_sz[0];
    }
    int getColumnCount(void) const {
        return this->_dims_sz[1];
    }

    // === VALUES ==============================================
    using MultiDim<T, 2U>::get;
    T get(int row, int col) const {
        return this->get( { { row, col } });
    }

    using MultiDim<T, 2U>::set;
    void set(const T &val, int row, int col) {
        this->set(val, { { row, col } });
    }

    // === OPERATIONS ==============================================

    using MultiDim<T, 2U>::operator *;
    std::vector<T> operator *(const std::vector<T> & x) const {

        if (this->_ncols != (int) x.size()) {
            throw InvalidDimensionsException(
                    "Cannot multiply: Matrix column count and vector size don't match.");
        }
        std::vector<T> result(this->_nrows, this->_defval);
        if (this->vallst.size() == 0) { // only if any value set
            return result;
        }
        for (int i = 0; i < this->_nrows; i++) {
            T sum = T();
            for (int j = this->rowlst.at(i); j < this->rowlst.at(i + 1); j++) {
                sum = sum
                        + this->vallst.at(j - 1)
                        * x[this->collst.at(j - 1) - 1];
            }

            result[i] = sum;
        }
        return result;
    }

    SparseMatrix<T> operator *(const SparseMatrix<T> & m) const {
        if (this->_ncols != m._nrows) {
            throw InvalidDimensionsException(
                    "Cannot multiply: Left matrix column count and right matrix row count don't match.");
        }

        SparseMatrix<T> result(this->_nrows, m._ncols);
        T a;

        // TODO: more efficient?
        // @see http://www.math.tamu.edu/~srobertp/Courses/Math639_2014_Sp/CRSDescription/CRSStuff.pdf

        for (int i = 1; i <= this->_nrows; i++) {
            for (int j = 1; j <= m._ncols; j++) {
                a = T();

                for (int k = 1; k <= this->_ncols; k++) {
                    a = a + this->get(i, k) * m.get(k, j);
                }

                result.set(a, i, j);
            }
        }

        return result;
    }

protected:

    // === HELPERS / VALIDATORS ==============================================
    void deepCopy(const SparseMatrix<T> & m);
    void validateCoordinates(int row, int col) const;

};

#endif
