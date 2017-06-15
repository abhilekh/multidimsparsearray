#include <cassert>
#include <ctime>       /* time */
#include "src/sparseDim.tpp"
#include "src/sultiDim.tpp"
#include "src/sparseMatrix.tpp"


using namespace std;

SparseMatrix<int> generateIntMatrix(int rows, int columns) {
    SparseMatrix<int> matrix(rows, columns);

    for (int i =0; i <rows; i++) {
        for (int j = 0; j< columns; j++)) {
            auto val = rand() % 101;
            if (val > 20)
                val = 0;
            matrix.set(val, i, j);
        }
    }
    return matrix;
}

void read_write2() {
    SparseMatrix<int> m1 = generateIntMatrix(500, 510);
    SparseMatrix<int> m2 = generateIntMatrix(510, 500);
    assert(m1 != m2);
    // Write _read test and tests equality
    m1.dump("test1.bin");
    int err;
    SparseMatrix<int> m11 = SparseMatrix<int>::load("test1.bin", err);
    assert(m11 == m1);

    m2.dump("test2.bin");
    SparseMatrix<int> m12 = SparseMatrix<int>::load("test2.bin", err);
    assert(m12 == m2);
}

struct foo7 {
    void operator()(int& i) {
        i *= 7;
    }
};

struct isOdd {
    void operator()(int &i) {
        i %= 2;
    }
};

void basic_operation2() {
    SparseMatrix<int> m1 = generateIntMatrix(500, 510);
    SparseMatrix<int> m2 = generateIntMatrix(500, 510);

    assert(m1 != m2);
    assert((m1 + m2) == (m2 + m1));
    assert((m1 * 2) == (m1 + m1));
    assert((m1 + (m2 * -1)) == (m1 - m2));
    assert(((m1 - m2) * -1) == (m2 - m1));
    assert(((m1 * 5) / 5) == m1);

    assert(m2.oper(foo7()) == m2 * 7);
    assert(m2.oper(isOdd()) == m2 % 2);
}

template<int N>
SparseDim<int, N> generateIntNDimArr(array<int, N> arr) {
    SparseDim<int, N> sd(arr);
    for (std::array<int, N> _arr : mdrange<N>(arr)) {
        int val = rand() % 101;
        if (val > 20)
            val = 0;
        sd.set(val, _arr);
    }
    return sd;
}

template<int N>
array<int, N> getRandomIndex(array<int, N> sz) {
    for (int i =0; i<N; i++) {
        sz[i] = (sz[i] != 1) ? rand() % sz[i] : 0;
    }
    return sz;
}

void set_getND() {
    {
        array<int, 5> arr1 = { { 40, 50, 30, 40, 70 } };
        SparseDim<int, 5> sd(arr1);
        for (int i=0 ; i< 500; i++)) {
            UNUSED(i);
            auto dim = getRandomIndex<5>(arr1);
            int val = rand() % 100;
            sd.set(val, dim);
            int orgval2 = sd.get(dim);
            assert(orgval2 == val);
        }
    }
    std::cout << "Case 1: SetGet done " << endl;

    {
        array<int, 10> arr2 = { { 4, 5, 30, 1, 9, 10, 4, 4, 12, 8 } };
        SparseDim<int, 10> sd(arr2);
        for (int i =0; i< 200; i++) {
            UNUSED(i);
            auto dim = getRandomIndex<10>(arr2);
            int val = rand() % 100;
            sd.set(val, dim);
            int orgval2 = sd.get(dim);
            assert(orgval2 == val);
        }
    }
    std::cout << "Case 2: SetGet done " << endl;

    {
        array<int, 6> arr3 = { { 8, 10, 5, 3, 3, 6 } };
        SparseDim<int, arr3.size()> sd = generateIntNDimArr<6>(arr3);
        for (int i =0; i< 200; i++) {
            UNUSED(i);
            auto dim = getRandomIndex<6>(arr3);
            int val = rand() % 100;
            sd.set(val, dim);
            int orgval2 = sd.get(dim);
            assert(orgval2 == val);
        }
    }
    std::cout << "Case 3: SetGet done " << endl;
}

void basic_operationND() {
    {
        array<int, 3> arr = { { 2, 4, 4 } };
        MultiDim<int, 3> m1 = generateIntNDimArr<3>(arr);
        MultiDim<int, 3> m2 = generateIntNDimArr<3>(arr);
        assert(m1 != m2);
        assert((m1 + m2) == (m2 + m1));
        assert((m1 * 2) == (m1 + m1));
        assert((m1 + (m2 * -1)) == (m1 - m2));
        assert(((m1 - m2) * -1) == (m2 - m1));
        assert(((m1 * 5) / 5) == m1);
        assert(m2.oper([](int &i) {i*=7;}) == m2 * 7);
        assert(m2.oper([](int &i) {i%=2;}) == m2 % 2);
    }
    cout << "Case 1: Operations Done "<< endl;

    {
        array<int, 5> arr = { { 20, 3, 7, 8, 9 } };
        MultiDim<int, 5> m1 = generateIntNDimArr<5>(arr);
        MultiDim<int, 5> m2 = generateIntNDimArr<5>(arr);
        assert(m1 != m2);
        assert((m1 + m2) == (m2 + m1));
        assert((m1 * 2) == (m1 + m1));
        assert((m1 + (m2 * -1)) == (m1 - m2));
        assert(((m1 - m2) * -1) == (m2 - m1));
        assert(((m1 * 5) / 5) == m1);
        assert(m2.oper([](int &i) {i*=7;}) == m2 * 7);
        assert(m2.oper([](int &i) {i%=2;}) == m2 % 2);
    }
    cout << "Case 2: Operations Done " << endl;
}

void read_writeND() {
    array<int, 5> arr = { { 20, 3, 7, 8, 9 } };
    SparseDim<int, 5> m1 = generateIntNDimArr<5>(arr);
    SparseDim<int, 5> m2 = generateIntNDimArr<5>(arr);
    assert(m1 != m2);
    // Write _read test and tests equality
    m1.dump("_test1.bin");
    int err;
    SparseDim<int, 5> m11 = SparseDim<int, 5>::load("_test1.bin", err);
    assert(m11 == m1);

    m2.dump("_test2.bin");
    SparseDim<int, 5> m12 = SparseDim<int, 5>::load("_test2.bin", err);
    assert(m12 == m2);
}

int main() {
    srand(time(NULL));
    read_write2();
    basic_operation2();
    set_getND();
    basic_operationND();
    read_writeND();
    return 0;
}