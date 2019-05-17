#include "Matrix.h"
#include "math.h"

Matrix::Matrix(unsigned int nRows,
               unsigned int nCols,
               const LinearAlgebra& linalg){

    // Calculate decomposed row block size
    float ratio = float(nRows) / float(linalg.size());
    unsigned int splitIndex = nRows - linalg.size() * floor( ratio );
    this->nLocalRows = (linalg.rank() < splitIndex) ? ceil( ratio ) : floor( ratio );

    this->nLocalColumns = nCols;
    this->nGlobalRows = nRows;
    this->nGlobalColumns = nCols;

    // Calculating global index range for local rows
    if (linalg.rank() < splitIndex){
        this->globalRowIndexRange.begin = linalg.rank() * this->nLocalRows;
        this->globalRowIndexRange.end = (linalg.rank() + 1) * this->nLocalRows;
    } else {
        unsigned int offset = splitIndex * (this->nLocalRows + 1);
        this->globalRowIndexRange.begin = (linalg.rank() - splitIndex) * this->nLocalRows + offset;
        this->globalRowIndexRange.end = (linalg.rank() - splitIndex +1) * this->nLocalRows + offset;
    }

    this->linalg = &linalg;
    this->matrixStorage = new double*[this->nLocalRows];
    for (unsigned int row = 0; row < this->nLocalRows; row++){
        this->matrixStorage[row] = new double[this->nLocalColumns];
        // Initialize entries to zero
        for (unsigned int col = 0; col < this->nLocalColumns; col++)
            this->matrixStorage[row][col] = 0.0;
    }
}

Matrix::~Matrix(){
    for (unsigned int row = 0; row < this->nLocalRows; row++)
        delete[] this->matrixStorage[row];
    delete[] this->matrixStorage;
}

std::vector<unsigned int> Matrix::getPartitionSize() const {
    std::vector<unsigned int> result(this->linalg->size());
    MPI_Allgather(&this->nLocalRows, 1, MPI_UNSIGNED, result.data(), 1, MPI_UNSIGNED, MPI_COMM_WORLD);
    return result;
}

void Matrix::setValues(const double &x) {
    for (unsigned int row = 0; row < this->nLocalRows; row++)
        for (unsigned int col = 0; col < this->nLocalColumns; col++)
            this->matrixStorage[row][col] = x;
}

void Matrix::setValues(const std::vector<std::vector<double>>& x){
    Nucleus_ASSERT_EQ(x.size(), this->nGlobalRows)
    for (unsigned int loc_row = 0, glob_row = this->globalRowIndexRange.begin; loc_row < this->nLocalRows; loc_row++, glob_row++){
        Nucleus_ASSERT_EQ(x[glob_row].size(), this->nGlobalColumns)
        for (unsigned int col = 0; col < this->nLocalColumns; col++)
            this->matrixStorage[loc_row][col] = x[glob_row][col];
    }
}

void Matrix::zeros() {
    for (unsigned int row = 0; row < this->nLocalRows; row++)
        for (unsigned int col = 0; col < this->nLocalColumns; col++)
            this->matrixStorage[row][col] = 0.0;
}

double Matrix::frobeniusNorm() const {
    double local_result = 0.0;
    for (unsigned int row = 0; row < this->nLocalRows; row++)
        for (unsigned int col = 0; col < this->nLocalColumns; col++)
            local_result += this->matrixStorage[row][col] * this->matrixStorage[row][col];

    double global_result = 0.0;
    MPI_Allreduce(&local_result, &global_result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return sqrt(global_result);
}
