#include "Matrix_BlockRowPartition.h"
#include "math.h"

Matrix_BlockRowPartition::Matrix_BlockRowPartition(unsigned int nRows,
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

Matrix_BlockRowPartition::Matrix_BlockRowPartition(const Matrix_BlockRowPartition &other){
    this->nLocalRows = other.nLocalRows;
    this->nLocalColumns = other.nLocalColumns;
    this->nGlobalRows = other.nGlobalRows;
    this->nGlobalColumns = other.nGlobalColumns;

    this->globalRowIndexRange = other.globalRowIndexRange;
    this->linalg = other.linalg;

    this->matrixStorage = new double*[this->nLocalRows];
    for (unsigned int row = 0; row < this->nLocalRows; row++){
        this->matrixStorage[row] = new double[this->nLocalColumns];
        for (unsigned int col = 0; col < this->nLocalColumns; col++)
            this->matrixStorage[row][col] = other.matrixStorage[row][col];
    }
}

Matrix_BlockRowPartition& Matrix_BlockRowPartition::operator=(const Matrix_BlockRowPartition &other){
    if (this != &other){
        this->nLocalRows = other.nLocalRows;
        this->nLocalColumns = other.nLocalColumns;
        this->nGlobalRows = other.nGlobalRows;
        this->nGlobalColumns = other.nGlobalColumns;

        this->globalRowIndexRange = other.globalRowIndexRange;
        this->linalg = other.linalg;

        this->matrixStorage = new double*[this->nLocalRows];
        for (unsigned int row = 0; row < this->nLocalRows; row++){
            this->matrixStorage[row] = new double[this->nLocalColumns];
            for (unsigned int col = 0; col < this->nLocalColumns; col++)
                this->matrixStorage[row][col] = other.matrixStorage[row][col];
        }
    }
    return *this;
}

std::unique_ptr<AbstractMatrix> Matrix_BlockRowPartition::clone() const {
    return std::make_unique<Matrix_BlockRowPartition>(Matrix_BlockRowPartition(*this));
}

Matrix_BlockRowPartition::~Matrix_BlockRowPartition(){
    for (unsigned int row = 0; row < this->nLocalRows; row++)
        delete[] this->matrixStorage[row];
    delete[] this->matrixStorage;
}

std::vector<unsigned int> Matrix_BlockRowPartition::getPartitionSize() const {
    std::vector<unsigned int> result(this->linalg->size());
    MPI_Allgather(&this->nLocalRows, 1, MPI_UNSIGNED, result.data(), 1, MPI_UNSIGNED, MPI_COMM_WORLD);
    return result;
}

void Matrix_BlockRowPartition::setValues(const double &x) {
    for (unsigned int row = 0; row < this->nLocalRows; row++)
        for (unsigned int col = 0; col < this->nLocalColumns; col++)
            this->matrixStorage[row][col] = x;
}

void Matrix_BlockRowPartition::setValues(const std::vector<std::vector<double>>& x){
    Nucleus_ASSERT_EQ(x.size(), this->nGlobalRows)
    for (unsigned int loc_row = 0, glob_row = this->globalRowIndexRange.begin; loc_row < this->nLocalRows; loc_row++, glob_row++){
        Nucleus_ASSERT_EQ(x[glob_row].size(), this->nGlobalColumns)
        for (unsigned int col = 0; col < this->nLocalColumns; col++)
            this->matrixStorage[loc_row][col] = x[glob_row][col];
    }
}

void Matrix_BlockRowPartition::zeros() {
    for (unsigned int row = 0; row < this->nLocalRows; row++)
        for (unsigned int col = 0; col < this->nLocalColumns; col++)
            this->matrixStorage[row][col] = 0.0;
}

unsigned int Matrix_BlockRowPartition::findRankWithIndex(const unsigned int index) const {
    float ratio = float(this->nGlobalRows) / float(this->linalg->size());
    unsigned int splitRank = this->nGlobalRows - this->linalg->size() * floor( ratio );
    unsigned int splitIndex = splitRank * ceil(ratio);
    unsigned int ownerProc;
    if (index < splitIndex)
        ownerProc = floor( index / ceil(ratio) );
    else {
        unsigned int offset = splitRank * ceil(ratio);
        ownerProc = floor( (index - offset) / floor(ratio) ) + splitRank;
    }

    return ownerProc;
}

double Matrix_BlockRowPartition::getValue(const unsigned int row, const unsigned int col) const {
    Nucleus_ASSERT_LT(row, this->nGlobalRows)
    Nucleus_ASSERT_LT(col, this->nGlobalColumns)

    unsigned int indexOnRank = this->findRankWithIndex(row);

    double result;
    if (indexOnRank == this->linalg->rank())
        result = this->matrixStorage[row - this->globalRowIndexRange.begin][col];

    MPI_Bcast(&result, 1, MPI_DOUBLE, indexOnRank, MPI_COMM_WORLD);

    return result;
}

double Matrix_BlockRowPartition::frobeniusNorm() const {
    double local_result = 0.0;
    for (unsigned int row = 0; row < this->nLocalRows; row++)
        for (unsigned int col = 0; col < this->nLocalColumns; col++)
            local_result += this->matrixStorage[row][col] * this->matrixStorage[row][col];

    double global_result = 0.0;
    MPI_Allreduce(&local_result, &global_result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return sqrt(global_result);
}

std::unique_ptr<AbstractVector> Matrix_BlockRowPartition::mult(const AbstractVector &other) const {
    Nucleus_ASSERT_EQ(this->nGlobalColumns, other.size())
    std::unique_ptr<AbstractVector> result = std::make_unique<Vector_BlockPartition>(this->nGlobalRows, *this->linalg);

    std::vector<double> global_values(other.size());
    other.getGlobalValues(global_values);

    std::vector<double> result_values(result->size());
    for (unsigned int row = 0, global_row = this->globalRowIndexRange.begin; row < this->nLocalRows; row++, global_row++){
        for (unsigned int col = 0; col < this->nLocalColumns; col++){
            result_values[global_row] += this->matrixStorage[row][col] * global_values[col];
        }
    }

    result->setValues(result_values);
  
    return result;
}

void Matrix_BlockRowPartition::setRowValues(const unsigned int row, const std::vector<double> &values){
    Nucleus_ASSERT_EQ(values.size(), this->nGlobalColumns)
    Nucleus_ASSERT_LT(row, this->nGlobalRows)

    unsigned int indexOnRank = this->findRankWithIndex(row);
    if (indexOnRank == this->linalg->rank())
        for (unsigned int i = 0; i < this->nLocalColumns; i++)
            this->matrixStorage[row-this->globalRowIndexRange.begin][i] = values[i];
}

void Matrix_BlockRowPartition::setLocalRowValues(const unsigned int row, const std::vector<double> &values){
    Nucleus_ASSERT_EQ(values.size(), this->nLocalColumns)
    Nucleus_ASSERT_LT(row, this->nLocalRows)

    for (unsigned int i = 0; i < this->nLocalColumns; i++)
        this->matrixStorage[row][i] = values[i];
}

std::vector<double> Matrix_BlockRowPartition::getRowValues(const unsigned int row) const {
    Nucleus_ASSERT_LT(row, this->nGlobalRows)
    std::vector<double> result(this->nLocalColumns);

    unsigned int indexOnRank = this->findRankWithIndex(row);

    if (this->linalg->rank() == indexOnRank)
        for (unsigned int i = 0; i < this->nLocalColumns; i++)
            result[i] = this->matrixStorage[row - this->globalRowIndexRange.begin][i];

    MPI_Bcast(result.data(), this->nLocalColumns, MPI_DOUBLE, indexOnRank, MPI_COMM_WORLD);

    return result;
}

double& Matrix_BlockRowPartition::getLocalRowValues(const unsigned int row) const {
    Nucleus_ASSERT_LT(row, this->nLocalRows)

    return *this->matrixStorage[row];
    // std::vector<double> result(this->nLocalColumns);

    // for (unsigned int i = 0; i < this->nLocalColumns; i++)
    //     result[i] = this->matrixStorage[row][i];

    // return result;
}
