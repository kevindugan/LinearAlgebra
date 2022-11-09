#include "SparseMatrix_BlockRowPartition.h"
#include "math.h"
#include <iomanip>
#include <algorithm>

SparseMatrix_BlockRowPartition::SparseMatrix_BlockRowPartition(unsigned int nRows,
                                                               unsigned int nCols,
                                                               const Parallel& linalg,
                                                               const double tol){

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
    this->zeroTol = tol;
    this->values.resize(0);
    this->columnIndices.resize(0);

    this->rowPtr.resize(this->nLocalRows+1);
}

SparseMatrix_BlockRowPartition::SparseMatrix_BlockRowPartition(const SparseMatrix_BlockRowPartition &other){
    this->nLocalRows = other.nLocalRows;
    this->nLocalColumns = other.nLocalColumns;
    this->nGlobalRows = other.nGlobalRows;
    this->nGlobalColumns = other.nGlobalColumns;

    this->globalRowIndexRange = other.globalRowIndexRange;
    this->linalg = other.linalg;
    this->zeroTol = other.zeroTol;

    this->values = other.values;
    this->columnIndices = other.columnIndices;
    this->rowPtr = other.rowPtr;
}

SparseMatrix_BlockRowPartition& SparseMatrix_BlockRowPartition::operator=(const SparseMatrix_BlockRowPartition &other){
    if (this != &other){
        this->nLocalRows = other.nLocalRows;
        this->nLocalColumns = other.nLocalColumns;
        this->nGlobalRows = other.nGlobalRows;
        this->nGlobalColumns = other.nGlobalColumns;

        this->globalRowIndexRange = other.globalRowIndexRange;
        this->linalg = other.linalg;
        this->zeroTol = other.zeroTol;

        this->values = other.values;
        this->columnIndices = other.columnIndices;
        this->rowPtr = other.rowPtr;
    }
    return *this;
}

SparseMatrix_BlockRowPartition::~SparseMatrix_BlockRowPartition(){}

std::unique_ptr<AbstractMatrix> SparseMatrix_BlockRowPartition::clone() const {
    return std::make_unique<SparseMatrix_BlockRowPartition>(SparseMatrix_BlockRowPartition(*this));
}

std::vector<unsigned int> SparseMatrix_BlockRowPartition::getPartitionSize() const {
    std::vector<unsigned int> result(this->linalg->size());
    MPI_Allgather(&this->nLocalRows, 1, MPI_UNSIGNED, result.data(), 1, MPI_UNSIGNED, MPI_COMM_WORLD);
    return result;
}

std::vector<double> SparseMatrix_BlockRowPartition::getRowValues(unsigned int row) const {
    Nucleus_ASSERT_LT(row, this->nGlobalRows)
    std::vector<double> result(this->nLocalColumns);

    unsigned int indexOnRank = this->findRankWithIndex(row);

    if (this->linalg->rank() == indexOnRank){
        unsigned int localRow = row - this->globalRowIndexRange.begin;
        auto col_begin = this->columnIndices.begin() + this->rowPtr[localRow];
        auto col_end = this->columnIndices.begin() + this->rowPtr[localRow+1];
        auto val_begin = this->values.begin() + this->rowPtr[localRow];
        auto val_end = this->values.begin() + this->rowPtr[localRow+1];
        auto col = col_begin;
        auto val = val_begin;
        for (; col != col_end; col++, val++){
            result[*col] = *val;
        }
    }

    MPI_Bcast(result.data(), this->nLocalColumns, MPI_DOUBLE, indexOnRank, MPI_COMM_WORLD);

    return result;
}

void SparseMatrix_BlockRowPartition::setValues(const double &x){
    // Reset storage
    this->values.clear();
    this->columnIndices.clear();
    this->rowPtr.clear();

    // Reserve space, set values the same
    unsigned int nsquare = this->nLocalColumns * this->nLocalRows;
    this->values.resize(nsquare, x);
    this->columnIndices.reserve(nsquare);
    this->rowPtr.reserve(this->nLocalRows+1);

    // Set indices
    for (unsigned int it = 0; it < nsquare; it++)
        this->columnIndices.push_back( it % this->nLocalColumns );
    for (unsigned int it = 0; it < this->nLocalRows+1; it++)
        this->rowPtr.push_back(it*this->nLocalColumns);
}

void SparseMatrix_BlockRowPartition::setValues(const std::vector<std::vector<double>> &x){
    // Reset storage
    this->values.clear();
    this->columnIndices.clear();
    this->rowPtr.clear();

    // Reserve space
    unsigned int nsquare = this->nLocalColumns * this->nLocalRows;
    this->values.reserve(nsquare);
    this->columnIndices.reserve(nsquare);
    this->rowPtr.reserve(this->nLocalRows+1);

    // Set Values
    this->rowPtr.push_back(0);
    Nucleus_ASSERT_EQ(x.size(), this->nGlobalRows)
    for (unsigned int loc_row = 0, glob_row = this->globalRowIndexRange.begin; loc_row < this->nLocalRows; loc_row++, glob_row++){
        Nucleus_ASSERT_EQ(x[glob_row].size(), this->nGlobalColumns)
        for (unsigned int col = 0; col < this->nLocalColumns; col++)
            if ( abs(x[glob_row][col]) > this->zeroTol){
                this->values.push_back( x[glob_row][col] );
                this->columnIndices.push_back( col );
            }
        this->rowPtr.push_back(this->columnIndices.size());
    }

    // Shrink Vectors
    this->values.shrink_to_fit();
    this->columnIndices.shrink_to_fit();
}

void SparseMatrix_BlockRowPartition::setRowValues(const unsigned int row,
                                                  const std::vector<double> &x){
    Nucleus_ASSERT_EQ(x.size(), this->nGlobalColumns)
    Nucleus_ASSERT_LT(row, this->nGlobalRows)

    unsigned int indexOnRank = this->findRankWithIndex(row);
    if (indexOnRank == this->linalg->rank()){
        unsigned int localRow = row - this->globalRowIndexRange.begin;
        auto col_begin = this->columnIndices.begin() + this->rowPtr[localRow];
        auto col_end = this->columnIndices.begin() + this->rowPtr[localRow+1];
        auto val_begin = this->values.begin() + this->rowPtr[localRow];
        auto val_end = this->values.begin() + this->rowPtr[localRow+1];

        // New Storage Vector
        std::vector<unsigned int> newCol, newRow(this->rowPtr.size());
        std::vector<double> newVal;
        // Copy upto current row
        newCol.resize(col_begin-this->columnIndices.begin());
        std::copy(this->columnIndices.begin(), col_begin, newCol.begin());
        newVal.resize(val_begin-this->values.begin());
        std::copy(this->values.begin(), val_begin, newVal.begin());
        std::copy(this->rowPtr.begin(), this->rowPtr.begin()+localRow, newRow.begin());

        // Set Current Row
        newRow[localRow] = newCol.size();
        for (unsigned int it = 0; it < x.size(); it++)
            if ( abs(x[it]) > this->zeroTol){
                newVal.push_back( x[it] );
                newCol.push_back( it );
            }
        newRow[localRow+1] = newCol.size();

        // Fill later row values
        newVal.resize(this->values.end()-val_end+newVal.size());
        newCol.resize(this->columnIndices.end()-col_end+newCol.size());
        std::copy(val_end, this->values.end(), newVal.begin()+newRow[localRow+1]);
        std::copy(col_end, this->columnIndices.end(), newCol.begin()+newRow[localRow+1]);

        // Modify row pointer
        int change =  newRow[localRow+1] - this->rowPtr[localRow+1];
        for (unsigned int it = localRow+2; it < newRow.size(); it++)
            newRow[it] = this->rowPtr[it] + change;

        // save new values
        this->columnIndices = newCol;
        this->values = newVal;
        this->rowPtr = newRow;
        this->rowPtr.resize(this->nLocalRows+1, newRow.back());
    }
}

void SparseMatrix_BlockRowPartition::zeros(){
    // Reset storage
    this->values.clear();
    this->columnIndices.clear();
    this->rowPtr.clear();

    // Set appropriate row pointer [0, 0, 0, ...]
    this->rowPtr.resize(this->nLocalRows+1);
}

unsigned int SparseMatrix_BlockRowPartition::findRankWithIndex(const unsigned int index) const {
    double ratio = float(this->nGlobalRows) / float(this->linalg->size());
    unsigned int leftPerRank = ceil(ratio);
    unsigned int rightPerRank = floor(ratio);
    unsigned int splitIndex = (this->nGlobalRows % this->linalg->size()) * leftPerRank;
    if (index < splitIndex)
        return floor( float(index)/float(leftPerRank) );
    else
        return floor( float(index-splitIndex)/float(rightPerRank) ) + (this->nGlobalRows % this->linalg->size());
}

double SparseMatrix_BlockRowPartition::getValue(const unsigned int row, const unsigned int col) const {
    Nucleus_ASSERT_LT(row, this->nGlobalRows)
    Nucleus_ASSERT_LT(col, this->nGlobalColumns)

    unsigned int indexOnRank = this->findRankWithIndex(row);

    double result;
    if (indexOnRank == this->linalg->rank()){
        unsigned int localRow = row - this->globalRowIndexRange.begin;
        auto col_beg = this->columnIndices.begin()+this->rowPtr[localRow];
        auto col_end = this->columnIndices.begin()+this->rowPtr[localRow+1];
        auto found = std::find(col_beg, col_end, col);
        result = found != col_end ? this->values[found-this->columnIndices.begin()] : 0.0;
    }

    MPI_Bcast(&result, 1, MPI_DOUBLE, indexOnRank, MPI_COMM_WORLD);

    return result;

}

double SparseMatrix_BlockRowPartition::frobeniusNorm() const {
    double local_result = 0.0;
    for (const auto e : this->values)
        local_result += e*e;

    double global_result = 0.0;
    MPI_Allreduce(&local_result, &global_result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return sqrt(global_result);
}

std::unique_ptr<AbstractVector> SparseMatrix_BlockRowPartition::mult(const AbstractVector &other) const {
    Nucleus_ASSERT_EQ(this->nGlobalColumns, other.size())
    std::unique_ptr<AbstractVector> result = std::make_unique<Vector_BlockPartition>(this->nGlobalRows, *this->linalg);

    std::vector<double> global_values(other.size());
    other.getGlobalValues(global_values);

    std::vector<double> result_values(result->size());
    for (unsigned int row_it = 0, global_row = this->globalRowIndexRange.begin; row_it < this->rowPtr.size()-1; row_it++, global_row++)
        for (unsigned int col_it = this->rowPtr[row_it]; col_it < this->rowPtr[row_it+1]; col_it++)
            result_values[global_row] += this->values[col_it] * global_values[this->columnIndices[col_it]];

    result->setValues(result_values);
    return result;
}

void SparseMatrix_BlockRowPartition::setLocalRowValues(const unsigned int row, const std::vector<double> &x){
    Nucleus_ASSERT_EQ(x.size(), this->nLocalColumns)
    Nucleus_ASSERT_LT(row, this->nLocalRows)

    unsigned int indexOnRank = this->findRankWithIndex(row);
    if (indexOnRank == this->linalg->rank()){
        auto col_begin = this->columnIndices.begin() + this->rowPtr[row];
        auto col_end = this->columnIndices.begin() + this->rowPtr[row+1];
        auto val_begin = this->values.begin() + this->rowPtr[row];
        auto val_end = this->values.begin() + this->rowPtr[row+1];

        // New Storage Vector
        std::vector<unsigned int> newCol, newRow(this->rowPtr.size());
        std::vector<double> newVal;
        // Copy upto current row
        newCol.resize(col_begin-this->columnIndices.begin());
        std::copy(this->columnIndices.begin(), col_begin, newCol.begin());
        newVal.resize(val_begin-this->values.begin());
        std::copy(this->values.begin(), val_begin, newVal.begin());
        std::copy(this->rowPtr.begin(), this->rowPtr.begin()+row, newRow.begin());

        // Set Current Row
        newRow[row] = newCol.size();
        for (unsigned int it = 0; it < x.size(); it++)
            if ( abs(x[it]) > this->zeroTol){
                newVal.push_back( x[it] );
                newCol.push_back( it );
            }
        newRow[row+1] = newCol.size();

        // Fill later row values
        newVal.resize(this->values.end()-val_end+newVal.size());
        newCol.resize(this->columnIndices.end()-col_end+newCol.size());
        std::copy(val_end, this->values.end(), newVal.begin()+newRow[row+1]);
        std::copy(col_end, this->columnIndices.end(), newCol.begin()+newRow[row+1]);

        // Modify row pointer
        int change =  newRow[row+1] - this->rowPtr[row+1];
        for (unsigned int it = row+2; it < newRow.size(); it++)
            newRow[it] = this->rowPtr[it] + change;

        // save new values
        this->columnIndices = newCol;
        this->values = newVal;
        this->rowPtr = newRow;
        this->rowPtr.resize(this->nLocalRows+1, newRow.back());
    }
}

std::vector<double> SparseMatrix_BlockRowPartition::getLocalRowValues(const unsigned int row) const {
    Nucleus_ASSERT_LT(row, this->nLocalRows)
    std::vector<double> result(this->nLocalColumns);

    auto col_begin = this->columnIndices.begin() + this->rowPtr[row];
    auto col_end = this->columnIndices.begin() + this->rowPtr[row+1];
    auto val_begin = this->values.begin() + this->rowPtr[row];
    auto val_end = this->values.begin() + this->rowPtr[row+1];
    auto col = col_begin;
    auto val = val_begin;
    for (; col != col_end; col++, val++){
        result[*col] = *val;
    }

    return result;
}

void SparseMatrix_BlockRowPartition::print(std::ostream &out) const {
    std::vector<size_t> storageSize(3*this->linalg->size());
    std::vector<size_t> localStorageSize = {this->columnIndices.size(), this->rowPtr.size(), this->nLocalColumns};
    MPI_Gather(localStorageSize.data(), 3, MPI_UNSIGNED_LONG, storageSize.data(), 3, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    if (this->linalg->rank() == 0){
        // Lambda for printing rows
        auto printLocal = [&out](const std::vector<double> &vals,
                                 const std::vector<unsigned int> &col,
                                 const std::vector<unsigned int> &rows,
                                 const unsigned int &nCols)
        {
            for (unsigned int row_it = 0; row_it < rows.size()-1; row_it++){
                unsigned int printCol = 0, crs_col = rows[row_it];
                // Pad Start
                while (crs_col < col.size() && printCol < col[crs_col]){
                    out << std::string(11, ' ') << " "; printCol++;
                }
                // Interior values
                while ( crs_col < rows[row_it+1])
                    if (printCol == col[crs_col]){
                        out << std::setprecision(4) << std::setw(11) << std::scientific << vals[crs_col] << " ";
                        crs_col++;
                        printCol++;
                    } else {
                        out << std::string(11, ' ') << " ";
                        printCol++;
                    }
                while( printCol < nCols){
                    out << std::string(11, ' ') << " "; printCol++;
                }
                out <<  std::endl;
            }
        };

        // Print Local values
        std::cout << std::string(this->nLocalColumns*12, '=') << std::endl;
        printLocal(this->values, this->columnIndices, this->rowPtr, this->nLocalColumns);
        for (unsigned int proc = 1; proc < this->linalg->size(); proc++){
            std::cout << std::string(this->nLocalColumns*12, '-') << std::endl;
            std::vector<double> otherVal(storageSize[3*proc]);
            std::vector<unsigned int> otherCol(storageSize[3*proc]),
                                      otherRow(storageSize[3*proc+1]);
            MPI_Recv(otherVal.data(), storageSize[3*proc], MPI_DOUBLE, proc, 1, MPI_COMM_WORLD, nullptr);
            MPI_Recv(otherCol.data(), storageSize[3*proc], MPI_UNSIGNED, proc, 2, MPI_COMM_WORLD, nullptr);
            MPI_Recv(otherRow.data(), storageSize[3*proc+1], MPI_UNSIGNED, proc, 3, MPI_COMM_WORLD, nullptr);
            printLocal(otherVal, otherCol, otherRow, storageSize[3*proc+2]);
        }
        std::cout << std::string(this->nLocalColumns*12, '=') << std::endl;
    } else {
        MPI_Send(this->values.data(), this->values.size(), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
        MPI_Send(this->columnIndices.data(), this->columnIndices.size(), MPI_UNSIGNED, 0, 2, MPI_COMM_WORLD);
        MPI_Send(this->rowPtr.data(), this->rowPtr.size(), MPI_UNSIGNED, 0, 3, MPI_COMM_WORLD);
    }
}

void SparseMatrix_BlockRowPartition::printSparsity(std::ostream &out) const {
    std::vector<size_t> storageSize(3*this->linalg->size());
    std::vector<size_t> localStorageSize = {this->columnIndices.size(), this->rowPtr.size(), this->nLocalColumns};
    MPI_Gather(localStorageSize.data(), 3, MPI_UNSIGNED_LONG, storageSize.data(), 3, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    if (this->linalg->rank() == 0){
        // Lambda for printing rows
        auto printLocal = [&out](const std::vector<unsigned int> &col,
                                 const std::vector<unsigned int> &rows,
                                 const unsigned int &nCols)
        {
            for (unsigned int row_it = 0; row_it < rows.size()-1; row_it++){
                unsigned int printCol = 0, crs_col = rows[row_it];
                // Pad Start
                while (crs_col < col.size() && printCol < col[crs_col]){
                    out <<  " "; printCol++;
                }
                // Interior values
                while ( crs_col < rows[row_it+1])
                    if (printCol == col[crs_col]){
                        out << "X";
                        crs_col++;
                        printCol++;
                    } else {
                        out << " ";
                        printCol++;
                    }
                while( printCol < nCols){
                    out << " "; printCol++;
                }
                out <<  std::endl;
            }
        };

        // Print Local values
        std::cout << "Sparsity:" << std::endl;
        // std::cout << std::string(this->nLocalColumns, '=') << std::endl;
        printLocal(this->columnIndices, this->rowPtr, this->nLocalColumns);
        for (unsigned int proc = 1; proc < this->linalg->size(); proc++){
            // std::cout << std::string(this->nLocalColumns, '-') << std::endl;
            std::vector<unsigned int> otherCol(storageSize[3*proc]),
                                      otherRow(storageSize[3*proc+1]);
            MPI_Recv(otherCol.data(), otherCol.size(), MPI_UNSIGNED, proc, 1, MPI_COMM_WORLD, nullptr);
            MPI_Recv(otherRow.data(), otherRow.size(), MPI_UNSIGNED, proc, 2, MPI_COMM_WORLD, nullptr);
            printLocal(otherCol, otherRow, storageSize[3*proc+2]);
        }
        // std::cout << std::string(this->nLocalColumns, '=') << std::endl;
    } else {
        MPI_Send(this->columnIndices.data(), this->columnIndices.size(), MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD);
        MPI_Send(this->rowPtr.data(), this->rowPtr.size(), MPI_UNSIGNED, 0, 2, MPI_COMM_WORLD);
    }
}