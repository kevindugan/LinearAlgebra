
#include "Matrix_CycleRowPartition.h"
#include "math.h"
#include <iomanip>

Matrix_CycleRowPartition::Matrix_CycleRowPartition(unsigned int nRows,
                                                   unsigned int nCols,
                                                   const LinearAlgebra &linalg){

    // Calculate decomposed row cycle size
    float ratio = float(nRows) / float(linalg.size());
    unsigned int splitIndex = nRows - linalg.size() * floor( ratio );
    this->nLocalRows = (linalg.rank() < splitIndex) ? ceil( ratio ) : floor( ratio );

    this->nLocalColumns = nCols;
    this->nGlobalRows = nRows;
    this->nGlobalColumns = nCols;

    // Calculating global index range for local rows
    this->globalRowIndexRange.begin = linalg.rank();
    this->globalRowIndexRange.end = this->nGlobalRows;
    this->globalRowIndexRange.skip = linalg.size();

    this->linalg = &linalg;
    this->matrixStorage = new double*[this->nLocalRows];
    for (unsigned int row = 0; row < this->nLocalRows; row++){
        this->matrixStorage[row] = new double[this->nLocalColumns];
        // Initialize entries to zero
        for (unsigned int col = 0; col < this->nLocalColumns; col++)
            this->matrixStorage[row][col] = 0.0;
    }
}

Matrix_CycleRowPartition::~Matrix_CycleRowPartition(){
    for (unsigned int i = 0; i < this->nLocalRows; i++)
        delete[] this->matrixStorage[i];
    delete[] this->matrixStorage;
}

std::vector<unsigned int> Matrix_CycleRowPartition::getPartitionSize() const {
    std::vector<unsigned int> result(this->linalg->size());
    MPI_Allgather(&this->nLocalRows, 1, MPI_UNSIGNED, result.data(), 1, MPI_UNSIGNED, MPI_COMM_WORLD);
    return result;
}

void Matrix_CycleRowPartition::setValues(const double &x) {
    for (unsigned int i = 0; i < this->nLocalRows; i++)
        for (unsigned int j = 0; j < this->nLocalColumns; j++)
            this->matrixStorage[i][j] = x;
}

void Matrix_CycleRowPartition::setValues(const std::vector<std::vector<double>> &x) {
    Nucleus_ASSERT_EQ(x.size(), this->nGlobalRows)
    
    for (unsigned int i = 0, j = this->globalRowIndexRange.begin; i < this->nLocalRows; i++, j+=this->globalRowIndexRange.skip){
        Nucleus_ASSERT_EQ(x[j].size(), this->nGlobalColumns)
        for (unsigned int k = 0; k < this->nLocalColumns; k++)
            this->matrixStorage[i][k] = x[j][k];
    }
}

void Matrix_CycleRowPartition::zeros() {
    for (unsigned int i = 0; i < this->nLocalRows; i++)
        for (unsigned int j = 0; j < this->nLocalColumns; j++)
            this->matrixStorage[i][j] = 0.0;
}

double Matrix_CycleRowPartition::frobeniusNorm() const {
    double local_result = 0.0;
    for (unsigned int i = 0; i < this->nLocalRows; i++)
        for (unsigned int j = 0; j < this->nLocalColumns; j++)
            local_result += this->matrixStorage[i][j] * this->matrixStorage[i][j];

    double result = 0.0;
    MPI_Allreduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return sqrt(result);
}

std::unique_ptr<AbstractVector> Matrix_CycleRowPartition::mult(const AbstractVector &other) const {
    Nucleus_ASSERT_EQ(other.size(), this->nGlobalColumns);
    std::unique_ptr<AbstractVector> result = std::make_unique<Vector_CyclePartition>(this->nGlobalRows, *this->linalg);

    std::vector<double> global_values(other.size());
    other.getGlobalValues(global_values);

    std::vector<double> result_values(result->size());
    for (unsigned int row = 0, global_row = this->globalRowIndexRange.begin; row < this->nLocalRows; row++, global_row+=this->globalRowIndexRange.skip)
        for (unsigned int col = 0; col < this->nLocalColumns; col++)
            result_values[global_row] += this->matrixStorage[row][col] * global_values[col];

    result->setValues(result_values);

    return result;
}

void Matrix_CycleRowPartition::print(std::ostream &out) const {

    
    std::vector<unsigned int> part = this->getPartitionSize();
    if (this->linalg->rank() == 0){

        std::vector<std::string> printRows(this->nGlobalRows);

        // Fill in rank 0 values
        for (unsigned int i = 0, j = 0; i < this->nLocalRows; i++, j+= this->linalg->size()){
            std::stringstream stream;
            for (unsigned int k = 0; k < this->nLocalColumns; k++){
                stream << std::setw(6) << std::setprecision(4) << this->matrixStorage[i][k] << ", ";
            }
            printRows[j] = stream.str() + "\n";
        }

        // Fill in remaining rank values
        for (unsigned int proc = 1; proc < this->linalg->size(); proc++){
            // Communicate
            double *dump = new double[part[proc] * this->nLocalColumns];
            MPI_Recv(dump, part[proc] * this->nLocalColumns, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for (unsigned int i = 0, j = proc; i < part[proc]; i++, j += this->linalg->size()){
                std::stringstream stream;
                for (unsigned int k = 0; k < this->nLocalColumns; k++){
                    stream << std::setw(6) << std::setprecision(4) << dump[i*this->nLocalColumns + k] << ", ";
                }
                printRows[j] = stream.str() + "\n";
            }
            delete[] dump;
        }

        for (const auto &line : printRows)
            out << line.c_str();

    } else {
        double* send = new double[this->nLocalRows * this->nLocalColumns];
        unsigned int index = 0;
        for (unsigned int row = 0; row < this->nLocalRows; row++)
            for (unsigned int col = 0; col < this->nLocalColumns; col++, index++)
                send[index] = this->matrixStorage[row][col];

        MPI_Send(send, this->nLocalRows * this->nLocalColumns, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        delete[] send;
    }
}