#include "Vector_BlockPartition.h"
#include <string>
#include "math.h"

Vector_BlockPartition::Vector_BlockPartition(unsigned int size, const LinearAlgebra& linalg){
    // Calculate the decomposed vector block size
    float ratio = float(size) / float(linalg.size());
    unsigned int splitIndex = size - linalg.size() * floor( ratio );
    this->local_size = (linalg.rank() < splitIndex) ? ceil( ratio ) : floor( ratio );

    // Calculate the global index range for local vector
    if (linalg.rank() < splitIndex){
        this->globalIndexRange.begin = linalg.rank() * this->local_size;
        this->globalIndexRange.end = (linalg.rank() + 1) * this->local_size;
    } else {
        unsigned int offset = splitIndex * (this->local_size + 1);
        this->globalIndexRange.begin = (linalg.rank() - splitIndex) * this->local_size + offset;
        this->globalIndexRange.end = (linalg.rank() - splitIndex +1) * this->local_size + offset;
    }

    this->global_size = size;
    this->values = new double[this->local_size];

    this->linalg = &linalg;

    for (unsigned int i = 0; i < this->local_size; i++)
        this->values[i] = 0.0;
}

Vector_BlockPartition::Vector_BlockPartition(const Vector_BlockPartition &other) {
    this->local_size = other.local_size;
    this->global_size = other.global_size;
    this->globalIndexRange = other.globalIndexRange;
    this->linalg = other.linalg;

    this->values = new double[this->local_size];
    for (unsigned int i = 0; i < this->local_size; i++)
        this->values[i] = other.values[i];
}

Vector_BlockPartition& Vector_BlockPartition::operator=(const Vector_BlockPartition &other) {
    if (this != &other){
        this->local_size = other.local_size;
        this->global_size = other.global_size;
        this->globalIndexRange = other.globalIndexRange;
        this->linalg = other.linalg;

        delete[] this->values;
        this->values = new double[this->local_size];
        for (unsigned int i = 0; i < this->local_size; i++)
            this->values[i] = other.values[i];
    }
    return *this;
}

Vector_BlockPartition::~Vector_BlockPartition(){
    delete[] this->values;
}

void Vector_BlockPartition::print() const {
    if (this->linalg->rank() == 0){
        std::cout << "=================\n";
        for (unsigned int i = 0; i < this->local_size; i++)
            std::cout << "  " << this->values[i] << "\n";
        std::cout << "-----------------\n";
        for (unsigned int proc = 1; proc < this->linalg->size(); proc++){
            double* buff = new double[this->local_size];
            MPI_Recv(buff, this->local_size, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (unsigned int index = 0; index < this->local_size; index++)
                std::cout << "  " << buff[index] << "\n";
            if ( proc != this->linalg->size()-1)
                std::cout << "-----------------\n";
            delete[] buff;
        }
        std::cout << "=================" << std::endl;

    } else {
        MPI_Send(this->values, this->local_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
}

std::vector<int> Vector_BlockPartition::getPartitionSize() const {
    std::vector<int> result(this->linalg->size());
    MPI_Allgather(&this->local_size, 1, MPI_UNSIGNED, result.data(), 1, MPI_UNSIGNED, MPI_COMM_WORLD);
    return result;
}

void Vector_BlockPartition::zeros() {
    for (unsigned int i = 0; i < this->local_size; i++)
        this->values[i] = 0.0;
}

void Vector_BlockPartition::setValues(const double &x){
    for (unsigned int i = 0; i < this->local_size; i++)
        this->values[i] = x;
}

void Vector_BlockPartition::setValues(const std::vector<double> &x){
    Nucleus_ASSERT_EQ(x.size(), this->global_size)
    for (unsigned int i = 0, j = this->globalIndexRange.begin; i < this->local_size; i++, j++)
        this->values[i] = x[j];
}

void Vector_BlockPartition::scale(const double &x){
    for (unsigned int i = 0; i < this->local_size; i++)
        this->values[i] *= x;
}

unsigned int Vector_BlockPartition::findRankWithIndex(const unsigned int index) const {
    float ratio = float(this->global_size) / float(this->linalg->size());
    unsigned int splitRank = this->global_size - this->linalg->size() * floor( ratio );
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

double Vector_BlockPartition::getValue(const unsigned int index) const {
    linalg_assert(index < this->global_size)

    // Calculate Rank containing index
    unsigned int indexOnRank = this->findRankWithIndex(index);

    double result;
    if (indexOnRank == this->linalg->rank())
        result = this->values[index - this->globalIndexRange.begin];

    MPI_Bcast(&result, 1, MPI_DOUBLE, indexOnRank, MPI_COMM_WORLD);

    return result;
}

void Vector_BlockPartition::getGlobalValues(std::vector<double> &global_values) const {
    Nucleus_ASSERT_EQ(this->global_size, global_values.size())
    std::vector<int> block_sizes = this->getPartitionSize();
    std::vector<int> offset(this->linalg->size());
    for (int i = 0, sum = 0; i < offset.size(); i++, sum+=block_sizes[i-1])
        offset[i] = sum;

    MPI_Allgatherv(this->values, this->local_size, MPI_DOUBLE, global_values.data(), block_sizes.data(), offset.data(), MPI_DOUBLE, MPI_COMM_WORLD);
}

Vector_BlockPartition Vector_BlockPartition::add(const Vector_BlockPartition &other) const {
    return this->add(1.0, other);
}

Vector_BlockPartition Vector_BlockPartition::add(const double &scale, const Vector_BlockPartition &other) const {
    Nucleus_ASSERT_EQ(this->global_size, other.size())
    Vector_BlockPartition result(this->global_size, *this->linalg);

    std::vector<double> vals(this->global_size);
    for (unsigned int i = 0, j = this->globalIndexRange.begin; i < this->local_size; i++, j++)
        vals[j] = this->values[i] + scale * other.values[i];

    result.setValues(vals);
    
    return result;
}

double Vector_BlockPartition::length() const {
    double local_result = 0.0;
    for (unsigned int i = 0; i < this->local_size; i++)
        local_result += this->values[i] * this->values[i];

    double global_result = 0.0;
    MPI_Allreduce(&local_result, &global_result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return sqrt(global_result);
}
