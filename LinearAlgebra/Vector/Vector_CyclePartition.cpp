#include "Vector_CyclePartition.h"
#include "math.h"

Vector_CyclePartition::Vector_CyclePartition(unsigned int size, const LinearAlgebra &linalg) {
    // Calculate the decomposed vector cycle size
    float ratio = float(size) / float(linalg.size());
    unsigned int splitIndex = size - linalg.size() * floor( ratio );
    this->local_size = (linalg.rank() < splitIndex) ? ceil( ratio ) : floor( ratio );
    this->global_size = size;

    // Calculate the global index range for local vector
    this->globalIndexRange.begin = linalg.rank();
    this->globalIndexRange.end = this->global_size;
    this->globalIndexRange.skip = linalg.size();

    this->values = new double[this->local_size];
    this->linalg = &linalg;

    for (unsigned int i = 0; i < this->local_size; i++)
        this->values[i] = 0.0;
}

Vector_CyclePartition::Vector_CyclePartition(const Vector_CyclePartition &other) {

}

Vector_CyclePartition& Vector_CyclePartition::operator=(const Vector_CyclePartition &other) {

}

Vector_CyclePartition::~Vector_CyclePartition() {
    delete[] this->values;
}

void Vector_CyclePartition::setValues(const double &x) {

}

void Vector_CyclePartition::setValues(const std::vector<double> &x) {

}

void Vector_CyclePartition::scale(const double &x) {

}

void Vector_CyclePartition::zeros() {

}

std::unique_ptr<AbstractVector> Vector_CyclePartition::add(const AbstractVector &other) const {
    return std::make_unique<Vector_CyclePartition>(this->global_size, *this->linalg);
}

std::unique_ptr<AbstractVector> Vector_CyclePartition::add(const double &scale, const AbstractVector &other) const {
    return std::make_unique<Vector_CyclePartition>(this->global_size, *this->linalg);
}

double Vector_CyclePartition::l2norm() const {
    return 0.0;
}

double Vector_CyclePartition::getValue(const unsigned int index) const {
    return 0.0;
}

void Vector_CyclePartition::getGlobalValues(std::vector<double> &global_values) const {

}

std::vector<int> Vector_CyclePartition::getPartitionSize() const {
    std::vector<int> result(this->linalg->size());
    MPI_Allgather(&this->local_size, 1, MPI_UNSIGNED, result.data(), 1, MPI_UNSIGNED, MPI_COMM_WORLD);
    return result;
}

unsigned int Vector_CyclePartition::findRankWithIndex(const unsigned int index) const {
    return 0;
}

double Vector_CyclePartition::getLocalValue(const unsigned int local_index) const {
    return 0.0;
}