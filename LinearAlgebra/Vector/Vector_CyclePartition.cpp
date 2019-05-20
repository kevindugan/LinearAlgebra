#include "Vector_CyclePartition.h"
#include "math.h"
#include <iomanip>

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
    this->local_size = other.local_size;
    this->global_size = other.global_size;
    this->globalIndexRange = other.globalIndexRange;
    this->linalg = other.linalg;

    this->values = new double[this->local_size];
    for (unsigned int i = 0; i < this->local_size; i++)
        this->values[i] = other.values[i];
}

Vector_CyclePartition& Vector_CyclePartition::operator=(const Vector_CyclePartition &other) {
    if (this != &other){
        this->local_size = other.local_size;
        this->global_size = other.global_size;
        this->globalIndexRange = other.globalIndexRange;
        this->linalg = other.linalg;

        this->values = new double[this->local_size];
        for (unsigned int i = 0; i < this->local_size; i++)
            this->values[i] = other.values[i];
    }
    return *this;
}

Vector_CyclePartition::~Vector_CyclePartition() {
    delete[] this->values;
}

std::unique_ptr<AbstractVector> Vector_CyclePartition::clone() const {
    return std::make_unique<Vector_CyclePartition>(*this);
}

void Vector_CyclePartition::setValues(const double &x) {
    for (unsigned int i = 0; i < this->local_size; i++)
        this->values[i] = x;
}

void Vector_CyclePartition::setValues(const std::vector<double> &x) {
    Nucleus_ASSERT_EQ(x.size(), this->global_size)
    for (unsigned int i = 0, j = this->globalIndexRange.begin; i < this->local_size; i++, j+=globalIndexRange.skip)
        this->values[i] = x[j];
}

void Vector_CyclePartition::scale(const double &x) {
    for (unsigned int i = 0; i < this->local_size; i++)
        this->values[i] *= x;
}

void Vector_CyclePartition::zeros() {
    for (unsigned int i = 0; i < this->local_size; i++)
        this->values[i] = 0.0;
}

std::unique_ptr<AbstractVector> Vector_CyclePartition::add(const AbstractVector &other) const {
    return this->add(1.0, other);
}

std::unique_ptr<AbstractVector> Vector_CyclePartition::add(const double &scale, const AbstractVector &other) const {
    Nucleus_ASSERT_EQ(this->global_size, other.size())
    std::unique_ptr<AbstractVector> result = std::make_unique<Vector_CyclePartition>(this->global_size, *this->linalg);

    std::vector<double> vals(this->global_size);
    for (unsigned int i = 0, j = this->globalIndexRange.begin; i < this->local_size; i++, j+=this->globalIndexRange.skip)
        vals[j] = this->values[i] + scale * other.getLocalValue(i);

    result->setValues(vals);

    return result;
}

double Vector_CyclePartition::l2norm() const {
    double local_result = 0.0;
    for (unsigned int i = 0; i < this->local_size; i++)
        local_result += this->values[i] * this->values[i];

    double global_result = 0.0;
    MPI_Allreduce(&local_result, &global_result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return sqrt(global_result);
}

double Vector_CyclePartition::getValue(const unsigned int index) const {
    Nucleus_ASSERT_LT(index, this->global_size)

    int indexOnRank = this->findRankWithIndex(index);
    double local_value;
    if (indexOnRank == this->linalg->rank())
        local_value = this->values[index / this->linalg->size()];

    MPI_Bcast(&local_value, 1, MPI_DOUBLE, indexOnRank, MPI_COMM_WORLD);

    return local_value;
}

void Vector_CyclePartition::getGlobalValues(std::vector<double> &global_values) const {
    Nucleus_ASSERT_EQ(this->global_size, global_values.size())
    std::vector<int> block_sizes = this->getPartitionSize();
    std::vector<int> offset(this->linalg->size());
    for (int i = 0, sum = 0; i < offset.size(); i++, sum+=block_sizes[i-1])
        offset[i] = sum;

    std::vector<double> block_data(this->global_size);

    MPI_Allgatherv(this->values, this->local_size, MPI_DOUBLE, block_data.data(), block_sizes.data(), offset.data(), MPI_DOUBLE, MPI_COMM_WORLD);

    // MPI_Allgather takes each ranks data and squishes it together in block_data
    // Now we need to distribute the values correctly in global_values.
    // [[R0.1 R0.2 R0.3 .. R0.N] [R1.1 R1.2 R1.3 .. R1.N] ... [Rp.1 Rp.2 Rp.3 .. Rp.N]] =>
    // [[R0.1 R1.1 R2.1 .. Rp.1] [R0.2 R1.2 R2.2 .. Rp.2] ... [R0.N R1.N R2.N .. Rp.N]]
    for (unsigned int i = 0; i < block_sizes.size(); i++)
        for (unsigned int j = 0; j < block_sizes[i]; j++)
            global_values[i + j * this->linalg->size()] = block_data[offset[i] + j];

}

std::vector<int> Vector_CyclePartition::getPartitionSize() const {
    std::vector<int> result(this->linalg->size());
    MPI_Allgather(&this->local_size, 1, MPI_UNSIGNED, result.data(), 1, MPI_UNSIGNED, MPI_COMM_WORLD);
    return result;
}

unsigned int Vector_CyclePartition::findRankWithIndex(const unsigned int index) const {
    Nucleus_ASSERT_LT(index, this->global_size)
    return index % this->linalg->size();
}

double Vector_CyclePartition::getLocalValue(const unsigned int local_index) const {
    return this->values[local_index];
}

void Vector_CyclePartition::print(std::ostream &out) const {
    std::vector<int> part = this->getPartitionSize();
    if (this->linalg->rank() == 0){
        std::vector<std::string> printValues(this->global_size);

        // Fill in Rank 0 values
        for (unsigned int local = 0, global = this->globalIndexRange.begin; local < this->local_size; local++, global+=this->globalIndexRange.skip){
            std::stringstream stream;
            stream << std::setw(9) << std::setprecision(2) << std::scientific << this->values[local] << "\n";
            printValues[global] = stream.str();
        }

        // Fill in other Rank values
        for (unsigned int proc = 1; proc < this->linalg->size(); proc++){
            double* result = new double[part[proc]];
            MPI_Recv(result, part[proc], MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (unsigned int local = 0, global = proc; local < part[proc]; local++, global+=this->linalg->size()){
                std::stringstream stream;
                stream << std::setw(9) << std::setprecision(2) << std::scientific << result[local] << "\n";
                printValues[global] = stream.str();
            }
            delete[] result;
        }

        for (const auto &line : printValues)
            out << line.c_str();

    } else {
        MPI_Send(this->values, this->local_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
}
