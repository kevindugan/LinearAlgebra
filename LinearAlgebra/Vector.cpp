#include "Vector.h"
#include <string>
#include "math.h"

Vector::Vector(unsigned int size, const LinearAlgebra& linalg){
    linalg_assert( (size % linalg.size()) == 0)

    this->global_size = size;
    this->local_size = size / linalg.size();
    this->global_starting_index = linalg.rank() * this->local_size;
    this->values = new double[this->local_size];

    this->linalg = &linalg;

    for (unsigned int i = 0; i < this->local_size; i++)
        this->values[i] = double(i + this->global_starting_index);
}

Vector::~Vector(){
    delete[] this->values;
}

void Vector::print() const {
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

void Vector::setValues(const double &x){
    for (unsigned int i = 0; i < this->local_size; i++)
        this->values[i] = x;
}

void Vector::scale(const double &x){
    for (unsigned int i = 0; i < this->local_size; i++)
        this->values[i] *= x;
}

double Vector::getLocal(const unsigned int index) const {
    return this->values[index];
}

void Vector::add(const Vector &x){
    linalg_assert(this->global_size == x.size())
    
    for (unsigned int i = 0; i < this->local_size; i++)
        this->values[i] += x.getLocal(i);
}

double Vector::length() const {
    double result = 0.0;
    for (unsigned int i = 0; i < this->local_size; i++)
        result += this->values[i] * this->values[i];

    if (this->linalg->rank() == 0){
        for (unsigned int proc = 1; proc < this->linalg->size(); proc++){
            double dump;
            MPI_Recv(&dump, 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            result += dump;
        }
    } else {
        MPI_Send(&result, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    return sqrt(result);
}