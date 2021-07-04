#include "GaussianElimination.h"

GaussianElimination::GaussianElimination(const LinearAlgebra &init){
    this->linalg = &init;
}

std::unique_ptr<AbstractVector> GaussianElimination::solve(const AbstractMatrix &A, const AbstractVector &b) const {
    std::unique_ptr<AbstractVector> result = A.mult(b);
    result->zeros();

    std::unique_ptr<AbstractMatrix> localA = A.clone();
    std::unique_ptr<AbstractVector> localb = b.clone();

    IndexRange range = localA->getGlobalRowIndexRange();
    std::vector<unsigned int> part = localA->getPartitionSize();

    // Upper Triangular Transformation
    std::vector<double> rhs_values(localb->size());
    localb->getGlobalValues(rhs_values);

    for (unsigned int col = 0; col < localA->nCols(); col++){
        std::vector<double> header = localA->getRowValues(col);
        double rhs_header = localb->getValue(col);

        for (unsigned int local_row = 0, global_row = range.begin; local_row < part[this->linalg->rank()]; local_row++, global_row+=range.skip){
            if (global_row <= col)
              continue;

            double* local_row_values = &localA->getLocalRowValues(local_row);
            double factor = -1.0 * local_row_values[col] / header[col];
            for (unsigned int entry = col; entry < localA->nCols(); entry++){
                local_row_values[entry] += factor * header[entry];
            }
            // localA->setLocalRowValues(local_row, local_row_values);

            rhs_values[global_row] += factor * rhs_header;
        }
        localb->setValues(rhs_values);
    }

    // Back Substitution
    std::vector<double> result_values(result->size());
    localb->getGlobalValues(result_values);
    for (int row = localb->size()-1; row >= 0; row--){
        unsigned int rowOnRank = result->findRankWithIndex(row);

        double factor = 0.0;
        if (this->linalg->rank() == rowOnRank){
            unsigned int local_row = (row - range.begin)/range.skip;
            double denom = ( &localA->getLocalRowValues(local_row) )[row];
            factor = result_values[row] / denom;
            result_values[row] = factor;
        }
        MPI_Bcast(&factor, 1, MPI_DOUBLE, rowOnRank, MPI_COMM_WORLD);
        for (unsigned int local_row = 0, global_row = range.begin; local_row < part[this->linalg->rank()]; local_row++, global_row+=range.skip){
            if (global_row >= row)
                continue;

            double* local_row_values = &localA->getLocalRowValues(local_row);
            result_values[global_row] -= local_row_values[row] * factor;
        }
    }
    result->setValues(result_values);

    return result;
}
