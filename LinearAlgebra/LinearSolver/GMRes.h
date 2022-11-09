#ifndef GMRES_H_CC09
#define GMRES_H_CC09

#include "LinearSolver.h"
#include "SparseMatrix_BlockRowPartition.h"

class GMRes : public LinearSolver {

    public:
        GMRes(const Parallel &init);
        std::unique_ptr<AbstractVector> solve(const AbstractMatrix &A, const AbstractVector &b) const override;

    private:
        const Parallel *linalg;
};

#endif // GMRES_H_CC09