#include "GMRes.h"

GMRes::GMRes(const Parallel &init){
    this->linalg = &init;
}

std::unique_ptr<AbstractVector> GMRes::solve(const AbstractMatrix &A, const AbstractVector &b) const {
    Nucleus_ASSERT_EQ(A.nRows(), b.size())
    std::unique_ptr<AbstractVector> result = A.mult(b);

    // If rhs is sufficiently small, return zero vector as solution
    double beta = b.l2norm();
    if (beta < 1e-8){
        result->zeros();
        return result;
    }

    // Initialize Subspace and Hessenberg Matrix
    std::unique_ptr<AbstractVector> v = b.clone();
    v->scale( 1.0/beta );
    // TODO BasisMatrix Space(v)
    // TODO HessenbergMatrix Hess(beta)
    double previousResidualNorm = beta;

    // Solution Loop
    bool done = false;
    unsigned int iterations = 1;
    do {

        // Apply RHS Preconditioner
        std::unique_ptr<AbstractVector> z = v->clone();

        // Perform Matrix Vector Multiply
        std::unique_ptr<AbstractVector> w = A.mult(*z);

        // Get the projection constants for the Grahm-Schmidt process. This is a vector of dot products
        // of the new Krylov vector with each vector from the existing space.
        // std::unique_ptr<AbstractVector> h = Space.Tmult(*w);

        // Remove the components of the previous space vectors from the new Krylov vector. Subtract
        // the dot product of the new vector with each space vector times each space vector from the
        // Krylov vector. Performs the following operation [v = w - Space * h]. This is the first
        // orthogonalization step.
        // std::unique_ptr<AbstractVector> v_prime = Space.mult(*h);
        // v_prime->scale(-1.0);
        // v_prime = v_prime->add(*w);

        // Perform the reorthogonalization by repeating the previous two steps, starting with on the vector obtained from
        // the previous steps.
        // std::unique_ptr<AbstractVector> h_prime = Space.Tmult(*v_prime);
        // v = Space.mult(*h_prime);
        // v->scale(-1.0);
        // v = v->add(*v_prime);
        // h = h->add(*h_prime);
        // v_prime.reset();

        // Normalize the new vector
        double len = v->l2norm();
        if (len > 1e-15)
            v->scale( 1.0/len );

        // Append the Hessenberg Matrix H|h,len
        // Vector<Real> H(h.size()+1);
        // for (unsigned int i = 0; i < h.size(); i++)
        //     H.set(i,h.get(i));
        // H.set(h.size(), len);
        // Hessenberg.updateMatrixAndVector(H);

        // Check that residual norm is decreasing monotonically.
        // Vector<Real> y = Hessenberg.solveLeastSquares();
        // Vector<Real> tmp(y.size()+1);
        // Hessenberg.mult(tmp, y);
        // tmp.set(0, tmp.get(0) - beta);
        // double residualNorm = tmp.l2norm();
        // double residualChange = fabs(residualNorm - previousResidualNorm);

        // Check for convergence
        // if (fabs(residualNorm) < tolerance ){
        //     done = true;
        //     converged = true;
        // }
        // Special case when the residual norm is exactly zero can mess up the GMRes solve routine.
        // if (residualNorm == 0.0){
        //     done = true;
        //     converged = true;
        // }
        // if (iterations > maxIterations){
        //     done = true;
        //     converged = false;
        // }
        // if (iterations == b.size()){
        //     done = true;
        //     Niterations = true;
        // }
        // previousResidualNorm = residualNorm;

        // If not converged, Append the Subspace Matrix S|v
        if (!done){
            // Space.appendColumn(*v);
            iterations++;
        }

        // Cleanup
        z.reset();
        w.reset();

    } while(!done);

    // Solve the minimization problem to get solution
    // Vector<Real> y = Hessenberg.solveLeastSquares();

    // Multiply solution to the minimization problem by the space vectors
    // and apply the preconditioner inverse to obtain the solution.
    // Vector<Real> xTmp(Space.getNumberRows());
    // Space.mult(xTmp, y);
    // M->returnInverseMultVec(solution, xTmp);
    v.reset();

    return result;
}