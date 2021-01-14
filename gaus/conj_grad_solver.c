#include "conj_grad_solver.h"
#include <math.h>

#define EPSILON 1e-10

int conj_grad_solver(matrix_t* eqs) {
    matrix_t* A = make_matrix(eqs->rn, eqs->cn - 1); // nightly
    for (int i = 0; i < A->rn; ++i)
        for (int j = 0; j < A->cn; ++j)
            *(A->e + i * A->cn + j) = *(eqs->e + i * eqs->cn + j);
    matrix_t* b = make_matrix(eqs->rn, 1); // nightly
    for (int i = 0; i < b->rn; ++i)
        *(b->e + i) = *(eqs->e + (i + 1) * eqs->cn - 1);
    matrix_t* x = make_matrix(eqs->rn, 1);
    matrix_t* r = copy_matrix(b); // nightly
    matrix_t* p = copy_matrix(r); // nightly
    matrix_t* newr = NULL;
    matrix_t* tmp = NULL;

    for (;;) {
        // calculate alpha_k
        matrix_t* Ap = mull_matrix(A, p);
        if (!Ap)
            return 1;
        double alpha = self_dot_matrix(r);
        if (alpha == NAN) {
            free_matrix(Ap);
            return 2;
        }
        tmp = transpose_matrix(p);
        alpha /= dot_product_matrix(tmp, Ap);
        if (alpha == NAN) {
            free_matrix(Ap);
            free_matrix(tmp);
            return 3;
        }
        free_matrix(tmp);

        // calculate x_k+1
        tmp = mull_scalar_matrix(p, alpha);
        if (!tmp) {
            free_matrix(Ap);
            return 4;
        }
        for (int i = 0; i < x->rn; ++i) // x += tmp, x = x + p*alpha
            for (int j = 0; j < x->cn; ++j)
                add_to_entry_matrix(x, i, j, get_entry_matrix(tmp, i, j));
        free_matrix(tmp);

        // calculate r_k+1 := newr
        newr = copy_matrix(r);
        tmp = mull_scalar_matrix(Ap, alpha);
        if (!tmp || !newr) {
            free_matrix(Ap);
            return 5;
        }
        for (int i = 0; i < newr->rn; ++i) // r -= tmp, r = r - A*p*alpha
            for (int j = 0; j < newr->cn; ++j)
                add_to_entry_matrix(newr, i, j, -get_entry_matrix(tmp, i, j));
        free_matrix(tmp);

        // check if solution is good enough
        if (self_dot_matrix(newr) < EPSILON) // NIGHTLY?
            break;

        // calculate beta
        double beta = self_dot_matrix(newr) / self_dot_matrix(r);
        if (beta == NAN) {
            return 6;
        }

        // calculate p_k+1
        tmp = mull_scalar_matrix(p, beta);
        for (int i = 0; i < tmp->rn; ++i) // p_k+1 = r_k+1 + beta*p_k
            for (int j = 0; j < tmp->cn; ++j)
                add_to_entry_matrix(tmp, i, j, get_entry_matrix(newr, i, j));
        free_matrix(p);
        p = tmp; // p_k becomes p_k+1

        // r_k becomes r_k+1
        free_matrix(r);
        r = newr;
    }
    for (int i = 0; i < x->rn; ++i) // przepisanie x do eqs zamiast b
        *(eqs->e + (i + 1) * eqs->cn - 1) = *(x->e + i);
    free_matrix(A);
    free_matrix(b);
    return 0;
}