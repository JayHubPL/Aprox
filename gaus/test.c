#include "matrix.h"
#include "piv_ge_solver.h"
#include "conj_grad_solver.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

int main()
{
    matrix_t* a = make_matrix(3, 4);
    //matrix_t* b = make_matrix(3, 1);


    put_entry_matrix(a, 0, 0, 1);
    put_entry_matrix(a, 0, 1, 4);
    put_entry_matrix(a, 0, 2, 3);
    put_entry_matrix(a, 0, 3, 1);

    put_entry_matrix(a, 1, 0, 4);
    put_entry_matrix(a, 1, 1, 2);
    put_entry_matrix(a, 1, 2, 5);
    put_entry_matrix(a, 1, 3, 2);

    put_entry_matrix(a, 2, 0, 3);
    put_entry_matrix(a, 2, 1, 5);
    put_entry_matrix(a, 2, 2, 3);
    put_entry_matrix(a, 2, 3, 3);



    write_matrix(a, stdout);

    int tmp = conj_grad_solver(a);



    printf("\n Wynik:%d \n", tmp);

    write_matrix(a, stdout);
    // write_matrix(b, stdout);

}