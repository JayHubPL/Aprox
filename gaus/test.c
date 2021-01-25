#include "matrix.h"
#include "conj_grad_solver.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

int main(int argc, char** argv)
{
    FILE* in = argc > 1 ? fopen(argv[1], "r") : stdin;
    FILE* b = argc > 1 ? fopen(argv[2], "r") : stdin;

    int size;
    fscanf(in, "%d", &size);
    //printf("%d", size);
    matrix_t* a = make_matrix(size, size + 1);

    double buf;

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            fscanf(in, " %lf", &buf);
            put_entry_matrix(a, i, j, buf);
        }
        fscanf(b, " %lf", &buf);
        put_entry_matrix(a, i, size, buf);

    }

    //write_matrix(a, stdout);

    int tmp = conj_grad_solver(a);
    printf("\n Wynik:%d \n", tmp);

    write_matrix(a, stdout);

}