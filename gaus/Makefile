libge: matrix.o pivot.o piv_ge_solver.o
	ar rvs libge.a $^

matrix.o: matrix.c matrix.h
pivot.o: pivot.c matrix.h
piv_ge_solver.o: piv_ge_solver.c piv_ge_solver.h matrix.h

.PHONY: clean

conj_grad_solver.o: conj_grad_solver.c conj_grad_solver.h matrix.o
	gcc -c conj_grad_solver.c conj_grad_solver.h matrix.o

test: matrix.o conj_grad_solver.o
	gcc matrix.o conj_grad_solver.o

clean:
	-rm *.o libge.a a.out *gch
