/*
 * Project1 Fys3150
 * Authors: Haakon Kristiansen and Stian Goplen
 * Program to solve the 1-dimensional Poisson equation
 * with Dirichlet boundary conditions.
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <armadillo>

using namespace std;
using namespace arma;

//Function declarations

double source_term(double);
double exact_solution(double);
void solver(double*,double*,double*,double*,double*,int);
void rel_error(double*,double*,double*,int);
void output(double*, double*,double*, int);

//Two output files may be needed, and are declared as global variables

ofstream ofile1;
ofstream ofile2;

int main(int argc, char* argv[])
{
    //Variable declarations

    int i,n;
    double *a, *b, *c, *g, *computed_solution,*exact,*relative_error;
    char *outfile1;
    char *outfile2;

    /* Read problem size and at least one outputfile from command line,
     * and abort if there are too few command-line arguments.
    */

    if(argc <= 2) {
        cout << "Bad usage: " << argv[0] << endl;
        cout << "Read problem size, " << endl;
        cout << "and at least one output file on command line" << endl;
        exit(1);
    } else {
        n = atoi(argv[1]);
        outfile1 = argv[2];
        ofile1.open(outfile1, ios::out);
        cout << n << endl;
    }

    //Initialize variables

    double h = 1.0/(n+1);

    a = new double[n+2];
    b = new double[n+2];
    c = new double[n+2];
    g =  new double[n+2];
    relative_error = new double[n+2];
    computed_solution = new double[n+2];
    exact = new double[n+2];

    /*
     * Compute right handside, initialize a,b,c, and relative_error
     * and finally compute the exact solution.
    */

    for(i=0; i < n+2; i++) {
        g[i] = pow(h,2)*source_term(i*h);
        a[i] = 2.0;
        b[i] = -1.0;
        c[i] = -1.0;
        relative_error[i] = 0.0;
        exact[i] = exact_solution(i*h);
    }

    //Set boundary conditions for the computed solution (Dirichlet).

    computed_solution[0] = 0.0;
    computed_solution[n+1] = 0.0;

    //Take time of the solver function, compute the solution, and extract the relative error,
    //and write the relative error at n/2, and the execution time to outfile1.

    clock_t begin1 = clock();

    solver(a,b,c,g,computed_solution,n);

    clock_t end1 = clock();
    double elapsed_secs1 = double(end1 - begin1) / CLOCKS_PER_SEC;

    rel_error(relative_error, exact, computed_solution,n);

    //cout << relative_error[n/2] << endl;

    ofile1 << setiosflags(ios::showpoint | ios::uppercase);
    ofile1 << setw(15) << setprecision(8) << relative_error[n/2];
    ofile1 << setw(15) << setprecision(8) << elapsed_secs1;

    if(argc >= 3) {

        /* If we want to write the computed solution, and exact solution to file
         * provide a second output file on the command line.
        */

        outfile2 = argv[3];
        ofile2.open(outfile2, ios::out);
        output(computed_solution,exact,relative_error, n);
    }

    //Clear arrays which are not longer needed, g is needed below

    delete [] a;
    delete [] b;
    delete [] c;
    delete [] relative_error;
    delete [] computed_solution;
    delete [] exact;

    //If the problem size is less then 10000, we also want to compute the solution,
    //using armadillos solve, and lu decomposition.

    if(n <= 10000) {

        /*Declare matrices, using armadillo.
         * L, U and P is used in the LU decomposition
         * A is the second-derivative matrix.
         * B is the right hand side in Ax=B
         * C is a randomally generated matrix, used to test the run time of armadillos solve
         * where the matrix is not necesserally tridiagonal
        */

        mat L = zeros<mat>(n,n);
        mat U = zeros<mat>(n,n);
        mat P = zeros<mat>(n,n);
        mat A = zeros<mat>(n,n);
        mat B = zeros<mat>(n,1);
        mat C = randu<mat>(n,n);

        for(i=0; i < n; i++) {
            A(i,i) = 2.0;
            B(i,0) = g[i+1];
            if(i < n-1) {
                A(i,i+1) = -1.0;
                A(i+1,i) = -1.0;
            }
        }

        /*Standard Gauss-Jordan elimination
         * X is the solution to AX=B.
        */

        clock_t begin2 = clock();

        mat X = solve(A,B);

        clock_t end2 = clock();

        double elapsed_secs2 = double(end2 - begin2) / CLOCKS_PER_SEC;


        /* LU decomposition
         * First decompose A into, st PA=LU, calling lu(L,U,P,A)
         * Then solve Ly=PB
         * Finally to obtain a solution solve Ux=y
         * We are only interested in the time it takes to solve the last system, so we only clock
         * the time it takes to solve Ux=y.
        */

        lu(L,U,P,A);

        mat y = solve(L,P*B);

        clock_t begin3 = clock();

        mat x = solve(U,y);

        clock_t end3 = clock();

        double elapsed_secs3 = double(end3 - begin3) / CLOCKS_PER_SEC;

        /*This second Gauss solution is of no significance to the actual problem,
          * but to test the exectuion time, as it seems that armadillo recognizes
          * the tridiagonal property of A, and reduces the exection time for GJ, from
          * O(n**3) to O(n**2). Hence we want to check that for a random generated matrix C, the
          * solution is actually O(n**3) for GJ as we would expect.
          * Y is the solution of CY=B.
        */


        clock_t begin4 = clock();

        mat Y = solve(C,B);

        clock_t end4 = clock();

        double elapsed_secs4 = double(end4 - begin4) / CLOCKS_PER_SEC;

        //Finally we write the exection times to outfile1.

        ofile1 << setw(15) << setprecision(8) << elapsed_secs2;
        ofile1 << setw(15) << setprecision(8) << elapsed_secs3;
        ofile1 << setw(15) << setprecision(8) << elapsed_secs4;



    }

    //Close files, and free memory.

    ofile1.close();
    ofile2.close();

    delete [] g;


    return 0;
}

double source_term(double x) {
    return 100*exp(-10*x);
}

double exact_solution(double x) {
    return 1.0-(1.0-exp(-10.0))*x - exp(-10*x);
}

void solver(double *a, double *b, double *c, double *g,double *sol,int n) {

    int k;

    for(k=2; k < n+1; k++) {
        a[k] = a[k] - (b[k-1]/(a[k-1]))*c[k-1];
        g[k] = g[k] - (b[k-1]/a[k-1])*g[k-1];
    }

    sol[n] = g[n]/a[n];

    for(k=n-1; k>=1; k--) {
        sol[k] = (g[k] - c[k]*sol[k+1])/a[k];
    }

}

void rel_error(double* relative_error,double* exact, double* comp_sol, int n) {
    int i;
    for(i=1; i < n+1;i++) {
        relative_error[i] = log10(abs((exact[i]-comp_sol[i])/(exact[i])));
    }
}

void output(double *v,double *u, double*rel, int num) {

    int i;
    ofile2 << setiosflags(ios::showpoint | ios::uppercase);

    for(i=0; i < num+2; i++) {
        ofile2 << setw(15) << setprecision(8) << v[i];
        ofile2 << setw(15) << setprecision(8) << u[i];
        ofile2 << setw(15) << setprecision(8) << rel[i] << endl;
    }
}

