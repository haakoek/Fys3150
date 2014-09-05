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

double source_term(double);
double exact_solution(double);
void solver(double*,double*,double*,double*,double*,int);
void rel_error(double*,double*,double*,int);
void output(double*, double*,double*, int);

ofstream ofile1;
ofstream ofile2;

int main(int argc, char* argv[])
{
    //Variable declarations

    int i,n;
    double *a, *b, *c, *g, *computed_solution,*exact,*relative_error;
    ofile1.open("output1.txt", ios::out);
    ofile2.open("output2.txt");

    //Read problem size from command line
    n = atoi(argv[1]);
    cout << n << endl;

    //Initialize variables

    double h = 1.0/(n+1);

    a = new double[n+2];
    b = new double[n+2];
    c = new double[n+2];
    g =  new double[n+2];
    relative_error = new double[n+2];
    computed_solution = new double[n+2];
    exact = new double[n+2];

    for(i=0; i < n+2; i++) {
        g[i] = pow(h,2)*source_term(i*h);
        a[i] = 2.0;
        b[i] = -1.0;
        c[i] = -1.0;
        relative_error[i] = 0.0;
        exact[i] = exact_solution(i*h);
    }

    //Set boundary conditions for the computed solution

    computed_solution[0] = 0.0;
    computed_solution[n+1] = 0.0;

    //Take time of the solver function, compute the solution, and extract the relative error,
    //and write to file.

    clock_t begin1 = clock();

    solver(a,b,c,g,computed_solution,n);

    clock_t end1 = clock();
    double elapsed_secs1 = double(end1 - begin1) / CLOCKS_PER_SEC;

    rel_error(relative_error, exact, computed_solution,n);

    ofile1 << setiosflags(ios::showpoint | ios::uppercase);
    ofile1 << setw(15) << setprecision(8) << relative_error[n/2];
    ofile1 << setw(15) << setprecision(8) << elapsed_secs1;

    //If the problem size is less then 10000, we also want to compute the solution,
    //using armadillos solve, and lu decomposition.

    if(n <= 10000) {

        mat L = zeros<mat>(n,n);
        mat U = zeros<mat>(n,n);
        mat P = zeros<mat>(n,n);
        mat A = zeros<mat>(n,n);
        mat B = zeros<mat>(n,1);
        //mat C = randu<mat>(n,n);

        for(i=0; i < n; i++) {
            A(i,i) = 2.0;
            B(i,0) = g[i+1];
            if(i < n-1) {
                A(i,i+1) = -1.0;
                A(i+1,i) = -1.0;
            }
        }
        /*
        for(i=0; i < n; i++) {
            B(i,0) = g[i+1];
        }
        */

        clock_t begin2 = clock();

        //Standard Gauss-Jordan elimination

        mat X = solve(A,B);


        clock_t end2 = clock();

        double elapsed_secs2 = double(end2 - begin2) / CLOCKS_PER_SEC;

        //LU decomposition

        lu(L,U,P,A);

        mat y = solve(L,P*B);

        clock_t begin3 = clock();

        mat x = solve(U,y);

        clock_t end3 = clock();

        double elapsed_secs3 = double(end3 - begin3) / CLOCKS_PER_SEC;

        //Output

        ofile1 << setw(15) << setprecision(8) << elapsed_secs2;
        ofile1 << setw(15) << setprecision(8) << elapsed_secs3;
    }

    //output(computed_solution,exact,relative_error, n);

    ofile1.close();
    ofile2.close();

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

