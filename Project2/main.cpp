#include <iostream>
#include <armadillo>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>

using namespace std;
using namespace arma;

double potential(double);
double max_offdiag(mat,int *, int *, int);
void JacobiRotation(mat &A, int k, int l, int n,mat &R);
void output_version1(vec first_three, double rho_max,double h, int iterations);
void output_version2(int n_step, int iterations, vec first_three, double time1, double time2);
void output_version3(mat R, int n);
double coulomb_potential(double rho, double omega_r);
void normalize(mat &R, int n);
void sortEigenVectors(mat &A, mat &R, int n);

//Explanation of const.

const double fixed_h = 0.05;

ofstream ofile1;

int main(int argc, char* argv[])
{

    double rho_min = 0.0;
    double rho_max;
    double h;
    int n_step;
    double omega_r = 0.0;
    string version = argv[1];

    /*
     * Version1 keeps h fixed in order to check the results for the first three
     * eigenvalues is stable with respect to rho_max, in the one-electron case, with
     * the harmonic oscillator potential V(rho) = rho**2
     *
     * In version2 h = (rho_max - rho_min)/n_step, so we can compare the time-usage of
     * Jacobi's rotation algorithm, and armadillo's eigensolver, with respect to different
     * matrix sizes (n_step x n_step).
     *
     * In version3 we use the coulomb potential, V(rho) = (omega_r**2)*(rho**2) + 1/rho,
     * to solve the radial schrodinger equation for two electrons.
     *
     * We require the specific version to be given as the first commandline argument; version1, version2 or version3.
     *
     * In addition version1 requires rho_max as the second commandline argument.
     * Version2 requires n_step as second, and rho_max as third commandline argument.
     * Version3 requires n_step as second, rho_max as third, and omega_r as fourth commandline argument.
     * */

    if (version.compare("version1") == 0) {
        if(argc <= 2) {
            cout << "Bad usage: " << argv[0] << endl;
            cout << "Give rho_max on command line " << endl;
            exit(1);
        } else {
            rho_max = atof(argv[2]);
            h = fixed_h;
            n_step = floor(rho_max/h);
        }

    } else if (version.compare("version2") == 0) {
        if(argc <= 3) {
            cout << "Bad usage: " << argv[0] << endl;
            cout << "Give n_step as second command line argument, and rho_max as third command line argument. " << endl;
            exit(1);
        } else {
            n_step = atoi(argv[2]);
            rho_max = atof(argv[3]);
            h = (rho_max - rho_min)/n_step;
        }

    } else if (version.compare("version3") == 0) {
        if(argc <= 4) {
            cout << "Bad usage: " << argv[0] << endl;
            cout << "Give n_step as second command line argument, rho_max as third command line argument, " << endl;
            cout << "and omega_r as fourth command line argument." << endl;
            exit(1);
        } else {
            n_step = atoi(argv[2]);
            rho_max = atof(argv[3]);
            omega_r = atof(argv[4]);
            h = (rho_max - rho_min)/n_step;
        }
    } else {
        cout << "Bad usage: " << argv[0] << endl;
        cout << "Specify version of program as first command line argument. " << endl;
        cout << "version1, version2 or version3" << endl;
        exit(1);
    }


    /*Initialize problem
     *
     * If version is version1 or version2 we use the harmonic oscillator potential V(rho) = rho**2
     * If version is version3 we use the coulomb potential.
     *
     * A is ...
     * R is a matrix which is used to hold the eigenvectors in it's columns.
     * */

    mat A = zeros<mat>(n_step+1,n_step+1);
    mat R = zeros<mat>(n_step+1,n_step+1);

    if(version.compare("version1") == 0 || version.compare("version2") == 0) {
        for(int i = 1; i < n_step; i++) {
            A(i,i) = (2.0/(h*h)) + potential(rho_min + i*h);
            R(i,i) = 1.0;
            if(i < n_step-1) {
                A(i,i+1) = -1.0/(h*h);
                A(i+1,i) = -1.0/(h*h);
            }
        }
    } else if (version.compare("version3") == 0) {
        for(int i = 1; i < n_step; i++) {
            A(i,i) = (2.0/(h*h)) + coulomb_potential(rho_min + i*h, omega_r);
            R(i,i) = 1.0;
            if(i < n_step-1) {
                A(i,i+1) = -1.0/(h*h);
                A(i+1,i) = -1.0/(h*h);
            }
        }
    }
    //End initialization


    //Compute the solution to the eigenproblem A*u=(lamba)*u, using Jacobi's Rotation Algorithm.

    int col;
    int row;
    double max_off = max_offdiag(A,&col,&row,n_step);
    double epsilon = 1.0e-8;
    int iterations = 0;
    double max_nr_iterations = (double) n_step * (double) n_step * (double) n_step;

    clock_t begin1 = clock();

    while (fabs(max_off) > epsilon && (double) iterations < max_nr_iterations ) {
        JacobiRotation(A,col,row,n_step,R);
        max_off = max_offdiag(A,&col,&row,n_step);
        iterations++;
    }

    clock_t end1 = clock();

    double elapsed_secs1 = double(end1 - begin1) / CLOCKS_PER_SEC;

    //Solve with Armadillo, used to compare the time usage between Jacobi's algrithm and armadillo's eigensolver

    clock_t begin2 = clock();

    mat Arma_Solution = eig_sym(A);

    clock_t end2 = clock();

    double elapsed_secs2 = double(end2 - begin2) / CLOCKS_PER_SEC;

    //end Solve with armadillo.

    //

    vec lambda = A.diag();

    lambda = sort(lambda);
    vec first_three = lambda.subvec(2,4);

    //

    //Write results to output files.

    if (version.compare("version1") == 0) {
        ofile1.open("version1_output.txt",ios::out);
        output_version1(first_three, rho_max, h, iterations);
        ofile1.close();
    } else if (version.compare("version2") == 0) {
        ofile1.open("version2_output.txt",ios::out);
        output_version2(n_step,iterations, first_three, elapsed_secs1, elapsed_secs2);
        ofile1.close();
    } else if (version.compare("version3") == 0) {

        ofile1.open("version3_output.txt", ios::out);
        sortEigenVectors(A,R,n_step);
        normalize(R,n_step);
        output_version3(R,n_step);
        ofile1.close();
    }



    return 0;
}


double potential(double rho) {
    return rho*rho;
}

double coulomb_potential(double rho, double omega_r) {

    return omega_r*omega_r*rho*rho + 1.0/(rho);
}

double max_offdiag(mat M, int *k, int *l, int n) {

    /*
     * max_offdiag computes the off-diagonal norm off(A).
     * */

    //Page 219 fys3150 notes.

    double max = 0.0;

    for(int i = 1; i < n; i++) {
        for(int j = i+1; j < n; j++) {
            if(fabs(M(i,j)) > max) {
                max = fabs(M(i,j));
                *l = i;
                *k = j;
            }
        }
    }

    return max;

}

void JacobiRotation(mat &A, int k, int l, int n,mat &R) {
    double s, c;

    if(A(k,l) != 0.0) {
        double t, tau;
        tau = (A(l,l)-A(k,k))/(2*A(k,l));
        if(tau > 0) {
            t = 1.0/(tau + sqrt(1 + tau*tau));
        } else {
            t = -1.0/(-tau + sqrt(1.0 + tau*tau));
        }
        c = 1.0/sqrt(1 + t*t);
        s = c*t;
    } else {
        c = 1.0;
        s = 0.0;
    }

    double a_kk, a_ll, a_ik, a_il,r_ik, r_il;
    a_kk = A(k,k);
    a_ll = A(l,l);
    //changing matrix elements
    A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
    A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
    A(k,l) = 0.0;
    A(l,k) = 0.0;

    for(int i = 1; i < n; i++) {

        if(i != k && i != l) {
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = c*a_ik - s*a_il;
            A(k,i) = A(i,k);
            A(i,l) = c*a_il + s*a_ik;
            A(l,i) = A(i,l);
        }

        r_ik = R(i,k);
        r_il = R(i,l);
        R(i,k) = c*r_ik - s*r_il;
        R(i,l) = c*r_il + s*r_ik;

    }
    return;
}

void output_version1(vec first_three, double rho_max,double h, int iterations){

    ofile1 << setiosflags(ios::showpoint | ios::uppercase);

    for(int i = 0; i < 3; i++) {
        ofile1 << setw(15) << setprecision(8) << first_three(i);
    }
}

void output_version2(int n_step, int iterations, vec first_three, double time1, double time2) {
    ofile1 << setiosflags(ios::showpoint | ios::uppercase);
    ofile1 << setw(15) << setprecision(8) << n_step;
    ofile1 << setw(15) << setprecision(8) << iterations;
    ofile1 << setw(15) << setprecision(8) << first_three(0);
    ofile1 << setw(15) << setprecision(8) << first_three(1);
    ofile1 << setw(15) << setprecision(8) << first_three(2);
    ofile1 << setw(15) << setprecision(8) << time1;
    ofile1 << setw(15) << setprecision(8) << time2;

}

void normalize(mat &R, int n) {

    /*
     * normalize, normalizes the columns/eigenvectors in R
     * */

    vec norm_factors = zeros<vec>(n+1);

    for(int j=1; j < n; j++) {
        double norm = 0.0;
        for(int i=1; i < n; i++) {
            norm = norm + R(i,j)*R(i,j);
        }

        norm_factors(j) = sqrt(norm);
    }

    for(int j=1; j < n; j++) {

        for(int i=1; i < n; i++) {
            R(i,j) = R(i,j)/norm_factors(j);
        }

    }

}

void output_version3(mat R, int n) {

    ofile1 << setiosflags(ios::showpoint | ios::uppercase);



    for(int i = 0; i < n+1; i++) {
        ofile1 << setw(20) << setprecision(8) << R(i,1);
        ofile1 << setw(20) << setprecision(8) << R(i,2);
        ofile1 << setw(20) << setprecision(8) << R(i,3) << endl;

    }

}

void sortEigenVectors(mat &A, mat &R, int n) {

    /*
     * sortEigenVectors sorts the columns in R, using bubble-sort, such that the column 1 corresponds to lamba1
     * column 2 to lambda2, and so on...
     * */

    for(int i=1; i < n; i++) {
        for(int j=1; j < n-1; j++) {
            if(A(j,j) > A(j+1,j+1)) {
                //std::swap(A(j,j), A(j+1, j+1));

                double temp = A(j,j);
                A(j,j) = A(j+1,j+1);
                A(j+1,j+1) = temp;
                R.swap_cols(j,j+1);
            }
        }
    }
}
