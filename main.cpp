#include <iostream>
#include <cmath>

using namespace std;

double source_term(double);
double exact_solution(double);
void solve(double*,double*,double*,double*,double*,int);
void rel_error(double*,double*,double*,int);

int main()
{
    int i;
    double *a;
    double *b;
    double *c;
    double *g;
    double *computed_solution;
    double *exact;
    double *relative_error;
    int n;

    cout << "Give dimension size (nxn): " << endl;
    cin >> n;

    double h = 1.0/(n+1);

    a = new double[n+2];
    b = new double[n+2];
    c = new double[n+2];
    g =  new double[n+2];
    relative_error = new double[n+2];
    computed_solution = new double[n+2];
    exact = new double[n+2];

    for(i=0; i < n+2; i++) {
        exact[i] = exact_solution(i*h);
    }


    for(i=0; i < n+2; i++) {
        g[i] = pow(h,2)*source_term(i*h);
        a[i] = 2.0;
        b[i] = -1.0;
        c[i] = -1.0;
        relative_error[i] = 0.0;
    }

    computed_solution[0] = 0.0;
    computed_solution[n+1] = 0.0;

    solve(a,b,c,g,computed_solution,n);
    rel_error(relative_error, exact, computed_solution,n);

    cout << relative_error[6] << endl;

    return 0;
}

double source_term(double x) {
    return 100*exp(-10*x);
}

double exact_solution(double x) {
    return 1.0-(1.0-exp(-10.0))*x - exp(-10*x);
}

void solve(double *a, double *b, double *c, double *g,double *sol,int n) {

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
