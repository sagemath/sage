/* sage_setup: distribution = sagemath-lcalc
 */
#include "lcalc/L.h"
int *new_ints(int l)
{
    return new int[l];
}


void del_ints(int *A)
{
    delete[] A;
}


double *new_doubles(int l)
{
    return new double[l];
}


void del_doubles(double *A)
{
    delete[] A;
}


Complex *new_Complexes(int l)
{
    return new Complex[l];
}


void del_Complexes(Complex *A)
{
    delete[] A;
}

Complex new_Complex(double r, double i)
{
    return  Complex(r,i);
}


void testL(L_function<Complex> *L)
{
    int i;
    cout << "number of coefficients " << L->number_of_dirichlet_coefficients << endl;
    cout << "Dirichlet coeffs"<< endl;
    for (i=0;i< min(30, L->number_of_dirichlet_coefficients +1); i++)
        cout << L->dirichlet_coefficient[i]<<endl;
    cout << "Q " << L->Q << endl;
    cout << "Omega " << L->OMEGA << endl;
    cout << "a " << L->a << endl;
    cout << "Period " << L->period << endl;
    cout << "Number of Poles " << L->number_of_poles << endl;
    cout << "What type " << L->what_type_L << endl;
    for (i=0;i< L->number_of_poles+1;i++)
    {
        cout<< "pole[" << i << "] =  " << L->pole[i] << endl;
        cout<< "residue[" << i << "] =  " << L->residue[i] << endl;
    }
    cout << "Value at .5 " << L->value(.5) <<endl;
    cout << "Value at 1"  << L->value(1.0) <<endl;
    cout << "Value at .5+I"  << L->value(.5+I) <<endl;
}
