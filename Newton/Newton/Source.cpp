#include <iostream>
#include <cmath>
#include <vector>
using namespace std;
void gaussElimination(vector<vector<double>> A, vector<double> B, vector<double>& x) {
    for (int i = 0; i < 2; i++)
    {
        B[i] *= -1;
    }
    for (int i = 0; i < 2; i++) {
        int maxRow = i;
        for (int k = i + 1; k < 2; k++) {
            if (abs(A[k][i]) > abs(A[maxRow][i])) {
                maxRow = k;
            }
        }
        swap(A[maxRow], A[i]);
        swap(B[maxRow], B[i]);
        for (int k = i + 1; k < 2; k++) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < 2; j++) {
                A[k][j] -= factor * A[i][j];
            }
            B[k] -= factor * B[i];
        }
    }
    for (int i = 0; i < 2; i++)
    {
        x.push_back(B[i]);
    }
    for (int i = 1; i >= 0; i--) {
        for (int j = i + 1; j < 2; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }
}
vector<double> f(vector<double> x) {
    vector<double> F;
    F.push_back(2 * x[0] * x[0] - x[0] * x[1] - 5 * x[0] + 1);
    F.push_back(x[0] + 3 * log10(x[0]) - x[1] * x[1]);
    return F;
}
void jac(double x1, double x2, vector<vector<double>>& jacobi) {
    jacobi[0][0] = 4 * x1 - x2 - 5;
    jacobi[0][1] = -x1;
    jacobi[1][0] = 1 + 3 / (x1 * log(10));
    jacobi[1][1] = -2 * x2;
}
void anotherJac(vector<double> x, vector<vector<double>>& jacobi, double m)
{
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            vector<double> xPlus = x;
            xPlus[j] += xPlus[j] * m;
            jacobi[i][j] = (f(xPlus)[i] - f(x)[i]) / (x[j] * m);
        }
    }
}
void newtonWithMyJacobi(double x1, double x2, double eps = 1e-9, int max_iter = 100) {
    vector<double> x = { x1,x2 };
    cout << "Начальное приближение: " << x1 << ", " << x2 << endl;
    cout << "ε1 = " << eps << "; ε2 = " << eps << "; max_iter = " << max_iter << endl;
    cout << "| Номер итерации | 1 | 2 |" << endl;
    int k;
    for (k = 0; k < max_iter; k++) {
        vector<vector<double>> jacobi = { {0,0},{0,0} };
        vector<double> F = f(x);
        vector<double> dx;
        jac(x[0], x[1], jacobi);
        anotherJac(x, jacobi, 0.1);
        gaussElimination(jacobi, F, dx);
        double maxF;
        double maxGap = 0;
        if (abs(F[0]) > abs(F[1]))
        {
            maxF = abs(F[0]);
        }
        else
        {
            maxF = abs(F[1]);
        }
        for (int i = 0; i < 2; i++)
        {
            double gap;
            if (abs(x[i] + dx[i]) < 1) { gap = abs(dx[i]); }
            else { gap = abs((dx[i]) / x[i] + dx[i]); }
            if (maxGap < gap) { maxGap = gap; }
        }
        for (int i = 0; i < 2; i++)
        {
            x[i] += dx[i];
        }
        cout << k << " " << maxF << " " << maxGap << endl;
        if (maxF <= eps || maxGap <= eps)
        {
            cout << "Решение: " << x[0] << "; " << x[1] << endl;
            break;
        }
    }
    if (k == max_iter)
    {
        cout << "IER = 2";
    }
}
void newtonWithGeneralJacobi(double x1, double x2, double m, double eps = 1e-9, int max_iter = 100) {
    vector<double> x = { x1,x2 };
    cout << "Начальное приближение: " << x1 << ", " << x2 << endl;
    cout << "ε1 = " << eps << "; ε2 = " << eps << "; max_iter = " << max_iter << endl;
    cout << "| Номер итерации | 1 | 2 |" << endl;
    int k;
    for (k = 0; k < max_iter; k++) {
        vector<vector<double>> jacobi = { {0,0},{0,0} };
        vector<double> F = f(x);
        vector<double> dx;
        anotherJac(x, jacobi, m);
        gaussElimination(jacobi, F, dx);
        double maxF;
        double maxGap = 0;
        if (abs(F[0]) > abs(F[1]))
        {
            maxF = abs(F[0]);
        }
        else
        {
            maxF = abs(F[1]);
        }
        for (int i = 0; i < 2; i++)
        {
            double gap;
            if (abs(x[i] + dx[i]) < 1) { gap = abs(dx[i]); }
            else { gap = abs((dx[i]) / x[i] + dx[i]); }
            if (maxGap < gap) { maxGap = gap; }
        }
        for (int i = 0; i < 2; i++)
        {
            x[i] += dx[i];
        }
        cout << k << " " << maxF << " " << maxGap << endl;
        if (maxF <= eps || maxGap <= eps)
        {
            cout << "Решение: " << x[0] << "; " << x[1] << endl;
            break;
        }
    }
    if (k == max_iter)
    {
        cout << "IER = 2";
    }
}
int main() {
    setlocale(LC_ALL, "RUS");
    newtonWithMyJacobi(3, 2);
    cout << endl;
    newtonWithMyJacobi(3, -2);
    cout << "Another jacobi" << endl;
    cout << "M = 0.01" << endl;
    newtonWithGeneralJacobi(3, 2, 0.01);
    cout << endl;
    newtonWithGeneralJacobi(3, -2, 0.01);
    cout << "M = 0.05" << endl;
    newtonWithGeneralJacobi(3, 2, 0.05);
    cout << endl;
    newtonWithGeneralJacobi(3, -2, 0.05);
    cout << "M = 0.1" << endl;
    newtonWithGeneralJacobi(3, 2, 0.1);
    cout << endl;
    newtonWithGeneralJacobi(3, -2, 0.1);
}