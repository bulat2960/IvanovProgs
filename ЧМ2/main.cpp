#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>

using namespace std;

const bool out = false;

struct Matrix
{
    vector<vector<double>> matrix;

    Matrix(int x, int number = 0)
    {
        for (int i = 0; i < x; i++)
        {
            matrix.push_back(vector<double>(x, number));
        }
    }

    Matrix(vector<vector<double>> m)
    {
        matrix = m;
    }

    Matrix operator+(const Matrix& m)
    {
        Matrix result(matrix.size());
        for (unsigned i = 0; i < matrix.size(); i++)
        {
            for (unsigned j = 0; j < matrix[i].size(); j++)
            {
                result.matrix[i][j] = matrix[i][j] + m.matrix[i][j];
            }
        }
        return result;
    }

    Matrix operator*(double number)
    {
        Matrix result(matrix.size());
        for (unsigned i = 0; i < matrix.size(); i++)
        {
            for (unsigned j = 0; j < matrix[i].size(); j++)
            {
                result.matrix[i][j] = matrix[i][j] * number;
            }
        }
        return result;
    }

    void transpose()
    {
        for (unsigned i = 0; i < matrix.size(); i++)
        {
            for (unsigned j = i + 1; j < matrix[i].size(); j++)
            {
                swap(matrix[i][j], matrix[j][i]);
            }
        }
    }

    Matrix getNewMatrix(const Matrix& m)
    {
        Matrix result(matrix.size());
        for (unsigned i = 0; i < matrix.size(); i++)
        {
            for (unsigned j = 0; j < matrix[i].size(); j++)
            {
                int number = 0;
                for (unsigned k = 0; k < matrix[i].size(); k++)
                {
                    number += matrix[k][i] * matrix[k][j];
                }
                result.matrix[i][j] = number;
            }
        }
        return result;
    }

};

void print(const Matrix& A, const vector<double>& b)
{
    for (unsigned i = 0; i < A.matrix.size(); i++)
    {
        cout << "(x" << i + 1 << ")" << ' ';
        for (unsigned j = 0; j < A.matrix[i].size(); j++)
        {
            cout << setprecision(4) << setw(9) << setfill(' ') << A.matrix[i][j];
        }
        cout << " (x" << i + 1 << ")";
        cout << setprecision(4) << setw(10) << setfill(' ') << b[i] << endl;
    }
    cout << endl;
}

void printResult(const vector<double>& result)
{
    for (unsigned i = 0; i < result.size(); i++)
    {
        cout << "x" << i + 1 << " = " << result[i] << endl;
    }
    cout << endl;
}

int findCorrectString(const Matrix& A, int i)
{
    for (unsigned j = i + 1; j < A.matrix.size(); j++)
    {
        if (A.matrix[j][i] != 0)
        {
            return j;
        }
    }
    return -1;
}

void gauss1(Matrix& A, vector<double>& b)
{
    double coef = 0.0;

    if (out) cout << "INITIAL SYSTEM:" << endl;
    for (unsigned i = 0; i < A.matrix.size(); i++)
    {
        if (out) print(A, b);
        double x = A.matrix[i][i];
        b[i] /= x;
        for (unsigned j = 0; j < A.matrix[i].size(); j++)
        {
            A.matrix[i][j] /= x;
        }

        for (unsigned j = i + 1; j < A.matrix[i].size(); j++)
        {
            if (A.matrix[i][i] == 0)
            {
                int x = findCorrectString(A, i);
                if (x == -1)
                {
                    break;
                }

                for (unsigned k = 0; k < A.matrix[i].size(); k++)
                {
                    swap(A.matrix[i][k], A.matrix[x][k]);
                }
            }
            coef = A.matrix[j][i];
            for (unsigned k = 0; k < A.matrix[i].size(); k++)
            {
                A.matrix[j][k] -= A.matrix[i][k] * coef;
                if (A.matrix[j][k] > -0.0001 && A.matrix[j][k] < 0.0001)
                {
                    A.matrix[j][k] = 0;
                }
            }
            b[j] -= b[i] * coef;
        }
    }
}

vector<double> gauss2(const Matrix& A, const vector<double>& b)
{
    vector<double> result(b.size());
    for (int i = b.size() - 1; i >= 0; i--)
    {
        if (i == (int)(b.size() - 1))
        {
            result[i] = b.back();
        }
        else
        {
            double x = b[i];
            for (unsigned j = i + 1; j < A.matrix[i].size(); j++)
            {
                x -= A.matrix[i][j] * result[j];
            }
            result[i] = x;
        }
    }
    return result;
}

void check(const Matrix B, vector<double> result, const vector<double> b)
{
    for (unsigned i = 0; i < B.matrix.size(); i++)
    {
        double x = 0;
        for (unsigned j = 0; j < B.matrix[i].size(); j++)
        {
            cout << setw(5) << setfill(' ') << B.matrix[i][j] << "*" << result[j];
            if (j < B.matrix[i].size() - 1) cout << " + ";
            x += B.matrix[i][j] * result[j];
        }
        cout << " = " << setw(5) << x << setw(10) << b[i] << endl;
    }
}

vector<double> zeidel(const Matrix& A, const vector<double>& b)
{
    Matrix C = A;
    vector<double> d = b;
    for (unsigned i = 0; i < d.size(); i++)
    {
        for (unsigned j = 0; j < d.size(); j++)
        {
            if (i != j)
            {
                C.matrix[i][j] = -1 * A.matrix[i][j] / A.matrix[i][i];
            }
            else
            {
                C.matrix[i][j] = 0;
            }
        }
        d[i] /= A.matrix[i][i];
    }
    if (out) print(C, d);


    vector<double> x = d;

    unsigned ITER_LIMIT = 100000;
    for (unsigned i = 0; i < ITER_LIMIT; i++)
    {
        vector<double> x_new(x.size());
        for (unsigned j = 0; j < x.size(); j++)
        {
            x_new[j] += d[j];
            for (unsigned k = 0; k < x.size(); k++)
            {
                if (k > j)
                {
                    x_new[j] += x[k] * C.matrix[j][k];
                }
                if (k < j)
                {
                    x_new[j] += x_new[k] * C.matrix[j][k];
                }
            }
        }

        unsigned counter = 0;
        for (unsigned i = 0; i < x.size(); i++)
        {
            if (fabs(x[i] - x_new[i]) < pow(10, -3))
            {
                counter++;
            }
        }
        if (counter == x.size())
        {
            return x;
        }

        x = x_new;
    }
    return vector<double>(0);
}

vector<double> getTransposedVector(const Matrix& m, const vector<double>& vec)
{
    vector<double> result(vec.size());
    for (unsigned i = 0; i < m.matrix.size(); i++)
    {
        int number = 0;
        for (unsigned j = 0; j < m.matrix[i].size(); j++)
        {
            number += m.matrix[i][j] * vec[j];
        }
        result[i] = number;
    }
    return result;
}

int main()
{
    const Matrix A({{12, 1, 1}, {1, 14, 1}, {1, 1, 16}});
    const vector<double> b = {14, 16, 18};

    cout << "-----------------GAUSS METHOD----------------" << endl;
    Matrix A_gauss = A;
    vector<double> b_gauss = b;

    gauss1(A_gauss, b_gauss);
    if (out) print(A_gauss, b_gauss);

    vector<double> result1 = gauss2(A_gauss, b_gauss);
    cout << "RESULT VECTOR:" << endl;
    printResult(result1);

    if (out) check(A, result1, b);
    cout << endl;

    cout << "-----------------ZEIDEL METHOD----------------" << endl;
    if (out) cout << "Matrix C and vector D" << endl;
    vector<double> result2 = zeidel(A, b);
    cout << "RESULT VECTOR:" << endl;
    printResult(result2);

    int k;
    cin >> k;
    Matrix A1(k, 0);
    Matrix A2(k, 1);

    double e = pow(10, -3);

    vector<double> bEpsilon(k, -1);
    bEpsilon.back() = 1 + e;

    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < k; j++)
        {

            if (j >= i)
            {
                A1.matrix[i][j] = 1;
            }
            if (j > i)
            {
                A1.matrix[i][j] = -1;
                A2.matrix[i][j] = -1;
            }
        }
    }

    Matrix AEpsilon = A1 + A2 * (pow(10, -3) * 10);

    cout << "-----------------GAUSS METHOD----------------" << endl;
    const Matrix AEpsilonOrig = AEpsilon;
    const vector<double> bEpsilonOrig = bEpsilon;

    gauss1(AEpsilon, bEpsilon);
    if (out) print(AEpsilon, bEpsilon);

    vector<double> result3 = gauss2(AEpsilon, bEpsilon);
    cout << "RESULT VECTOR:" << endl;
    printResult(result3);

    if (out) check(AEpsilonOrig, result3, bEpsilonOrig);
    cout << endl;

    cout << "-----------------ZEIDEL METHOD----------------" << endl;
    Matrix AEpsilonTransposed = AEpsilonOrig;
    Matrix temp = AEpsilonOrig;
    AEpsilonTransposed.transpose();

    Matrix newMatrix = temp.getNewMatrix(AEpsilonTransposed);
    vector<double> newVector = getTransposedVector(AEpsilonTransposed, bEpsilonOrig);

    if (out) cout << "Matrix A and vector B - AFTER TRANSPOSE:" << endl;
    if (out) print(newMatrix, newVector);

    if (out) cout << "Matrix C and vector D:" << endl;
    vector<double> result4 = zeidel(newMatrix, newVector);
    cout << "RESULT VECTOR:" << endl;
    printResult(result4);



    return 0;
}
