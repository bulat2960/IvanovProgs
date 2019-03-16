#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

const double eps = pow(10, -4);
const double c1 = 0.5;
const double c2 = 0.3;
const double pi = 3.1415926;

namespace math
{
    int fact(int x)
    {
        if (x == 1) return 1;
        return x * fact(x - 1);
    }

    double sin(double x)
    {
        double eps1 = pow(10, -3) / (3 * c1);
        int k = 0;
        double pAnswer = 1;
        double cAnswer = 0;
        while (fabs(pAnswer - cAnswer) > eps1)
        {
            pAnswer = cAnswer;
            cAnswer += (pow(-1, k) * pow(x, 2 * k + 1)) / fact(2 * k + 1);
            k += 1;
        }
        return cAnswer;
    }

    double atan(double x)
    {
        double eps2 = pow(10, -3) / (3 * c2);
        int k = 0;
        double pAnswer = 1;
        double cAnswer = 0;
        while (fabs(pAnswer - cAnswer) > eps2)
        {
            pAnswer = cAnswer;
            cAnswer += (pow(-1, k) * pow(x, 2 * k + 1)) / (2 * k + 1);
            k += 1;
        }
        return cAnswer;
    }

    double sqrt(double x)
    {
        double eps3 = pow(10, -3) / 3;
        int k = 1;
        double pAnswer = x;
        double cAnswer = x + eps3 + 0.1;
        while (fabs(pAnswer - cAnswer) > eps3)
        {
            pAnswer = cAnswer;
            cAnswer = 0.5 * (pAnswer + x / pAnswer);
            k += 1;
        }
        return cAnswer;
    }

    double z(double x)
    {
        return math::sin(pi / 2 - 2.8 * x - math::sqrt(1 + x)) * math::atan(1.5 * x + 0.2);
    }
}

double z(double x)
{
    return sin(pi / 2 - 2.8 * x - sqrt(1 + x)) * atan(1.5 * x + 0.2);
}

int main()
{
    double x1 = 0.1;
    double step = 0.01;
    double x2 = 0.2;

    while (x1 <= x2 + 0.0001)
    {
        double x = x1;
        double z1 = math::z(x);
        double z2 = z(x);
        double zAbs = fabs(z1 - z2);
        cout << fixed << setprecision(2) << setw(5)
             << setfill(' ')    << x  << ' '
             << setw(15) << setprecision(4) << z1 << ' '
             << setw(15) << setprecision(4) << z2 << ' ';
        cout << setw(15) << scientific << setprecision(4) << zAbs << endl;
        x1 += step;
    }
    return 0;
}
