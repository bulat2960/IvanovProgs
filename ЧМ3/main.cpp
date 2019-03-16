#include <iostream>
#include <cmath>

using namespace std;

double f(double x)
{
    return pow(x, 2) + 4 * sin(x) - 2;
}

double df_dx(double x)
{
    return 2 * x + 4 * cos(x);
}

double newton(double x0)
{
    while (true)
    {
        cout << "x0 = " << x0 << endl;
        double x1 = x0 - (f(x0) / df_dx(x0));
        if (abs(x1 - x0) < pow(10, -5))
        {
            return x1;
        }
        x0 = x1;
    }
    return 0;
}

double det(double a11, double a12, double a21, double a22)
{
    return a11 * a22 - a12 * a21;
}

double f1(double x, double y)
{
    return sin(y - 1) + x - 1.3;
}

double f2(double x, double y)
{
    return y - sin(x + 1) - 0.8;
}

double df1_dx(double x, double y)
{
    return 1;
}

double df1_dy(double x, double y)
{
    return cos(y - 1);
}

double df2_dx(double x, double y)
{
    return -cos(x + 1);
}

double df2_dy(double x, double y)
{
    return 1;
}

pair<double, double> newton2(double x0, double y0)
{
	double x1 = x0;
	double y1 = y0;
	while (true)
    {
		cout << x1 << ' ' << y1 << endl;
		double dx_num = det(df1_dy(x1,y1), f1(x1,y1), df2_dy(x1,y1), f2(x1,y1));
		double dx_den = det(df1_dx(x1,y1), df1_dy(x1,y1), df2_dx(x1,y1), df2_dy(x1,y1));
		double dx = dx_num / dx_den;
		double dy_num = det(f1(x1,y1), df1_dx(x1,y1), f2(x1,y1), df2_dx(x1,y1));
		double dy_den = det(df1_dx(x1,y1), df1_dy(x1,y1), df2_dx(x1,y1), df2_dy(x1,y1));
		double dy = dy_num / dy_den;
		if (abs(dx) < pow(10, -5) && abs(dy) < pow(10, -5))
        {
			return pair<double, double>(x1, y1);
        }
		x1 += dx;
		y1 += dy;
    }
}

int main()
{
    cout << "NEWTON1" << endl;
    double x1 = newton(0.5);
    double x2 = newton(-2);

    cout << "Result:" << endl;
    cout << "x1 = " << x1 << endl;
    cout << "x2 = " << x2 << endl;

    cout << "f:" << endl;
    cout << f(x1) << endl;
    cout << f(x2) << endl;

    cout << endl;

    cout << "NEWTON2" << endl;
    pair<double, double> p = newton2(0.5, 2);

    cout << "Result:" << endl;
    cout << p.first << ' ' << p.second << endl;

    cout << "f:" << endl;
    cout << f1(p.first, p.second) << endl;
    cout << f2(p.first, p.second) << endl;
    return 0;
}
