#ifndef QUARTIC_EQUATION_H
#define QUARTIC_EQUATION_H

//求解一元四次方程 a*x^4+b*x^3+c*x^2+d*x+e=0；
//输入任意系数a,b,c,d,e; 输出四个解分别为x1,x2,x3,x4；

#include  <complex>
#include  <cmath>

using namespace std;

struct Quartic_Equation {
    complex<double>  x1, x2, x3, x4;

    Quartic_Equation() = default;
    ~Quartic_Equation() = default;

    void quadratic_equation(double, double, double, complex<double> &, complex<double> &) const;
    double cubic_equation(double, double, double, double);
    void quartic_equation(double, double, double, double, double);
};

// 1. 三次方程、四次方程中所调用的二次方程；
inline
void Quartic_Equation::quadratic_equation(double a4, double b4, double c4, complex<double> &z1, complex<double> &z2) const {

    double delta;
    complex<double> temp1, temp2;

    delta = b4*b4 - 4 * a4*c4;

    complex<double> temp(delta, 0);

    temp1 = (-b4) / (2 * a4);
    temp2 = sqrt(temp) / (2 * a4);

    z1 = temp1 + temp2;
    z2 = temp1 - temp2;
}

// 2. 四次方程中所调用的三次方程；
inline
double Quartic_Equation::cubic_equation(double a3, double b3, double c3, double d3) {

    double p, q, delta;
    double M, N;
    double y0;

    complex<double> temp1, temp2;
    complex<double> y1, y2;

    if (a3 == 0) {
        quadratic_equation(b3, c3, d3, y1, y2);

        x1 = y1;
        x2 = y2;
        x3 = 0;
    }
    else {
        p = -1.0 / 3 * pow((b3*1.0 / a3), 2.0) + c3*1.0 / a3;
        q = 2.0 / 27 * pow((b3*1.0 / a3), 3.0) - 1.0 / 3 * b3*c3 / (a3*a3) + d3*1.0 / a3;

        delta = pow((q / 2.0), 2.0) + pow((p / 3.0), 3.0);

        complex<double>  omega1(-1.0 / 2, sqrt(3.0) / 2.0);
        complex<double>  omega2(-1.0 / 2, -sqrt(3.0) / 2.0);

        complex<double>  yy(b3 / (3.0*a3), 0.0);

        M = -q / 2.0;

        if (delta<0) {

            N = sqrt(fabs(delta));
            complex<double>  s1(M, N);
            complex<double>  s2(M, -N);

            x1 = (pow(s1, (1.0 / 3)) + pow(s2, (1.0 / 3))) - yy;
            x2 = (pow(s1, (1.0 / 3))*omega1 + pow(s2, (1.0 / 3))*omega2) - yy;
            x3 = (pow(s1, (1.0 / 3))*omega2 + pow(s2, (1.0 / 3))*omega1) - yy;

        }
        else {
            N = sqrt(delta);

            complex<double>  f1(M + N, 0);
            complex<double>  f2(M - N, 0);

            if (M + N >= 0)
                temp1 = pow((f1), 1.0 / 3);
            else
                temp1 = -norm(pow(sqrt(f1), 1.0 / 3));


            if (M - N >= 0)
                temp2 = pow((f2), 1.0 / 3);
            else
                temp2 = -norm(pow(sqrt(f2), 1.0 / 3));


            x1 = temp1 + temp2 - yy;
            x2 = omega1*temp1 + omega2*temp2 - yy;
            x3 = omega2*temp1 + omega1*temp2 - yy;

        }

    }

    // y0 为所调用的三次方程的返回值；
    y0 = real(x1);

    return y0;
}

// 3.main函数所调用的四次方程；
inline
void Quartic_Equation::quartic_equation(double a, double b, double c, double d, double e) {

    double a2, b2, c2;
    double a4, b4, c4;
    double a3, b3, c3, d3;
    double y;
    complex<double>  y1, y2, y3, y4;

    if (b == 0 && c == 0 && d == 0 && e == 0) {
        x1 = 0; x2 = 0; x3 = 0; x4 = 0;
    }
    else if (b == 0 && d == 0 && e == 0) {
        b3 = a; c3 = 0; d3 = c;
        quadratic_equation(b3, c3, d3, y1, y2);
        x1 = y1; x2 = y2;  x3 = 0; x4 = 0;
    }
    else {
        //把任意系数的四次方程化为首项系数为1的四次方程；
        b = b / a;  c = c / a;  d = d / a;  e = e / a;

        //所调用的三次方程的系数;
        a3 = 8.0;
        b3 = -4.0*c;
        c3 = 2.0*b*d - 8.0*e;
        d3 = e*(4.0*c - b*b) - d*d;

        y = cubic_equation(a3, b3, c3, d3);
        //把三次方程的返回值赋给 y ；

        //第一次调用的二次方程的系数；
        a2 = 1.0;
        b2 = b / 2.0 - sqrt(8.0*y + b*b - 4 * c) / 2.0;
        c2 = y - (b*y - d) / sqrt(8.0*y + b*b - 4 * c);
        quadratic_equation(a2, b2, c2, y1, y2);

        x1 = y1;
        x2 = y2;

        //第二次调用的二次方程的系数；
        a4 = 1.0;
        b4 = b / 2.0 + sqrt(8.0*y + b*b - 4.0*c) / 2.0;
        c4 = y + (b*y - d) / sqrt(8.0*y + b*b - 4.0*c);
        quadratic_equation(a4, b4, c4, y3, y4);

        x3 = y3;
        x4 = y4;
    }
}

#endif // QUARTIC_EQUATION_H
