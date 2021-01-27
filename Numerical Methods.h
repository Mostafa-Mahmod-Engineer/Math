#pragma once
#include"Polynomial.h"
#define fileName "Numerical Methods.h"


//******************************************** DEFFERENCIATION *********************************************************
///===================  FORWARD,BACKWARD DIFFERENCE METHODS:
double forward(double(*f)(double), double x, double h)
{
	return (f(x + h) - f(x)) / h;
}
double backward(double(*f)(double), double x, double h)
{
	return (f(x) - f(x - h)) / h;
}
double central(double(*f)(double), double x, double h)
{
	return (f(x + h) - f(x - h)) / (2 * h);
}
double forward2(double(*f)(double), double x, double h)
{
	return (f(x) - 2 * f(x + h) + f(x + 2 * h)) / (h*h);
}
double backward2(double(*f)(double), double x, double h)
{
	return (f(x) - 2 * f(x - h) + f(x - 2 * h)) / (h*h);
}
double central2(double(*f)(double), double x, double h)
{
	return (f(x + h) - 2 * f(x) + f(x - h)) / (h*h);
}
/// n-th derivative, with m-th order of accuracy
double derive(double(*f)(double x), double x, int n = 1, int m = 0, double h = 0.004)
{
	if (n < 1)
		error("An Invalide Derivative Order Has Been Entered\nError Call from derive()");
	if (m == 0)
		m = n + n % 2;
	if (m < 2)
		error("An Invalide Derivative Accuracy Order Has Been Entered\nError Call from derive()");
	if (m % 2)
		error("Accuracy Order Can't Be Odd\nError Call from derive()");

	int p = std::floor((n + 1) / 2) - 1 + m / 2;

	real::mat P(2 * p + 1);
	P._setrows(0, 0, 1);
	for (register int i = 1; i <= 2 * p; i++)
	{
		register int k = -p;
		for (register int j = 0; j <= 2 * p; j++, k++)
			P.setelement(i, j, std::pow(k, i));
	}
	real::vec coefficients(2 * p + 1);
	coefficients.set(n, Factorial(n));
	coefficients = (!P)*coefficients;
	double res = 0;
	for (register int i = -p; i <= p; i++)
	{
		res += coefficients[i + p] * f(x + h * i);
	}
	return res / std::pow(h, n);
}

///==================================================== PARTIAL DEFFERENCIATION =================================================================
double rk4(double(*f)(double, double), double x0, double y0, double x1, double h)
{
	int n = std::round(std::abs(x1 - x0) / h);
	double k1, k2, k3, k4, kav;
	for (register int i = 0; i < n; i++)
	{
		k1 = h * f(x0, y0);
		k2 = h * f(x0 + h / 2, y0 + k1 / 2);
		k3 = h * f(x0 + h / 2, y0 + k2 / 2);
		k4 = h * f(x0 + h, y0 + k3);
		kav = (k1 + 2 * k2 + 2 * k3 + k4) / 6;
		x0 += h;
		y0 += kav;
	}
	return y0;
}
double rk2(double(*f)(double, double), double x0, double y0, double x1, double h)
{
	int n = std::round(std::abs(x1 - x0) / h);
	double k1, k2, kav;
	for (register int i = 0; i < n; i++)
	{
		k1 = h * f(x0, y0);
		k2 = h * f(x0 + h, y0 + k1);
		kav = (k1 + k2) / 2;
		x0 += h;
		y0 += kav;
	}
	return y0;
}
pair<double, double> rk2multi(double(*f)(double, double, double), double(*g)(double, double, double),
	double x0, double y0, double z0, double x1, double h)
{
	int n = std::round(std::abs(x1 - x0) / h);
	double ky1, ky2, kyav, kz1, kz2, kzav;
	for (register int i = 0; i < n; i++)
	{
		ky1 = h * f(x0, y0, z0);
		kz1 = h * g(x0, y0, z0);
		ky2 = h * f(x0 + h, y0 + ky1, z0 + kz1);
		kz2 = h * g(x0 + h, y0 + ky1, z0 + kz1);
		kyav = (ky1 + ky2) / 2;
		kzav = (kz1 + kz2) / 2;
		x0 += h;
		y0 += kyav;
		z0 += kzav;
	}
	return pair<double, double>(y0, z0);
}

double partialx(double(*f)(double x, double y), double x, double y, double h = 0.0001) { return (f(x + h, y) - f(x - h, y)) / (2 * h); }

double partialy(double(*f)(double x, double y), double x, double y, double h = 0.0001) { return (f(x, y + h) - f(x, y - h)) / (2 * h); }

double partialxy(double(*f)(double x, double y), double x, double y, double hx = 0.0001, double hy = 0.0001)
{
	return (f(x + hx, y + hy) - f(x + hx, y - hy) - f(x - hx, y + hy) + f(x - hx, y - hy)) / (4 * hx*hy);
}


//*********************************************************** ROOT FINDING *********************************************************************
///====================================== AUXILIARY FUNCTIONS ====================================
void findRootP(double(*f)(double x), double& a, double& b)
{
	a = 0;
	short sign = abs(f(0)) / f(0);
	while (abs(f(a)) / f(a) == sign)
		if (a == 100000000)
			printf("No Positive Root Exists");
		else
			a++;
	b = a;
	a--;
}
void findRootN(double(*f)(double x), double& a, double& b)
{
	a = 0;
	short sign = abs(f(0)) / f(0);
	while (abs(f(a)) / f(a) == sign)
		if (a == -100000000)
			printf("No Negative Root Exists");
		else
			a--;
	b = a;
	a++;
}
void findRootP(double(*f)(double x), double& a, double& b, double a0, double step = 1, double limit = 1000000)
{
	a = a0;
	short sign = abs(f(a0)) / f(a0);
	while (abs(f(a)) / f(a) == sign)
		if (a == limit)
			printf("No Positive Root Exists");
		else
			a += step;
	b = a;
	a -= step;
}
void findRootN(double(*f)(double x), double& a, double& b, double a0, double step = 1, double limit = -1000000)
{
	a = a0;
	short sign = abs(f(a0)) / f(a0);
	while (abs(f(a)) / f(a) == sign)
		if (a == limit)
			printf("No Negative Root Exists");
		else
			a -= step;
	b = a;
	a += step;
}
///=============================================  NON-LINEAR EQUATIONS  ============================================================
double FPI(double(*phai)(double x), double x0, int n, int k)
{
	if (derive(phai, x0) > 1)
	{
		printf("This Function Doesnot Converge");
		system("pause");
	}
	double x = x0;
	double e = 1;
	int counter = 0;
	if (n == 0)
		while (e > 0.5*pow(10, -k))
		{
			x = phai(x0);
			counter++;
			e = abs(x - x0);
			printf("x%d = %f\n", counter, x);
			x0 = x;
		}
	if (k == 0)
		for (register int i = 1; i <= n; i++)
		{
			x = phai(x0);
			printf("x%d = %f\n", i, x);
			x0 = x;
		}
	return x;
}

double bisection(double(*f)(double x), double a, double b, int n, int k)
{
	double x = (a + b) / 2;
	double x0;
	if (a > b)
		std::swap(a, b);
	//   a is the lower-bound <==> b is the upper-bound
	double e = 1;
	array<double> A(10), B(10), X(10);
	A.push_back(a);
	B.push_back(b);
	X.push_back(x);
	if (n == 0)
		while (e > 0.5*pow(10, -k))
		{
			x0 = x;
			if (f(a)*f(x) > 0)
				a = x;
			else
				b = x;
			x = (a + b) / 2;
			A.push_back(a);
			B.push_back(b);
			X.push_back(x);
			e = abs((x - x0) / x);
		}
	if (k == 0)
		for (register int i = 1; i <= n; i++)
		{
			x0 = x;
			if (a*f(x) > 0)
				a = x;
			else
				b = x;
			x = (a + b) / 2;
			A.push_back(a);
			B.push_back(b);
			X.push_back(x);
		}
	std::string stra = "", strb = "", strx = "", line = "";
	array<std::string> lines(A.Size());
	array<std::string> aCol(A.Size()), bCol(A.Size()), xCol(A.Size());

	for (register int i = 0; i < A.Size(); i++)
	{
		stra = std::to_string(A[i]);
		strb = std::to_string(B[i]);
		strx = std::to_string(X[i]);
		for (register int j = stra.size() - 1; j >= 0; j--)
			if (A[i] == floor(A[i]))
			{
				stra = std::to_string(int(A[i]));
				break;
			}
			else if (stra[j] == '0' && stra.find('.') != -1 || stra[j] == '.')
				stra.pop_back();
			else
				break;
		for (register int j = strx.size() - 1; j >= 0; j--)
			if (X[i] == floor(X[i]))
			{
				strx = std::to_string(int(X[i]));
				break;
			}
			else if (strx[j] == '0' && strx.find('.') != -1 || strx[j] == '.')
				strx.pop_back();
			else
				break;
		for (register int j = strb.size() - 1; j >= 0; j--)
			if (B[i] == floor(B[i]))
			{
				strb = std::to_string(int(B[i]));
				break;
			}
			else if (strb[j] == '0'&&strb.find('.') != -1 || strb[j] == '.')
				strb.pop_back();
			else
				break;
		if (stra.size() > 11)
			stra = stra.substr(0, 10);
		if (strx.size() > 11)
			strx = strx.substr(0, 10);
		if (strb.size() > 11)
			strb = strb.substr(0, 10);
		aCol.push_back(stra);
		xCol.push_back(strx);
		bCol.push_back(strb);
	}
	line += "N";
	for (register int i = 0; i < int(log10(A.Size())) + 1; i++)
		line += " ";
	line += "|";
	for (register int i = 0; i < 11; i++)
		if (i == 5)
			line += "a";
		else
			line += " ";
	line += "|";
	for (register int i = 0; i < 11; i++)
		if (i == 5)
			line += "x";
		else
			line += " ";
	line += "|     b\n";
	lines.push_back(line);
	for (register int i = 1; i <= A.Size(); i++)
	{
		line = std::to_string(i);
		for (register int j = 0; j < int(log10(A.Size())) - int(log10(i)); j++)
			line += " ";
		line += " |";
		int a1 = (11 - aCol[i - 1].size()) / 2, a2 = ceil((11 - aCol[i - 1].size()) / 2.0);
		int x1 = (11 - xCol[i - 1].size()) / 2, x2 = ceil((11 - xCol[i - 1].size()) / 2.0);
		int b1 = (11 - bCol[i - 1].size()) / 2;
		for (register int j = 0; j < a1; j++)
			line += " ";
		line += aCol[i - 1];
		for (register int j = 0; j < a2; j++)
			line += " ";
		line += "|";
		for (register int j = 0; j < x1; j++)
			line += " ";
		line += xCol[i - 1];
		for (register int j = 0; j < x2; j++)
			line += " ";
		line += "|";
		for (register int j = 0; j < b1; j++)
			line += " ";
		line += bCol[i - 1];
		line += "\n";
		lines.push_back(line);
	}
	for (register int i = 0; i < lines.Size(); i++)
		std::cout << lines[i];
	for (register int i = 0; i < lines[lines.Size() - 1].size() - 1; i++)
		std::cout << "=";
	std::cout << "|\n";
	return x;
}

double falsePosition(double(*f)(double x), double a, double b, int n, int k)
{
	auto F = [&a, &b, f] {return b + (f(b)*(a - b)) / (f(b) - f(a)); };
	if (a > b)
		std::swap(a, b);
	double x = F();
	double x0, e = 1;
	array<double> A(10), B(10), X(10);
	A.push_back(a);
	B.push_back(b);
	X.push_back(x);
	if (n == 0)
		while (e > 0.5*pow(10, -k))
		{
			x0 = x;
			if (f(a)*f(x) > 0)
				a = x;
			else
				b = x;
			x = F();
			A.push_back(a);
			B.push_back(b);
			X.push_back(x);
			e = abs((x - x0) / x0);
		}
	if (k == 0)
		for (register int i = 1; i <= n; i++)
		{
			x0 = x;
			if (f(a)*f(x) > 0)
				a = x;
			else
				b = x;
			x = F();
			A.push_back(a);
			B.push_back(b);
			X.push_back(x);
		}
	std::string stra = "", strb = "", strx = "", line = "";
	array<std::string> lines(A.Size());
	array<std::string> aCol(A.Size()), bCol(A.Size()), xCol(A.Size());

	for (register int i = 0; i < A.Size(); i++)
	{
		stra = std::to_string(A[i]);
		strb = std::to_string(B[i]);
		strx = std::to_string(X[i]);
		for (register int j = stra.size() - 1; j >= 0; j--)
			if (A[i] == std::floor(A[i]))
			{
				stra = std::to_string(int(A[i]));
				break;
			}
			else if (stra[j] == '0' && stra.find('.') != -1 || stra[j] == '.')
				stra.pop_back();
			else
				break;
		for (register int j = strx.size() - 1; j >= 0; j--)
			if (X[i] == std::floor(X[i]))
			{
				strx = std::to_string(int(X[i]));
				break;
			}
			else if (strx[j] == '0' && strx.find('.') != -1 || strx[j] == '.')
				strx.pop_back();
			else
				break;
		for (register int j = strb.size() - 1; j >= 0; j--)
			if (B[i] == std::floor(B[i]))
			{
				strb = std::to_string(int(B[i]));
				break;
			}
			else if (strb[j] == '0'&&strb.find('.') != -1 || strb[j] == '.')
				strb.pop_back();
			else
				break;
		if (stra.size() > 11)
			stra = stra.substr(0, 10);
		if (strx.size() > 11)
			strx = strx.substr(0, 10);
		if (strb.size() > 11)
			strb = strb.substr(0, 10);
		aCol.push_back(stra);
		xCol.push_back(strx);
		bCol.push_back(strb);
	}
	line += "N";
	for (register int i = 0; i < int(log10(A.Size())) + 1; i++)
		line += " ";
	line += "|";
	for (register int i = 0; i < 11; i++)
		if (i == 5)
			line += "a";
		else
			line += " ";
	line += "|";
	for (register int i = 0; i < 11; i++)
		if (i == 5)
			line += "x";
		else
			line += " ";
	line += "|     b\n";
	lines.push_back(line);
	for (register int i = 1; i <= A.Size(); i++)
	{
		line = std::to_string(i);
		for (register int j = 0; j < int(log10(A.Size())) - int(log10(i)); j++)
			line += " ";
		line += " |";
		int a1 = (11 - aCol[i - 1].size()) / 2, a2 = ceil((11 - aCol[i - 1].size()) / 2.0);
		int x1 = (11 - xCol[i - 1].size()) / 2, x2 = ceil((11 - xCol[i - 1].size()) / 2.0);
		int b1 = (11 - bCol[i - 1].size()) / 2;
		for (register int j = 0; j < a1; j++)
			line += " ";
		line += aCol[i - 1];
		for (register int j = 0; j < a2; j++)
			line += " ";
		line += "|";
		for (register int j = 0; j < x1; j++)
			line += " ";
		line += xCol[i - 1];
		for (register int j = 0; j < x2; j++)
			line += " ";
		line += "|";
		for (register int j = 0; j < b1; j++)
			line += " ";
		line += bCol[i - 1];
		line += "\n";
		lines.push_back(line);
	}
	for (register int i = 0; i < lines.Size(); i++)
		std::cout << lines[i];
	for (register int i = 0; i < lines[lines.Size() - 1].size() - 1; i++)
		std::cout << "=";
	std::cout << "|" << std::endl << std::endl;
	return x;
}

double secant(double(*f)(double x), double x0, double x00, int n, int k)
{
	double x;
	double e = 1;
	int counter = 0;
	if (n == 0)
		while (e > 0.5*pow(10, -k))
		{
			x = x0 - f(x0)*(x0 - x00) / (f(x0) - f(x00));
			counter++;
			e = abs((x - x0) / x00);
			printf("x%d = %f\n", counter, x);
			x00 = x0; x0 = x;
		}
	if (k == 0)
		for (register int i = 1; i <= n; i++)
		{
			x = x0 - f(x0)*(x0 - x00) / (f(x0) - f(x00));
			e = abs((x - x0) / x00);
			printf("x%d = %f\n", counter, x);
			x00 = x0; x0 = x;
		}
	return x;
}

double newton(double(*F)(double x), double(*f)(double x), double x0, int n, int k)
{
	double x = x0;
	double e = 1;
	int counter = 0;
	if (n == 0)
		while (e > 0.5*pow(10, -k))
		{
			x = x0 - F(x0) / f(x0);
			counter++;
			e = abs(x - x0);
			printf("x%d = %f\n", counter, x);
			x0 = x;
		}
	if (k == 0)
		for (register int i = 1; i <= n; i++)
		{
			x = x0 - F(x0) / f(x0);
			printf("x%d = %f\n", i, x);
			x0 = x;
		}
	return x;
}

double newton(double(*f)(double), double x0 = 0)
{
	double e = 1, x;
	while (e > 5.0e-5)
	{
		x = x0 - f(x0) / derive(f, x0);
		e = std::fabs(x - x0);
		x0 = x;
	}
	return x;
}

///=========================================  NON-LINEAR SYSTEM OF EQUATIONS  ======================================================
void NewtonRaphson(double(*f)(double x, double y), double(*g)(double x, double y), double& x0, double& y0, int n, int k)
{
	double f0;
	double fx;
	double fy;
	double g0;
	double gx;
	double gy;
	double ex = 1, ey = 1, x, y;
	printf("0: (%g,%g)\n", x0, y0);
	if (n == 0)
		while (ex > 0.5*pow(10, -k) || ey > 0.5*pow(10, -k))
		{
			static int i = 0;
			i++;
			f0 = f(x0, y0);
			fx = partialx(f, x0, y0);
			fy = partialy(f, x0, y0);
			g0 = g(x0, y0);
			gx = partialx(g, x0, y0);
			gy = partialy(g, x0, y0);
			x = x0 + (fy * g0 - f0 * gy) / (fx * gy - fy * gx);
			y = y0 + (f0 * gx - fx * g0) / (fx * gy - fy * gx);
			printf("%d: (%g, %g)\n", i, x, y);
			ex = abs(x - x0);
			ey = abs(y - y0);
			x0 = x;  y0 = y;
		}
	if (k == 0)
		for (register int i = 1; i <= n; i++)
		{
			f0 = f(x0, y0);
			fx = partialx(f, x0, y0);
			fy = partialy(f, x0, y0);
			g0 = g(x0, y0);
			gx = partialx(g, x0, y0);
			gy = partialy(g, x0, y0);
			x = x0 + (fy * g0 - f0 * gy) / (fx * gy - fy * gx);
			y = y0 + (f0 * gx - fx * g0) / (fx * gy - fy * gx);
			printf("%d: (%g, %g)\n", i, x, y);
			x0 = x;  y0 = y;
		}
}

//*************************************************************** INTEGRATION **********************************************************************
double rectangular(double(*f)(double x), double a, double b, int n, double h)
{
	double I, sum = 0;
	if (n == 0)
		n = (b - a) / h;
	else
		h = (b - a) / n;
	for (register int i = 0; i < n; i++)
		sum += f(a + (i + 0.5)*h);
	I = h * sum;
	return I;
}

double trapezoidal(double(*f)(double x), double a, double b, int n, double h)
{
	double I, sum = 0;
	if (n == 0)
		n = (b - a) / h;
	else
		h = (b - a) / n;
	for (register int i = 1; i < n; i++)
		sum += f(a + i * h);
	I = (f(a) + f(b) + 2 * sum) * h / 2;
	return I;
}

double simpsons(double(*f)(double x), double a, double b, int n, double h)
{
	double I, Esum = 0, Osum = 0;
	if (n == 0)
		n = (b - a) / h;
	else
		h = (b - a) / n;
	for (register int i = 1; i < n / 2; i++)
		Esum += f(a + 2 * i*h);
	for (register int i = 1; i <= n / 2; i++)
		Osum += f(a + (2 * i - 1) * h);
	I = (f(a) + f(b) + 2 * Esum + 4 * Osum) * h / 3;
	return I;
}
/*  Applications of definite integration:
	 -average value of a function
	 -arc length
	 -area under curve
	 -surface area
	 -volume
 */
//*********************************************************** FUNCTION APPROXIMATION  ***********************************************************
pair<double, double> BFL(const real::vec& x, const real::vec& y);

real::vec fit(const real::vec& x, const real::vec& y, int d = -1)
{
	chk_eq(x.Size(), y.Size(), "The Number of x-points is not Equal to the Number of y-points");
	if (d == 0)
		return y.mean();
	if (d == -1)
		d = x.Size() - 1;
	chk_rng(d, x.Size(), "The Degree is Greater than the Number of Points");
	real::mat A(d + 1, 1), Y(d + 1, 1);
	Y.setelement(0, 0, y.sum());
	real::mat N(d + 1);
	N.setelement(0, 0, y.Size());
	N.setelement(d, d, x.Psum(2 * d));
	for (register int i = 1; i <= d; i++)
		Y.setelement(i, 0, (y.Wsum(x^i)));
	for (register int i = 1; i <= d; i++)
	{
		double val1 = x.Psum(i), val2 = x.Psum(2 * d - i);
		for (register int j = 0; j <= i; j++)
		{
			N.setelement(j, i - j, val1);
			N.setelement(d - i + j, d - j, val2);
		}
	}
	A = (!N)*Y;
	return A.unroll();
}

polynomial Lagrange();

polynomial* CubicSplines();

//*********************************************************** Eigen ***********************************************************
real::vec power(const real::mat& A, real::vec& x)
{
	real::mat res = A * A * A;
	double x0 = 0, x1, e = 1;
	while (e > 5.0e-7)
	{
		x1 = x0;
		x = res * x;
		x0 = x.abs_max();
		if (x0 < 0)
			x0 = -std::pow(-x0, 1.0 / 3);
		else
			x0 = std::pow(x0, 1.0 / 3);
		x._norm();
		e = x0 - x1;
	}
	return x0*x;
};

real::vec inversePower(const real::mat& A, real::vec& x);

real::vec shiftedInversePower(const real::mat& A, real::vec& x);
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

// Least Square Error:
double LSE(real::vec& x, real::vec& y, real::vec& a)
{
	double sum = 0, h;
	for (register int i = 0; i < x.Size(); i++)
	{
		h = 0;
		for (register int j = 0; j < a.Size(); j++)
			h += a[j] * pow(x[i], j);
		sum += pow(y[i] - h, 2);
	}
	return sum;
}
double Sum(double start, double end, double(*F)(double))
{
	double sum = 0;
	for (register int i = start; i <= end; i++)
		sum += F(i);
	return sum;
}
//===================  FILE HANDELING  =======================
void fprint(const real::vec& v, std::string file_Name)
{
	std::ofstream file(file_Name + ".txt");
	int size = v.Size();
	file << "real::vector\n" << size << "\n";
	file << "( " << v[0];
	for (register int i = 1; i < size; i++)
		file << ", " << v[i];
	file << ")";
	file.close();
}
void fprint(std::ofstream& file, const real::vec& v, std::string vecName)
{
	int size = v.Size();
	file << "real::vector\n" << size << "\n";
	file << vecName << " : ( " << v[0];
	for (register int i = 1; i < size; i++)
		file << ", " << v[i];
	file << ")";
	file.close();
}
//===================  <<---------->>  =======================


#undef fileName