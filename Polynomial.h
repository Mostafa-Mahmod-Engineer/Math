#pragma once
#include"Mat.h"
#define fileName "Polynomial.h"

struct term
{
	int power; double coef;
	term() { power = 0; coef = 0; }
	term(double c, int p) : coef(c), power(p) {  }
	term operator*(term& t) { return term(coef*t.coef, power + t.power); }
	term operator/(term& t) { zero_check(t.coef); return term(coef / t.coef, power - t.power); }

	term operator*(term& t) const { return term(coef*t.coef, power + t.power); }
	term operator/(term& t) const { zero_check(t.coef); return term(coef / t.coef, power - t.power); }

	term operator*(double val) { return term(val*coef, power); }
	term operator/(double val) { zero_check(val); return term(coef / val, power); }

	term operator*(double val) const { return term(val*coef, power); }
	term operator/(double val) const { zero_check(val); return term(coef / val, power); }
	void print()
	{
		printf(" %gX^%d%\n", coef, power);
	}
	~term() {  }
};

class polynomial
{
	double* coef = nullptr; int degree;
public:
	polynomial() {}
	explicit polynomial(int d) { degree = d; if (d != -1) coef = new double[d + 1]{}; }
	explicit polynomial(double a0)
	{
		if (a0 == 0)
			degree = -1;
		else
		{
			degree = 0;
			coef = new double{ a0 };
		}
	}
	template<class T>
	void assign(const T& t) { coef[0] = t; }
	template<class H, class ... T>
	void assign(const H& h, const T& ... t)
	{
		coef[sizeof...(t)] = h;
		assign(t...);
	}
	template<class H, class ... T>
	polynomial(const H& h, const T& ... t)
	{
		degree = sizeof...(T);
		coef = new double[degree + 1]{};
		assign(h, t...);
	}
	polynomial(const polynomial& other)
	{
		degree = other.degree;
		coef = new double[degree + 1];
		for (register int i = 0; i <= degree; i++)
			coef[i] = other.coef[i];
	}
	polynomial(polynomial&& other)
	{
		degree = other.degree;
		coef = other.coef;
		other.coef = nullptr;
	}
	polynomial(polynomial* other)
	{
		degree = other->degree;
		coef = new double[degree + 1];
		for (register int i = 0; i <= degree; i++)
			coef[i] = other->coef[i];
	}
	polynomial(const polynomial* other)
	{
		degree = other->degree;
		coef = new double[degree + 1];
		for (register int i = 0; i <= degree; i++)
			coef[i] = other->coef[i];
	}
	polynomial(term& t)
	{
		degree = t.power;
		coef = new double[degree + 1]{};
		coef[degree] = t.coef;
	}
	polynomial(const real::vec& v)
	{
		degree = v.Size() - 1;
		int k;
		for (k = 0; k <= degree; k++)
			if (v[k] != 0)
				break;
		degree -= k;
		coef = new double[degree];
		for (register int i = 0; i <= degree; i++)
			coef[degree - i] = v[i + k];
	}
	polynomial& operator =(const real::vec& v)
	{
		if (coef != nullptr)
			delete[] coef;
		degree = v.Size() - 1;
		int k;
		for (k = 0; k <= degree; k++)
			if (v[k] != 0)
				break;
		degree -= k;
		coef = new double[degree];
		for (register int i = 0; i <= degree; i++)
			coef[degree - i] = v[i + k];
		return *this;
	}
	//=============================== << G E T T E R S >> ====================================
	int Degree() { return degree; }
	double Coef(int i) { chk_rng(i, degree + 1); return coef[i]; }
	double lead() { return coef[degree]; }
	double free() { return coef[0]; }
	bool is_zero()
	{
		if (degree == -1)
			return true;
		if (degree == 0)
			return coef[0] == 0;
		for (register int i = 0; i <= degree; i++)
			if (coef[i] != 0)
				return false;
		return true;
	}

	int Degree() const { return degree; }
	double Coef(int i) const { chk_rng(i, degree + 1); return coef[i]; }
	double lead() const { return coef[degree]; }
	double free() const { return coef[0]; }
	bool is_zero() const
	{
		if (degree == -1)
			return true;
		if (degree == 0)
			return coef[0] == 0;
		for (register int i = 0; i <= degree; i++)
			if (coef[i] != 0)
				return false;
		return true;
	}
	//=============================== << G E T T E R S \E.N.D>> ===============================

	//=============================== << F U N C T I O N S >> =================================
	double sub(double x)
	{
		if (degree == -1)
			return 0;
		double res = coef[degree];
		for (register int i = degree; i > 0; i--)
		{
			res *= x;
			res += coef[i - 1];
		}
		return res;
	}
	polynomial reduce()
	{
		if (degree == -1 || coef[degree] != 0)
			return *this;
		int d = degree;
		for (register int i = degree; i >= 0; i--)
			if (coef[i] == 0)
				d--;
			else
				break;
		polynomial res(d);
		for (register int i = 0; i <= d; i++)
			res.coef[i] = coef[i];
		return res;
	}
	polynomial differentiate()
	{
		if (is_zero() || degree == 0) return polynomial(-1);
		polynomial res(degree - 1);
		for (register int i = 0; i < degree; i++)
			res.coef[i] = coef[i + 1] * (i + 1);
		return res;
	}
	polynomial integrate(double c = 0)
	{
		if (is_zero()) return polynomial(c);
		polynomial res(degree + 1);
		res.coef[0] = c;
		for (register int i = 1; i < degree + 2; i++)
			res.coef[i] = coef[i - 1] / i;
		return res;
	}
	polynomial shift(double val)
	{
		polynomial res(degree);
		for (register int i = 1; i <= degree; i++)
			res += (polynomial(1, -val) ^ i)*coef[i];
		res.coef[0] += coef[0];
		return res;
	}
	polynomial monic() { return *this / lead(); }

	void _reduce() { degree = reduce().degree; }
	void _differentiate() { coef = differentiate().coef; degree--; }
	void _integrate(double c = 0) { coef = integrate(c).coef; degree++; }
	void _shift(double val) { coef = shift(val).coef; }
	void _monic() { for (register int i = 0; i <= degree; i++) coef[i] /= coef[degree]; }

	double sub(double x) const
	{
		double res = coef[degree];
		for (register int i = degree; i > 0; i--)
		{
			res *= x;
			res += coef[i - 1];
		}
		return res;
	}
	polynomial reduce() const
	{
		if (degree == -1 || coef[degree] != 0)
			return *this;
		int d = degree;
		for (register int i = degree; i >= 0; i--)
			if (coef[i] == 0)
				d--;
			else
				break;
		polynomial res(d);
		for (register int i = 0; i <= d; i++)
			res.coef[i] = coef[i];
		return res;
	}
	polynomial differentiate() const
	{
		if (is_zero() || degree == 0) return polynomial(-1);
		polynomial res(degree - 1);
		for (register int i = 0; i < degree; i++)
			res.coef[i] = coef[i + 1] * (i + 1);
		return res;
	}
	polynomial integrate(double c = 0) const
	{
		if (is_zero()) return polynomial(c);
		polynomial res(degree + 1);
		res.coef[0] = c;
		for (register int i = 1; i < degree + 2; i++)
			res.coef[i] = coef[i - 1] / i;
		return res;
	}
	polynomial shift(double val) const
	{
		polynomial res(degree);
		for (register int i = 1; i <= degree; i++)
			res += (polynomial(1, -val) ^ i)*coef[i];
		res.coef[0] += coef[0];
		return res;
	}
	polynomial monic() const { return *this / lead(); }
	//================================== << S T A T I C >> ====================================
	static polynomial zeroes(const real::vec& v)
	{
		polynomial res(1, -v[0]); int size = v.Size();
		for (register int i = 1; i < size; i++)
			res *= polynomial(1, -v[i]);
		return res;
	}
	static polynomial Czeroes(const real::vec& v)
	{
		polynomial res(2, -v[0]); int size = v.Size();
		for (register int i = 1; i < size; i++)
			res *= polynomial(2, -v[i]);
		return res;
	}
	static polynomial R(int n)
	{
		polynomial res(n);
		for (register int i = 0; i < n; i++)
			res.coef[i] = std::rand() % 10;
		res.coef[n] = (std::rand() % 9) + 1;
		return res;
	}
	static polynomial R(int n, double upperlimit)
	{
		polynomial res(n);
		for (register int i = 0; i < n; i++)
			res.coef[i] = std::fmodf(std::rand(), upperlimit);
		res.coef[n] = std::fmodl(double(std::rand()), upperlimit - 1) + 1;
		return res;
	}
	static polynomial R(int n, double lowerlimit, double upperlimit)
	{
		polynomial res(n);
		for (register int i = 0; i <= n; i++)
			res.coef[i] = std::fmodl(double(std::rand()), upperlimit - lowerlimit) + lowerlimit;
		if (res.coef[n] == 0)
			res.coef[n]++;
		return res;
	}
	//=============================== << F U N C T I O N S \E.N.D>> ===========================

	//=============================== << O P E R A T O R S >> =================================
	polynomial operator-()
	{
		polynomial res(degree);
		for (register int i = 0; i <= degree; i++)
			res.coef[i] = -coef[i];
		return res;
	}
	polynomial operator-() const
	{
		polynomial res(degree);
		for (register int i = 0; i <= degree; i++)
			res.coef[i] = -coef[i];
		return res;
	}
	//===================================    POLY X POLY    ===================================
	polynomial operator+(const polynomial& other)
	{
		if (degree > other.degree)
		{
			polynomial res(degree);
			for (register int i = 0; i <= other.degree; i++)
				res.coef[i] = coef[i] + other.coef[i];
			for (register int i = other.degree + 1; i <= degree; i++)
				res.coef[i] = coef[i];
			return res;
		}
		else if (degree < other.degree)
		{
			polynomial res(other.degree);
			for (register int i = 0; i <= degree; i++)
				res.coef[i] = coef[i] + other.coef[i];
			for (register int i = degree + 1; i <= other.degree; i++)
				res.coef[i] = other.coef[i];
			return res;
		}
		else
		{
			int d = degree;
			for (register int i = degree; i >= 0; i--)
				if (coef[i] == -other.coef[i])
					d--;
				else
					break;
			polynomial res(d);
			for (register int i = 0; i <= d; i++)
				res.coef[i] = coef[i] + other.coef[i];
			return res;
		}
	}
	polynomial operator-(const polynomial& other)
	{
		if (degree > other.degree)
		{
			polynomial res(degree);
			for (register int i = 0; i <= other.degree; i++)
				res.coef[i] = coef[i] - other.coef[i];
			for (register int i = other.degree + 1; i <= degree; i++)
				res.coef[i] = coef[i];
			return res;
		}
		else if (degree < other.degree)
		{
			polynomial res(other.degree);
			for (register int i = 0; i <= degree; i++)
				res.coef[i] = coef[i] - other.coef[i];
			for (register int i = degree + 1; i <= other.degree; i++)
				res.coef[i] = -other.coef[i];
			return res;
		}
		else
		{
			int d = degree;
			for (register int i = degree; i >= 0; i--)
				if (coef[i] == other.coef[i])
					d--;
				else
					break;
			polynomial res(d);
			for (register int i = 0; i <= d; i++)
				res.coef[i] = coef[i] - other.coef[i];
			return res;
		}
	}
	polynomial operator*(const polynomial& other)
	{
		polynomial res(degree + other.degree);
		for (register int i = 0; i < res.degree + 1; i++)
		{
			if (degree >= other.degree)
			{
				int m = degree, n = other.degree;
				if (i <= n)
					for (register int j = 0; j <= i; j++)
						res.coef[i] += coef[j] * other.coef[i - j];
				else if (i <= m)
					for (register int j = i - n; j <= i; j++)
						res.coef[i] += coef[j] * other.coef[i - j];
				else
					for (register int j = i - n; j <= m; j++)
						res.coef[i] += coef[j] * other.coef[i - j];
			}
			else
			{
				int m = other.degree, n = degree;
				if (i <= n)
					for (register int j = 0; j <= i; j++)
						res.coef[i] += other.coef[j] * coef[i - j];
				else if (i <= m)
					for (register int j = i - n; j <= i; j++)
						res.coef[i] += other.coef[j] * coef[i - j];
				else
					for (register int j = i - n; j <= m; j++)
						res.coef[i] += other.coef[j] * coef[i - j];
			}
		}
		return res;
	}
	polynomial operator/(const polynomial& other)
	{
		if (other.is_zero())
			zero_check(0);
		if (other.degree > degree)
			return polynomial(-1);
		else if (other.degree == degree)
			return polynomial(coef[degree] / other.coef[degree]);
		else if (other.degree == 1)
		{
			double pivot = -other.coef[0] / other.coef[1];
			real::mat coefTable(3, degree + 1);
			for (register int i = 0; i <= degree; i++)
				coefTable.setelement(0, i, coef[degree - i] / other.coef[1]);
			coefTable.setelement(2, 0, coefTable.get(0, 0));
			for (register int i = 1; i <= degree; i++)
			{
				coefTable.setelement(1, i, pivot*coefTable.get(2, i - 1));
				coefTable.setelement(2, i, coefTable.get(0, i) + coefTable.get(1, i));
			}
			polynomial res(degree - 1);
			for (register int i = 0; i < degree; i++)
				res.coef[degree - 1 - i] = coefTable.get(2, i);
			return res;
		}
		else
		{
			double *pivots = new double[other.degree];
			for (register int i = 0; i < other.degree; i++)
				pivots[i] = -other.coef[i] / other.coef[other.degree];
			//______________________________________________________________________
			real::mat coefTable(other.degree + 2, degree + 1);
			// first row for the top P - last row for the result 
			// middle rows for the pivots of the bottom
			for (register int i = 0; i <= degree; i++)
				coefTable.setelement(0, i, coef[degree - i] / other.coef[other.degree]);
			//______________________________________________________________________
			coefTable.setelement(other.degree + 1, 0, coefTable.get(0, 0)); double sum = 0;
			for (register int i = 1; i <= degree - other.degree + 1; i++)
			{
				for (register int j = i + other.degree - 1, k = 1; j >= i; j--, k++)
					coefTable.setelement(k, j, pivots[k - 1] * coefTable.get(other.degree + 1, i - 1));
				sum = 0;
				for (register int j = 0; j <= other.degree; j++)
					sum += coefTable.get(j, i);
				coefTable.setelement(other.degree + 1, i, sum);
			}
			polynomial res(degree - other.degree);
			for (register int i = 0; i <= res.degree; i++)
				res.coef[res.degree - i] = coefTable.get(other.degree + 1, i);
			return res;
		}
	}
	polynomial operator%(const polynomial& other)
	{
		if (other.is_zero())
			zero_check(0);
		if (other.degree > degree)
			return *this;
		else if (other.degree == 1)
		{
			double pivot = -other.coef[0] / other.coef[1];
			real::mat coefTable(3, degree + 1);
			// first row for the top P - last row for the result 
			// middle rows for the pivots of the bottom
			for (register int i = 0; i <= degree; i++)
				coefTable.setelement(0, i, coef[degree - i] / other.coef[1]);
			coefTable.setelement(2, 0, coefTable.get(0, 0));
			for (register int i = 1; i <= degree; i++)
			{
				coefTable.setelement(1, i, pivot*coefTable.get(2, i - 1));
				coefTable.setelement(2, i, coefTable.get(0, i) + coefTable.get(1, i));
			}
			polynomial res(0);
			coefTable.get(2, degree) > 0.00001 ? res.coef[0] = coefTable.get(2, degree) : res.coef[0] = 0;
			return (res * other.lead()).reduce();
		}
		else
		{
			double *pivots = new double[other.degree];
			for (register int i = 0; i < other.degree; i++)
				pivots[i] = -other.coef[i] / other.coef[other.degree];

			real::mat coefTable(other.degree + 2, degree + 1);
			// first row for the top P - last row for the result 

			// middle rows for the pivots of the bottom
			for (register int i = 0; i <= degree; i++)
				coefTable.setelement(0, i, coef[degree - i] / other.coef[other.degree]);

			coefTable.setelement(other.degree + 1, 0, coefTable.get(0, 0)); double sum = 0;
			for (register int i = 1; i <= degree - other.degree + 1; i++)
			{
				for (register int j = i + other.degree - 1, k = 1; j >= i; j--, k++)
					coefTable.setelement(k, j, pivots[k - 1] * coefTable.get(other.degree + 1, i - 1));
				sum = 0;
				for (register int j = 0; j <= other.degree; j++)
					sum += coefTable.get(j, i);
				coefTable.setelement(other.degree + 1, i, sum);
			}
			for (register int i = degree - other.degree + 2; i <= degree; i++)
			{
				sum = 0;
				for (register int j = 0; j <= other.degree; j++)
					sum += coefTable.get(j, i);
				coefTable.setelement(other.degree + 1, i, sum);
			}
			//system("pause");
			int d = other.degree - 1;
			for (register int i = degree - other.degree + 1; i <= degree; i++)
				if (coefTable.get(other.degree + 1, i) == 0)
					d--;
				else
					break;
			if (d == -1 || d == 0 && coefTable.get(other.degree + 1, degree) == 0)
				return polynomial(-1);
			else if (d == 0)
				return polynomial(coefTable.get(other.degree + 1, degree));
			polynomial res(d);

			for (register int i = 0; i <= d; i++)
				coefTable.get(other.degree + 1, degree - i) > 0.00001 ? res.coef[i] = coefTable.get(other.degree + 1, degree - i) : res.coef[i] = 0;
			return (res * other.lead()).reduce();
		}
	}
	bool operator==(const polynomial& other)
	{
		chk_eq(degree, other.degree);
		for (register int i = 0; i <= degree; i++)
			if (coef[i] != other.coef[i])
				return false;
		return true;
	}
	bool operator!=(const polynomial& other) { return !(*this == other); }

	polynomial operator+(const polynomial& other) const
	{
		if (degree > other.degree)
		{
			polynomial res(degree);
			for (register int i = 0; i <= other.degree; i++)
				res.coef[i] = coef[i] + other.coef[i];
			for (register int i = other.degree + 1; i <= degree; i++)
				res.coef[i] = coef[i];
			return res;
		}
		else if (degree < other.degree)
		{
			polynomial res(other.degree);
			for (register int i = 0; i <= degree; i++)
				res.coef[i] = coef[i] + other.coef[i];
			for (register int i = degree + 1; i <= other.degree; i++)
				res.coef[i] = other.coef[i];
			return res;
		}
		else
		{
			int d = degree;
			for (register int i = degree; i >= 0; i--)
				if (coef[i] == -other.coef[i])
					d--;
				else
					break;
			polynomial res(d);
			for (register int i = 0; i <= d; i++)
				res.coef[i] = coef[i] + other.coef[i];
			return res;
		}
	}
	polynomial operator-(const polynomial& other) const
	{
		if (degree > other.degree)
		{
			polynomial res(degree);
			for (register int i = 0; i <= other.degree; i++)
				res.coef[i] = coef[i] - other.coef[i];
			for (register int i = other.degree + 1; i <= degree; i++)
				res.coef[i] = coef[i];
			return res;
		}
		else if (degree < other.degree)
		{
			polynomial res(other.degree);
			for (register int i = 0; i <= degree; i++)
				res.coef[i] = coef[i] - other.coef[i];
			for (register int i = degree + 1; i <= other.degree; i++)
				res.coef[i] = -other.coef[i];
			return res;
		}
		else
		{
			int d = degree;
			for (register int i = degree; i >= 0; i--)
				if (coef[i] == other.coef[i])
					d--;
				else
					break;
			polynomial res(d);
			for (register int i = 0; i <= d; i++)
				res.coef[i] = coef[i] - other.coef[i];
			return res;
		}
	}
	polynomial operator*(const polynomial& other) const
	{
		polynomial res(degree + other.degree);
		for (register int i = 0; i < res.degree + 1; i++)
		{
			if (degree >= other.degree)
			{
				int m = degree, n = other.degree;
				if (i <= n)
					for (register int j = 0; j <= i; j++)
						res.coef[i] += coef[j] * other.coef[i - j];
				else if (i <= m)
					for (register int j = i - n; j <= i; j++)
						res.coef[i] += coef[j] * other.coef[i - j];
				else
					for (register int j = i - n; j <= m; j++)
						res.coef[i] += coef[j] * other.coef[i - j];
			}
			else
			{
				int m = other.degree, n = degree;
				if (i <= n)
					for (register int j = 0; j <= i; j++)
						res.coef[i] += other.coef[j] * coef[i - j];
				else if (i <= m)
					for (register int j = i - n; j <= i; j++)
						res.coef[i] += other.coef[j] * coef[i - j];
				else
					for (register int j = i - n; j <= m; j++)
						res.coef[i] += other.coef[j] * coef[i - j];
			}
		}
		return res;
	}
	polynomial operator/(const polynomial& other) const
	{
		if (other.is_zero())
			zero_check(0);
		if (other.degree > degree)
			return polynomial(-1);
		else if (other.degree == degree)
			return polynomial(coef[degree] / other.coef[degree]);
		else if (other.degree == 1)
		{
			double pivot = -other.coef[0] / other.coef[1];
			real::mat coefTable(3, degree + 1);
			for (register int i = 0; i <= degree; i++)
				coefTable.setelement(0, i, coef[degree - i] / other.coef[1]);
			coefTable.setelement(2, 0, coefTable.get(0, 0));
			for (register int i = 1; i <= degree; i++)
			{
				coefTable.setelement(1, i, pivot*coefTable.get(2, i - 1));
				coefTable.setelement(2, i, coefTable.get(0, i) + coefTable.get(1, i));
			}
			polynomial res(degree - 1);
			for (register int i = 0; i < degree; i++)
				res.coef[degree - 1 - i] = coefTable.get(2, i);
			return res;
		}
		else
		{
			double *pivots = new double[other.degree];
			for (register int i = 0; i < other.degree; i++)
				pivots[i] = -other.coef[i] / other.coef[other.degree];
			//______________________________________________________________________
			real::mat coefTable(other.degree + 2, degree + 1);
			// first row for the top P - last row for the result 
			// middle rows for the pivots of the bottom
			for (register int i = 0; i <= degree; i++)
				coefTable.setelement(0, i, coef[degree - i] / other.coef[other.degree]);
			//______________________________________________________________________
			coefTable.setelement(other.degree + 1, 0, coefTable.get(0, 0)); double sum = 0;
			for (register int i = 1; i <= degree - other.degree + 1; i++)
			{
				for (register int j = i + other.degree - 1, k = 1; j >= i; j--, k++)
					coefTable.setelement(k, j, pivots[k - 1] * coefTable.get(other.degree + 1, i - 1));
				sum = 0;
				for (register int j = 0; j <= other.degree; j++)
					sum += coefTable.get(j, i);
				coefTable.setelement(other.degree + 1, i, sum);
			}
			polynomial res(degree - other.degree);
			for (register int i = 0; i <= res.degree; i++)
				res.coef[res.degree - i] = coefTable.get(other.degree + 1, i);
			return res;
		}
	}
	polynomial operator%(const polynomial& other) const
	{
		if (other.is_zero())
			zero_check(0);
		if (other.degree > degree)
			return *this;
		else if (other.degree == 1)
		{
			double pivot = -other.coef[0] / other.coef[1];
			real::mat coefTable(other.degree + 2, degree + 1);
			// first row for the top P - last row for the result 
			// middle rows for the pivots of the bottom
			for (register int i = 0; i <= degree; i++)
				coefTable.setelement(0, i, coef[degree - i] / other.coef[1]);
			coefTable.setelement(2, 0, coefTable.get(0, 0));
			for (register int i = 1; i <= degree; i++)
			{
				coefTable.setelement(1, i, pivot*coefTable.get(2, i - 1));
				coefTable.setelement(2, i, coefTable.get(0, i) + coefTable.get(1, i));
			}
			polynomial res(0);
			coefTable.get(2, degree) > 0.00001 ? res.coef[0] = coefTable.get(2, degree) : res.coef[0] = 0;
			return (res * other.lead()).reduce();
		}
		else
		{
			double *pivots = new double[other.degree];
			for (register int i = 0; i < other.degree; i++)
				pivots[i] = -other.coef[i] / other.coef[other.degree];

			real::mat coefTable(other.degree + 2, degree + 1);
			// first row for the top P - last row for the result 

			// middle rows for the pivots of the bottom
			for (register int i = 0; i <= degree; i++)
				coefTable.setelement(0, i, coef[degree - i] / other.coef[other.degree]);

			coefTable.setelement(other.degree + 1, 0, coefTable.get(0, 0)); double sum = 0;
			for (register int i = 1; i <= degree - other.degree + 1; i++)
			{
				for (register int j = i + other.degree - 1, k = 1; j >= i; j--, k++)
					coefTable.setelement(k, j, pivots[k - 1] * coefTable.get(other.degree + 1, i - 1));
				sum = 0;
				for (register int j = 0; j <= other.degree; j++)
					sum += coefTable.get(j, i);
				coefTable.setelement(other.degree + 1, i, sum);
			}
			for (register int i = degree - other.degree + 2; i <= degree; i++)
			{
				sum = 0;
				for (register int j = 0; j <= other.degree; j++)
					sum += coefTable.get(j, i);
				coefTable.setelement(other.degree + 1, i, sum);
			}
			double d = other.degree - 1;
			for (register int i = degree - other.degree + 1; i <= degree; i++)
				if (coefTable.get(other.degree + 1, i) == 0)
					d--;
				else
					break;
			if (d == -1 || d == 0 && coefTable.get(other.degree + 1, degree) == 0)
				return polynomial(-1);
			else if (d == 0)
				return polynomial(coefTable.get(other.degree + 1, degree));
			polynomial res(d);

			for (register int i = 0; i <= d; i++)
				coefTable.get(other.degree + 1, degree - i) > 0.00001 ? res.coef[i] = coefTable.get(other.degree + 1, degree - i) : res.coef[i] = 0;
			return (res * other.lead()).reduce();
		}
	}
	bool operator==(const polynomial& other) const { chk_eq(degree, other.degree); for (register int i = 0; i <= degree; i++) if (coef[i] != other.coef[i])return false; return true; }
	bool operator!=(const polynomial& other) const { return !(*this == other); }

	polynomial& operator =(polynomial&& other)
	{
		degree = other.degree;
		coef = other.coef;
		other.coef = nullptr;
		return *this;
	}
	polynomial& operator =(const polynomial& other)
	{
		this->~polynomial();
		degree = other.degree;
		coef = new double[degree + 1];
		for (register int i = 0; i <= degree; i++)
			coef[i] = other.coef[i];
		return *this;
	}
	polynomial& operator+=(const polynomial& other) { *this = *this + other; return *this; }
	polynomial& operator-=(const polynomial& other) { *this = *this - other; return *this; }
	polynomial& operator*=(const polynomial& other) { *this = *this * other; return *this; }
	polynomial& operator/=(const polynomial& other) { *this = *this / other; return *this; }
	polynomial& operator%=(const polynomial& other) { *this = *this % other; return *this; }
	//===================================    POLY X TERM    ==================================
	polynomial operator+(term& other)
	{
		polynomial res(std::fmax(other.power, degree));
		for (register int i = 0; i <= degree; i++)
			res.coef[i] = coef[i];
		res.coef[other.power] += other.coef;
	}
	polynomial operator-(term& other)
	{
		polynomial res(std::fmax(other.power, degree));
		for (register int i = 0; i <= degree; i++)
			res.coef[i] = coef[i];
		res.coef[other.power] -= other.coef;
	}
	polynomial operator*(term& other)
	{
		polynomial res(degree + other.power);
		for (register int i = other.power; i <= res.degree; i++)
			res.coef[i] = coef[i - other.power] * other.coef;
		return res;
	}
	polynomial operator/(term& other)
	{
		if (other.power > degree)
			return polynomial(-1);
		polynomial res(degree - other.power);
		for (register int i = 0; i < res.degree + 1; i++)
			res.coef[i] = coef[i + other.power] / other.coef;
		return res;
	}
	polynomial operator%(term& other)
	{
		if (other.power > degree)
			return *this;
		polynomial res(other.power - 1);
		for (register int i = 0; i < other.power; i++)
			res.coef[i] = coef[i];
		return res;
	}
	bool operator==(term& other)
	{
		if (!(other.power == degree && other.coef == coef[degree]))
			return false;
		else
			for (register int i = 0; i < degree; i++)
				if (coef[i] != 0)
					return false;
		return true;
	}
	bool operator!=(term& other) { return !(*this == other); }

	polynomial operator+(term& other) const
	{
		polynomial res(std::fmax(other.power, degree));
		for (register int i = 0; i <= degree; i++)
			res.coef[i] = coef[i];
		res.coef[other.power] += other.coef;
	}
	polynomial operator-(term& other) const
	{
		polynomial res(std::fmax(other.power, degree));
		for (register int i = 0; i <= degree; i++)
			res.coef[i] = coef[i];
		res.coef[other.power] -= other.coef;
	}
	polynomial operator*(term& other) const
	{
		polynomial res(degree + other.power);
		for (register int i = other.power; i <= res.degree; i++)
			res.coef[i] = coef[i - other.power] * other.coef;
		return res;
	}
	polynomial operator/(term& other) const
	{
		if (other.power > degree)
			return polynomial(-1);
		polynomial res(degree - other.power);
		for (register int i = 0; i < res.degree + 1; i++)
			res.coef[i] = coef[i + other.power] / other.coef;
		return res;
	}
	polynomial operator%(term& other) const
	{
		if (other.power > degree)
			return *this;
		polynomial res(other.power - 1);
		for (register int i = 0; i < other.power; i++)
			res.coef[i] = coef[i];
		return res;
	}
	bool operator==(term& other) const
	{
		if (!(other.power == degree && other.coef == coef[degree]))
			return false;
		else
			for (register int i = 0; i < degree; i++)
				if (coef[i] != 0)
					return false;
		return true;
	}
	bool operator!=(term& other) const { return !(*this == other); }

	polynomial& operator+=(term& other) { *this = *this + other; return *this; }
	polynomial& operator-=(term& other) { *this = *this - other; return *this; }
	polynomial& operator*=(term& other) { *this = *this * other; return *this; }
	polynomial& operator/=(term& other) { *this = *this / other; return *this; }
	polynomial& operator%=(term& other) { *this = *this % other; return *this; }
	//===================================    POLY X NUM    ===================================
	polynomial operator+(double val) { polynomial res(this);  res.coef[0] += val; return res; }
	polynomial operator-(double val) { polynomial res(this); res.coef[0] -= val; return res; }
	polynomial operator*(double val)
	{
		polynomial res(degree);
		for (register int i = 0; i <= degree; i++)
			res.coef[i] = coef[i] * val;
		return res;
	}
	polynomial operator/(double val)
	{
		zero_check(val);
		polynomial res(degree);
		for (register int i = 0; i <= degree; i++)
			res.coef[i] = coef[i] / val;
		return res;
	}
	polynomial operator^(int p)
	{
		polynomial res(this);
		for (register int i = 1; i < p; i++)
			res *= *this;
		return res;
	}
	bool operator==(double val) { if (degree != 0)return false; return coef[0] == val; }
	bool operator!=(double val) { if (degree != 0)return true; return coef[0] != val; }

	polynomial operator+(double val) const { polynomial res(this); res.coef[0] += val; return res; }
	polynomial operator-(double val) const { polynomial res(this); res.coef[0] -= val; return res; }
	polynomial operator*(double val) const
	{
		polynomial res(degree);
		for (register int i = 0; i <= degree; i++)
			res.coef[i] = coef[i] * val;
		return res;
	}
	polynomial operator/(double val) const
	{
		zero_check(val);
		polynomial res(degree);
		for (register int i = 0; i <= degree; i++)
			res.coef[i] = coef[i] / val;
		return res;
	}
	polynomial operator^(int p) const
	{
		polynomial res(this);
		for (register int i = 1; i < p; i++)
			res *= *this;
		return res;
	}
	bool operator==(double val) const { if (degree != 0)return false; return coef[0] == val; }
	bool operator!=(double val) const { if (degree != 0)return true; return coef[0] != val; }

	polynomial& operator+=(double val) { *this = *this + val; return *this; }
	polynomial& operator-=(double val) { *this = *this - val; return *this; }
	polynomial& operator*=(double val) { *this = *this * val; return *this; }
	polynomial& operator/=(double val) { *this = *this / val; return *this; }
	polynomial& operator^=(int p) { *this = *this ^ p; return *this; }
	//=============================== << O P E R A T O R S \E.N.D>> ===================================
	void set(int i, double val) { chk_rng(i, degree + 1); coef[i] = val; }
	void print()
	{
		if (degree == 0 && coef[0] == 0 || degree == -1)
			printf(" 0");
		else
		{
			if (degree == 1)
				coef[1] == 1 ? printf(" X ") : coef[1] == -1 ? printf(" -X ") : printf(" %gX ", coef[1]);
			else if (degree > 1)
			{
				coef[degree] == 1 ? printf(" X^%d ", degree) : coef[degree] == -1 ? printf(" -X^%d ", degree) : printf(" %gX^%d ", coef[degree], degree);
				for (register int i = degree - 1; i >= 2; i--)
					if (coef[i] == 0)
						continue;
					else
						coef[i] == 1 ? printf("+X^%d ", i) : coef[i] == -1 ? printf("-X^%d ", i) : printf("%+gX^%d ", coef[i], i);
			}
			if (degree > 1 && coef[1] != 0)
				coef[1] == 1 ? printf("+X ") : coef[1] == -1 ? printf("-X ") : printf("%+gX ", coef[1]);
			if (coef[0] != 0) printf("%+g ", coef[0]);
		}
		printf("\n");
	}
	//friend std::ostream& operator<<(std::ostream& C, const polynomial& poly);
	~polynomial() { if (coef != nullptr) delete[] coef; }
};
std::ostream& operator<<(std::ostream& C, polynomial& p) { p.print(); return C; }


#undef fileName