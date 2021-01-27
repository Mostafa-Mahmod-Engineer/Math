#pragma once
#include"Data Structures.h"
#define fileName "Complex.h"

struct Com
{
	double a, b;
	Com() : a(0), b(0) {}
	Com(double a) : a(a), b(0) {}
	Com(double a, double b) : a(a), b(b) {}
	Com(std::string str)
	{
		if (str == "i")
		{
			a = 0;
			b = 1;
		}
		if (str == "-i")
		{
			a = 0;
			b = -1;
		}
		if (str.find('i') == -1)
		{
			a = std::atof(str.c_str());
			b = 0;
		}
		else
		{
			int size = str.size();
			bool flag = false; int signIndex;
			for (register int i = 1; i < size; i++)
				if (str[i] == '+' || str[i] == '-')
				{
					flag = true;
					signIndex = i;
					break;
				}
			if (flag)
			{
				a = std::atof(str.c_str());
				str = str.substr(signIndex);
				if (str == "+i")
					b = 1;
				else if (str == "-i")
					b = -1;
				else
					b = std::atof(str.c_str());
			}
			else
			{
				a = 0;
				b = std::atof(str.c_str());
			}
		}
	}
	double mag() const { return std::sqrt(a*a + b * b); }
	double theta() const
	{
		if (a == 0 && b == 0)
			return 0;
		if (a <= 0 && b >= 0)
			return pi - std::atan(std::fabsf(b) / std::fabsf(a));
		if (a <= 0 && b <= 0)
			return -pi + std::atan(std::fabsf(b) / std::fabsf(a));
		if (a >= 0 && b <= 0)
			return -std::atan(std::fabsf(b) / std::fabsf(a));
		return std::atan(std::fabsf(b) / std::fabsf(a));
	}
	bool isReal() const { return b == 0; }
	Com conj() const { return Com(a, -b); }

	Com operator -() const { return Com(-a, -b); }
	Com operator ~() const { return Com(a, -b); }

	Com operator +(Com other) const { return Com(a + other.a, b + other.b); }
	Com operator -(Com other) const { return Com(a - other.a, b - other.b); }
	Com operator *(Com other) const { return Com(a*other.a - b * other.b, a*other.b + other.a*b); }
	Com operator /(Com other) const { return *this*other.conj() / (other.a*other.a + other.b*other.b); }
	Com operator ^(Com other) const
	{
		double r = pow(this->mag(), other.a)*std::exp(-other.b*this->theta());
		double ang = other.a*this->theta() + other.b*std::log(this->mag());
		return Com::make_Com(r, ang);
	}
	Com operator +(double other) const { return Com(a + other, b); }
	Com operator -(double other) const { return Com(a - other, b); }
	Com operator *(double other) const { return Com(a * other, b * other); }
	Com operator /(double other) const { return Com(a / other, b / other); }
	Com operator ^(double other) const
	{
		if (other == 0)
			return 1;
		else if (other == 1)
			return *this;
		else
			return Com::make_Com(std::pow(this->mag(), other), this->theta()*other);
	}

	Com& operator =(Com other) { this->a = other.a; this->b = other.b; return *this; }
	Com& operator+=(Com other) { this->a += other.a; this->b += other.b; return *this; }
	Com& operator-=(Com other) { this->a -= other.a; this->b -= other.b; return *this; }
	Com& operator*=(Com other) { *this = *this * other; return *this; }
	Com& operator/=(Com other) { *this = *this / other; return *this;; }
	Com& operator^=(Com other) { *this = *this ^ other; return *this;; }
	Com& operator =(double other) { this->a = other; this->b = 0; return *this; }
	Com& operator+=(double other) { this->a += other; return *this; }
	Com& operator-=(double other) { this->a -= other; return *this; }
	Com& operator*=(double other) { this->a *= other; this->b *= other; return *this; }
	Com& operator/=(double other) { this->a /= other; this->b /= other; return *this; }
	Com& operator^=(double other) { *this = *this ^ other; return *this; }

	bool operator<=(Com other) const { return mag() <= other.mag(); }
	bool operator< (Com other) const { return mag() <  other.mag(); }
	bool operator>=(Com other) const { return mag() >= other.mag(); }
	bool operator> (Com other) const { return mag() >  other.mag(); }
	bool operator==(Com other) const { return a == other.a && b == other.b; }
	bool operator!=(Com other) const { return a != other.a && b != other.b; }
	bool operator<=(double other) const { return mag() <= other; }
	bool operator< (double other) const { return mag() <  other; }
	bool operator>=(double other) const { return mag() >= other; }
	bool operator> (double other) const { return mag() >  other; }
	bool operator==(double other) const { return a == other && b == 0; }
	bool operator!=(double other) const { return a != other && b != 0; }

	static Com make_Com(double r, double ang) { return Com(r * std::cos(ang), r * std::sin(ang)); }
	friend std::ostream& operator<<(std::ostream& C, Com z);
	void print() const
	{
		if (std::fabsf(a) <= 10e-7 && std::fabsf(b) <= 10e-7)
			printf("0");
		else if (std::fabsf(a) <= 10e-7)
			if (b == 1)
				printf("i");
			else if (b == -1)
				printf("-i");
			else
				printf("%gi", b);
		else if (std::fabsf(b) <= 10e-7)
			printf("%g", a);
		else if (b == 1)
			printf("%g+i", a);
		else if (b == -1)
			printf("%g-i", a);
		else
			printf("%g%+gi", a, b);
		printf(" ");
	}
	~Com() {}
};
const Com i(0, 1);
std::ostream& operator<<(std::ostream& C, Com z) { z.print(); return C; }
std::fstream& operator<<(std::fstream& F, Com z) 
{
	if (std::fabsf(z.a) <= 10e-7 && std::fabsf(z.b) <= 10e-7)
		F << "0";
	else if (std::fabsf(z.a) <= 10e-7)
		if (z.b == 1)
			F << "i";
		else if (z.b == -1)
			F << "-i";
		else
			F << std::to_string(z.b).c_str() << "i";
	else if (std::fabsf(z.b) <= 10e-7)
		F << std::to_string(z.a).c_str();
	else if (z.b == 1)
		F << std::to_string(z.a).c_str() << "+i";
	else if (z.b == -1)
		F << std::to_string(z.a).c_str() << "-i";
	else if (z.b > 0)
		F << std::to_string(z.a).c_str() << "+" << std::to_string(z.b).c_str() << "i";
	else
		F << std::to_string(z.a).c_str() << "-" << std::to_string(z.b).c_str() << "i";
	return F;
}

Com operator+ (double val, const Com z) { return Com(val + z.a, z.b); }
Com operator- (double val, const Com z) { return Com(val - z.a, -z.b); }
Com operator* (double val, const Com z) { return Com(val * z.a, val * z.b); }
Com operator/ (double val, const Com z) { return Com(val) / z; }
Com operator^ (double val, const Com z)
{
	if (val > 0)
		return Com::make_Com(std::pow(val, z.a), z.b*std::log(val));
	if (val == 0)
		return 0;
	return Com(val) ^ z;
}
bool operator<=(double val, const Com z) { return z.mag() <= val; }
bool operator< (double val, const Com z) { return z.mag() <  val; }
bool operator>=(double val, const Com z) { return z.mag() >= val; }
bool operator> (double val, const Com z) { return z.mag() >  val; }
bool operator==(double val, const Com z) { return (z.a == val && z.b == 0); }
bool operator!=(double val, const Com z) { return !(z == val); }
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
namespace com
{
	std::string to_string(Com z)
	{
		std::string a, b = " ";
		if (std::fabsf(z.a) <= 10e-7 && std::fabsf(z.b) <= 10e-7)
			return "0";
		else if (std::fabsf(z.a) <= 10e-7)
			if (z.b == 1)
				return "i";
			else if (z.b == -1)
				return "-i";
			else
				b = std::to_string(z.b);
		else if (std::fabsf(z.b) <= 10e-7)
			a = std::to_string(z.a);
		else if (z.b == 1)
		{
			a = std::to_string(z.a);
			b = "+i";
		}
		else if (z.b == -1)
		{
			a = std::to_string(z.a);
			b = "-i";
		}
		else if (z.b > 0)
		{
			a = std::to_string(z.a);
			b = "+" + std::to_string(z.b);
		}
		else
		{
			a = std::to_string(z.a);
			b = std::to_string(z.b);
		}
		for (register int j = a.length() - 1; j > 0; j--)
			if (a.find('.') != -1 && (a.at(j) == '0' || a.at(j) == '.'))
				a.pop_back();
			else
				break;
		if (!(b == "i" || b == "+i" || b == "-i"))
		{
			for (register int j = b.length() - 1; j > 0; j--)
				if (b.find('.') != -1 && (b.at(j) == '0' || b.at(j) == '.'))
					b.pop_back();
				else
					break;
			b += "i";
		}
		if (b[0] == ' ')
			return a;
		return a + b;
	}

	Com mod(Com z1, Com z2) { return z2.b == 0 ? Com(std::fmodl(z1.a, z2.a), 0) : Com(std::fmodl(z1.a, z2.a), std::fmodl(z1.b, z2.b)); }
	Com abs(Com z) { return Com(std::labs(z.a), std::labs(z.b)); }
	Com ceil(Com z) { return Com(std::ceil(z.a), std::ceil(z.b)); }
	Com floor(Com z) { return Com(std::floor(z.a), std::floor(z.b)); }

	Com exp(Com z) { return Com(std::exp(z.a), z.b); }
	Com pow(Com base, Com exp)
	{
		double r = std::pow(base.mag(), exp.a)*std::exp(-exp.b*base.theta());
		double ang = exp.a*base.theta() + exp.b*std::log(base.mag());
		return Com::make_Com(r, ang);
	}
	Com pow(Com base, double exp) { return Com::make_Com(std::pow(base.mag(), exp), exp * base.theta()); }
	Com pow(double base, Com exp) { return Com::make_Com(std::pow(base, exp.a), exp.b*std::log(base)); }
	Com pow(double base, double exp)
	{
		if (base < 0)
			return pow(Com(base), exp);
		return std::pow(base, exp);
	}
	Com sqrt(Com z) { return pow(z, 0.5); }

	Com ln(Com z) { return Com::make_Com(std::log(z.theta()), z.theta()); }
	Com log(Com z, Com base)
	{
		double den = base.mag() * base.mag() + base.theta() * base.theta();
		double re = std::log(z.mag()) * std::log(base.mag()) + z.theta() * base.theta();
		double im = z.theta() * std::log(base.mag()) - base.theta() * std::log(z.mag());
		return Com(re / den, im / den);
	}
	Com log(Com z, double base = 10)
	{
		double den = base * base;
		double re = std::log(z.mag()) * std::log(base);
		double im = z.theta() * std::log(base);
		return Com(re / den, im / den);
	}
	Com log(double x, Com base)
	{
		double den = base.mag() * base.mag() + base.theta() * base.theta();
		double re = std::log(x) * std::log(base.mag());
		double im = -base.theta() * std::log(x);

		return Com(re / den, im / den);
	}
	Com log(double x, double base = 10)
	{
		if (x > 0 && base > 0)
			return std::log(x) / std::log(base);
		if (x < 0 && base < 0)
			return (pi*i + std::log(-x)) / (pi*i + std::log(-base));
		if (x < 0)
			return (pi*i + std::log(-x)) / std::log(base);
		if (base < 0)
			return std::log(x) / (pi*i + std::log(-base));
		return -INFINITY;
	}

	Com sin(Com z) { return (exp(z*i) - exp(-z * i)) / (2 * i); }
	Com cos(Com z) { return (exp(z*i) + exp(-z * i)) / 2; }
	Com tan(Com z) { return -i * (exp(2 * z*i) - 1) / (exp(2 * z*i) + 1); }
	Com sec(Com z) { return 1 / sin(z); }
	Com csc(Com z) { return 1 / cos(z); }
	Com cot(Com z) { return 1 / tan(z); }

	Com arcsin(Com z) { return -i * ln(i*z + ((1 - (z ^ 2)) ^ 0.5)); }
	Com arccos(Com z) { return -i * ln(z + i * ((1 - (z ^ 2)) ^ 0.5)); }
	Com arctan(Com z) { return i / 2 * ln((i + z) / (i - z)); }
	Com arcsec(Com z) { return arccos(1 / z); }
	Com arccsc(Com z) { return arcsin(1 / z); }
	Com arccot(Com z) { return arctan(1 / z); }

	Com sinh(Com z) { return (exp(z) - exp(-z)) / 2; }
	Com cosh(Com z) { return (exp(z) + exp(-z)) / 2; }
	Com tanh(Com z) { return (exp(2 * z) - 1) / (exp(2 * z) + 1); }
	Com sech(Com z) { return 1 / cosh(z); }
	Com csch(Com z) { return 1 / sinh(z); }
	Com coth(Com z) { return (exp(2 * z) + 1) / (exp(2 * z) - 1); }

	Com arcsinh(Com z) { return ln(z + (((z ^ 2) + 1) ^ 0.5)); }
	Com arccosh(Com z) { return ln(z + (((z ^ 2) - 1) ^ 0.5)); }
	Com arctanh(Com z) { return 0.5*ln((1 + z) / (1 - z)); }
	Com arcsech(Com z) { return arccosh(1 / z); }
	Com arccsch(Com z) { return arcsinh(1 / z); }
	Com arccoth(Com z) { return arctanh(1 / z); }
}
#undef zero_check(x)
void zero_check(int line, const std::string& file_name, Com x, const std::string& errorCode = "Zero Division")
{
	if (x == 0)
	{
		printf("Error\n%s\nError Call from Line No.%d in %s\n", errorCode.c_str(), line, file_name.c_str());
		system("pause");
	}
}
#define zero_check(x) zero_check(__LINE__, fileName, x)

#undef fileName