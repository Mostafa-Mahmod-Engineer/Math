#pragma once
#include<iostream>
#include<fstream>
#include<string>
#define fileName "Utilities.h"
// C O N S T A N T S
const double pi = 3.14159265358979;
const double e = 2.718281828459;
const double phi = 1.61803398874989;
const double phi_bar = -1 / phi;

// F U N C T I O N S

// E R R O R   H A N D L I N G:
void chk_rng(int line, const std::string& file_name, int i, int range, const std::string& errorCode = "Index Exceeds the Range")
{
	if (i >= range)
	{
		printf("Error\n%s\nError Call from Line No.%d in %s\n", errorCode.c_str(), line, file_name.c_str());
		system("pause");
	}
}
void chk_eq(int line, const std::string& file_name, int a, int b, const std::string& errorCode = "Dimensions don't Match")
{
	if (a != b)
	{
		printf("Error\n%s\nError Call from Line No.%d in %s\n", errorCode.c_str(), line, file_name.c_str());
		system("pause");
	}
}
void zero_check(int line, const std::string& file_name, double x, const std::string& errorCode = "Zero Division")
{
	if (x == 0)
	{
		printf("Error\n%s\nError Call from Line No.%d in %s\n", errorCode.c_str(), line, file_name.c_str());
		system("pause");
	}
}
void error(int line, const std::string& file_name, const std::string& errorCode)
{
	std::cout << "Error\n" << errorCode << "\n";
	printf("Error Call from Line No.%d in %s\n", line, file_name.c_str());
	system("pause");
}
#define chk_rng(x1, x2, x3) chk_rng(__LINE__, fileName, x1, x2, x3)
#define chk_rng(x1, x2) chk_rng(__LINE__, fileName, x1, x2)
#define chk_eq(x1, x2, x3) chk_eq(__LINE__, fileName, x1, x2, x3)
#define chk_eq(x1, x2) chk_eq(__LINE__, fileName, x1, x2)
#define zero_check(x1, x2) zero_check(__LINE__, fileName, x1, x2)
#define zero_check(x) zero_check(__LINE__, fileName, x)
#define error(x) error(__LINE__, fileName, x)

// M A T H:
long double Factorial(int x)
{
	if (x == 0)
		return 1;
	long double f = 1;
	for (register int i = 1; i <= x; i++)
		f *= i;
	return f;
}
long double Permutations(int x,  int y)
{
	if (x < y) std::swap(x, y);
	if (x == 0) return 1;
	long double p = 1;
	for (register int i = 0; i < y; i++)
		p *= (x - i);
	return p;
}
long double Combinations(int x, int y)
{
	if (x < y) std::swap(x, y);
	if (y > x / 2) y = x - y;
	if (y == 0 || y == x)return 1;
	if (y == 1) return x;
    return Permutations(x, y) / Factorial(y);
}

float GCD(int first, int second)
{
	if (first == second)
		return first;
	if (second > first)
	{
		second ^= first;
		first ^= second;
		second ^= first;
	}
	int res = first % second;
	while (res != 0)
	{
		first = second;
		second = res;
		res = first % second;
	} return res;
}
bool isPrime(int n)
{
	if (n < 2 && n > -2)
		return false;
	if (n < 0)
		n = -n;
	if (n == 2)
		return true;
	if (n % 2 == 0)
		return false;
	for (register int i = 3; i <= sqrt(n); i += 2)
		if (n%i == 0)
			return false;
	return true;
}

// F O R M A T T I N G:
void pad(std::string& str, int width, char p = ' ')
{
	if (str[0] != '-')
	{
		str = ' ' + str;
		width--;
	}
	for (register int i = 0; i < width; i++)
		str += p;
}


#undef FileName