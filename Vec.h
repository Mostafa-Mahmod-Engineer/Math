#pragma once
#include"Complex.h"
#define fileName "Vec.h"

namespace real
{
	class vec
	{
		double* Vec = nullptr; int size;
	public:
		vec() {}
		explicit vec(int s) : size(s) { Vec = new double[s] {}; }
		template<class T>
		void _assign(T t) { Vec[size - 1] = t; }
		template<class H, class ... T>
		void _assign(H h, T ... t)
	{
		Vec[size - sizeof...(t) - 1] = h;
		_assign(t...);
	}
		template<class H, class ... T>
		vec(H h, T ... t)
	{
		size = sizeof...(T) + 1;
		Vec = new double[size];
		_assign(h, t...);
	}
		vec(const vec& other)
	{
		size = other.size;
		Vec = new double[size];
		for (register int i = 0; i < size; i++)
			Vec[i] = other.Vec[i];
	}
		vec(vec&& other)
	{
		size = other.size;
		Vec = other.Vec;
		other.Vec = nullptr;
	}
		vec(const vec* other)
	{
		size = other->size;
		Vec = new double[size];
		for (register int i = 0; i < size; i++)
			Vec[i] = other->Vec[i];
	}
		vec(vec* other)
	{
		size = other->size;
		Vec = new double[size];
		for (register int i = 0; i < size; i++)
			Vec[i] = other->Vec[i];
	}
		//=============================== << G E T T E R S >> ====================================
		int Size() const { return size; }
		double*begin() const { return Vec; }
		double*end() const { return Vec + sizeof(double)*(size - 1); }
		double operator[](int i) const { chk_rng(i, size); return Vec[i]; }
		double sum() const
	{
		double sum = 0;
		for (register int i = 0; i < size; i++)
			sum += Vec[i];
		return sum;
	}
		double Psum(double p) const
	{
		double sum = 0;
		for (register int i = 0; i < size; i++)
			sum += std::pow(Vec[i], p);
		return sum;
	}
		double Wsum(const vec& other) const
	{
		chk_eq(other.size, size);
		double sum = 0;
		for (register int i = 0; i < size; i++)
			sum += Vec[i] * other.Vec[i];
		return sum;
	}
		double min() const
	{
		double min = Vec[0];
		for (register int i = 0; i < size; i++)
			if (Vec[i] < min)
				min = Vec[i];
		return min;
	}
		double max() const
	{
		double max = Vec[0];
		for (register int i = 0; i < size; i++)
			if (Vec[i] > max)
				max = Vec[i];
		return max;
	}
		double mini() const
	{
		double min = Vec[0]; int ind = 0;
		for (register int i = 0; i < size; i++)
			if (Vec[i] < min)
			{
				min = Vec[i];
				ind = i;
			}
		return ind;
	}
		double maxi() const
	{
		double max = Vec[0]; int ind = 0;
		for (register int i = 0; i < size; i++)
			if (Vec[i] > max)
			{
				max = Vec[i];
				ind = i;
			}
		return ind;
	}
		double abs_min() const
	{
		double min = Vec[0];
		for (register int i = 0; i < size; i++)
			if (std::fabs(Vec[i]) < std::fabs(min))
				min = Vec[i];
		return min;
	}
		double abs_max() const
	{
		double max = Vec[0];
		for (register int i = 0; i < size; i++)
			if (std::fabs(Vec[i]) > std::fabs(max))
				max = Vec[i];
		return max;
	}
		double abs_mini() const
	{
		double min = Vec[0]; int ind = 0;
		for (register int i = 0; i < size; i++)
			if (std::fabs(Vec[i]) < std::fabs(min))
			{
				min = Vec[i];
				ind = i;
			}
		return ind;
	}
		double abs_maxi() const
	{
		double max = Vec[0]; int ind = 0;
		for (register int i = 0; i < size; i++)
			if (std::fabs(Vec[i]) > std::fabs(max))
			{
				max = Vec[i];
				ind = i;
			}
		return ind;
	}
		double mag() const
	{
		double sum = 0;
		for (register int i = 0; i < size; i++)
			sum += Vec[i] * Vec[i];
		return std::sqrt(sum);
	}
		double magsqr() const
	{
		double sum = 0;
		for (register int i = 0; i < size; i++)
			sum += Vec[i] * Vec[i];
		return sum;
	}
		double Norm() const
	{
		double sum = 0;
		for (register int i = 0; i < size; i++)
			sum += std::fabs(Vec[i]);
		return sum;
	}
		double E_Norm() const
	{
		double sum = 0;
		for (register int i = 0; i < size; i++)
			sum += Vec[i] * Vec[i];
		return std::sqrt(sum);
	}
		double P_Norm(int p) const
	{
		double sum = 0;
		for (register int i = 0; i < size; i++)
			sum += std::pow(Vec[i], p);
		return std::pow(sum, 1.0 / p);
	}
		double inf_Norm() const { return abs_max(); }
		double range() const { return max() - min(); }
		double mean() const { return sum() / size; }
		double median() const
	{
		vec res = sort();
		return size % 2 == 1 ? res.Vec[size / 2] : (res.Vec[size / 2 - 1] + res.Vec[size / 2]) / 2;
	}
		double variance() const
	{
		double miu = mean(), sum = 0;
		for (register int i = 0; i < size; i++)
			sum += (Vec[i] - miu)*(Vec[i] - miu);
		return sum / size;
	}
		double standard_deviation() const { return std::sqrt(variance()); }
		double Q1() const { return median() - min(); }
		double Q3() const { return max() - median(); }
		double iqr()/*InterQuartile Range*/ const { return max() - 2 * median() + min(); }
		bool is_sorted() const
	{
		int c = 0;
		while (Vec[c] == Vec[c + 1] && c < size - 1)
			c++;
		if (Vec[c] > Vec[c + 1])
			for (register int i = c + 2; i < size; i++)
				if (Vec[i - 1] < Vec[i])
					return false;
				else
					continue;
		else
			for (register int i = c + 2; i < size; i++)
				if (Vec[i - 1] > Vec[i])
					return false;
		return true;
	}
		//=============================== << G E T T E R S \E.N.D>> ===============================
	
		//=============================== << F U N C T I O N S >> =================================
		/* FIX THE COPYING PROBLEM */
		void _merge(vec& res, int start, int mid, int end)
	{
		double* temp;
		int i = start, j = mid + 1, k = start;
		while (i <= mid && j <= end)
			if (Vec[i] > Vec[j])
			{
				res.Vec[k] = Vec[j];
				j++; k++;
			}
			else
			{
				res.Vec[k] = Vec[i];
				i++; k++;
			}
		while (i <= mid)
		{
			res.Vec[k] = Vec[i];
			i++; k++;
		}
		while (j <= end)
		{
			res.Vec[k] = Vec[j];
			j++; k++;
		}
		// swap the 2 vectors to avoid unnecessary copy
		temp = res.Vec;
		res.Vec = Vec;
		Vec = temp;
	}
		void _sort(int start = 0, int end = -1)
	{
		vec res(this);
		if (end == -1)
			end = size - 1;
		if (start < end)
		{
			int mid = (start + end) / 2;
			_sort(start, mid);
			_sort(mid + 1, end);
			_merge(res, start, mid, end);
		}
		return;
	}
		void _reverse()
	{
		double temp;
		for (register int i = 0; i < size / 2; i++)
		{
			temp = Vec[i];
			Vec[i] = Vec[size - i - 1];
			Vec[size - i - 1] = temp;
		}
	}
		void _shift_right(int n) { Vec = shift_right(n).Vec; }
		void _shift_left(int n) { Vec = shift_left(n).Vec; }
		void _hat()
	{
		double magnitude = mag();
		for (register int i = 0; i < size; i++)
			Vec[i] /= magnitude;
	}
		void _norm()
	{
		double Max = abs_max();
		zero_check(Max, "Normalizing a Zero Vector");
		for (register int i = 0; i < size; i++)
			Vec[i] /= Max;
	}
	
		vec trunc(int start, int end) const
	{
		chk_rng(end, size);
		vec res(end - start);
		for (register int i = start; i < end; i++)
			res.Vec[i - start] = Vec[i];
		return res;
	}
		vec trunc(int end) const
	{
		chk_rng(end, size);
		vec res(end);
		for (register int i = 0; i < end; i++)
			res.Vec[i] = Vec[i];
		return res;
	}
		vec hat() const { return *this / mag(); }
		vec norm() const { return *this / abs_max(); }
		vec deviations() const
	{
		vec res(size); double miu = mean();
		for (register int i = 0; i < size; i++)
			res.Vec[i] = Vec[i] - miu;
		return res;
	}
		vec adjacent_differences() const
	{
		vec res(size - 1);
		for (register int i = 1; i < size; i++)
			res.Vec[i - 1] = Vec[i] - Vec[i - 1];
		return res;
	}
		vec adjacent_ratios() const
	{
		vec res(size - 1);
		for (register int i = 1; i < size; i++)
		{
			zero_check(Vec[i - 1], "Call from adjacent_ratios()");
			res.Vec[i - 1] = Vec[i] / Vec[i - 1];
		}
		return res;
	}
		vec sort() const
	{
		vec res(this);
		res._sort();
		return res;
	}
		vec reverse() const
	{
		vec res(this);
		res._reverse();
		return res;
	}
		vec shift_right(int n) const
	{
		vec res(size);
		for (register int i = 0; i < size; i++)
			res.Vec[(i + n) % size] = Vec[i];
		return res;
	}
		vec shift_left(int n) const
	{
		vec res(size);
		for (register int i = 0; i < size; i++)
			res.Vec[i] = Vec[(i + n) % size];
		return res;
	}
		int find(double val) const
	{
		for (register int i = 0; i < size; i++)
			if (Vec[i] == val)
				return i;
		return -1;
	}
	
		static vec O(int n)
	{
		vec res(n);
		for (register int i = 0; i < n; i++)
			res.Vec[i] = 1;
		return res;
	}
		static vec R(int n)
	{
		vec res(n);
		for (register int i = 0; i < n; i++)
			res.Vec[i] = (std::rand() % 20) - 10;
		return res;
	}
		static vec R(int n, double upperlimit)
	{
		vec res(n);
		for (register int i = 0; i < n; i++)
			res.Vec[i] = std::fmodf(std::rand(), upperlimit);
		return res;
	}
		static vec R(int n, double lowerlimit, double upperlimit)
	{
		vec res(n);
		for (register int i = 0; i < n; i++)
			res.Vec[i] = std::fmodf(std::rand(), upperlimit - lowerlimit) + lowerlimit;
		return res;
	}
		//=============================== << F U N C T I O N S \E.N.D>> ===========================
	
		//=============================== << O P E R A T O R S >> =================================
		vec operator-() const
	{
		vec res(size);
		for (register int i = 0; i < size; i++)
			res.Vec[i] = -Vec[i];
		return res;
	}
		//================================    VEC X VEC    =====================================
		vec operator+(const vec& other) const
	{
		chk_eq(size, other.size);
		vec res(size);
		for (register int i = 0; i < size; i++)
			res.Vec[i] = Vec[i] + other.Vec[i];
		return res;
	}
		vec operator-(const vec& other) const
	{
		chk_eq(size, other.size);
		vec res(size);
		for (register int i = 0; i < size; i++)
			res.Vec[i] = Vec[i] - other.Vec[i];
		return res;
	}
		vec operator*(const vec& other) const
	{
		chk_eq(size, other.size);
		vec res(size);
		for (register int i = 0; i < size; i++)
			res.Vec[i] = Vec[i] * other.Vec[i];
		return res;
	}
		vec operator/(const vec& other) const
	{
		chk_eq(size, other.size);
		vec res(size);
		for (register int i = 0; i < size; i++)
		{
			zero_check(other.Vec[i]);
			res.Vec[i] = Vec[i] / other.Vec[i];
		}
		return res;
	}
		array<bool> operator==(const vec& other) const
	{
		chk_eq(size, other.size);
		array<bool> res(size);
		for (register int i = 0; i < size; i++)
			res.push_back(Vec[i] == other.Vec[i]);
		return res;
	}
		array<bool> operator!=(const vec& other) const
	{
		chk_eq(size, other.size);
		array<bool> res(size);
		for (register int i = 0; i < size; i++)
			res.push_back(Vec[i] != other.Vec[i]);
		return res;
	}
		array<bool> operator> (const vec& other) const
	{
		chk_eq(size, other.size);
		array<bool> res(size);
		for (register int i = 0; i < size; i++)
			res.push_back(Vec[i] > other.Vec[i]);
		return res;
	}
		array<bool> operator< (const vec& other) const
	{
		chk_eq(size, other.size);
		array<bool> res(size);
		for (register int i = 0; i < size; i++)
			res.push_back(Vec[i] < other.Vec[i]);
		return res;
	}
		array<bool> operator>=(const vec& other) const
	{
		chk_eq(size, other.size);
		array<bool> res(size);
		for (register int i = 0; i < size; i++)
			res.push_back(Vec[i] >= other.Vec[i]);
		return res;
	}
		array<bool> operator<=(const vec& other) const
	{
		chk_eq(size, other.size);
		array<bool> res(size);
		for (register int i = 0; i < size; i++)
			res.push_back(Vec[i] <= other.Vec[i]);
		return res;
	}
	
		vec& operator=(const vec& other)
	{
		if (Vec == other.Vec)
			return *this;
		this->~vec();
		size = other.size;
		Vec = new double[size];
		for (register int i = 0; i < size; i++)
			Vec[i] = other.Vec[i];
		return *this;
	}
		vec& operator=(vec&& other)
	{
		size = other.size;
		Vec = other.Vec;
		other.Vec = nullptr;
		return *this;
	}
		vec& operator+=(const vec& other)
	{
		*this = *this + other;
		return *this;
	}
		vec& operator-=(const vec& other)
	{
		*this = *this - other;
		return *this;
	}
		vec& operator*=(const vec& other)
	{
		*this = *this * other;
		return *this;
	}
		vec& operator/=(const vec& other)
	{
		*this = *this / other;
		return *this;
	}
		//===================================    VEC X NUM    =====================================
		vec  operator+(double val) const
	{
		vec res(size);
		for (register int i = 0; i < size; i++)
			res.Vec[i] = Vec[i] + val;
		return res;
	}
		vec  operator-(double val) const
	{
		vec res(size);
		for (register int i = 0; i < size; i++)
			res.Vec[i] = Vec[i] - val;
		return res;
	}
		vec  operator*(double val) const
	{
		vec res(size);
		for (register int i = 0; i < size; i++)
			res.Vec[i] = Vec[i] * val;
		return res;
	}
		vec  operator/(double val) const
	{
		zero_check(val);
		vec res(size);
		for (register int i = 0; i < size; i++)
			res.Vec[i] = Vec[i] / val;
		return res;
	}
		vec  operator^(double val) const
	{
		vec res(size);
		for (register int i = 0; i < size; i++)
			res.Vec[i] = std::pow(Vec[i], val);
		return res;
	}
		array<bool> operator==(double val) const
	{
		array<bool> res(size);
		for (register int i = 0; i < size; i++)
			res.push_back(Vec[i] == val);
		return res;
	}
		array<bool> operator!=(double val) const
	{
		array<bool> res(size);
		for (register int i = 0; i < size; i++)
			res.push_back(Vec[i] != val);
		return res;
	}
		array<bool> operator> (double val) const
	{
		array<bool> res(size);
		for (register int i = 0; i < size; i++)
			res.push_back(Vec[i] > val);
		return res;
	}
		array<bool> operator< (double val) const
	{
		array<bool> res(size);
		for (register int i = 0; i < size; i++)
			res.push_back(Vec[i] < val);
		return res;
	}
		array<bool> operator>=(double val) const
	{
		array<bool> res(size);
		for (register int i = 0; i < size; i++)
			res.push_back(Vec[i] >= val);
		return res;
	}
		array<bool> operator<=(double val) const
	{
		array<bool> res(size);
		for (register int i = 0; i < size; i++)
			res.push_back(Vec[i] <= val);
		return res;
	}
	
		vec& operator+=(double val)
	{
		*this = *this + val;
		return *this;
	}
		vec& operator-=(double val)
	{
		*this = *this - val;
		return *this;
	}
		vec& operator*=(double val)
	{
		*this = *this * val;
		return *this;
	}
		vec& operator/=(double val)
	{
		*this = *this / val;
		return *this;
	}
		vec& operator^=(double val)
	{
		*this = *this ^ val;
		return *this;
	}
		//=============================== << O P E R A T O R S \E.N.D>> ===================================
		void set(int i, double val)
	{
		chk_rng(i, size);
		if (i < 0)
			error("Negative Index");
		Vec[i] = val;
	}
		void set(double val)
	{
		for (register int i = 0; i < size; i++)
			Vec[i] = val;
	}
		void reset(int s)
	{
		if (Vec != nullptr) delete[] Vec;
		size = s;
		Vec = new double[size] {};
	}
		void reset()
	{
		delete[] Vec;
		Vec = new double[size] {};
	}
		void print() const
	{
		printf("\n(");
		for (register int i = 0; i < size - 1; i++)
			printf("%g,", Vec[i]);
		printf("%g)\n", Vec[size - 1]);
	}
		void erase() { Vec = nullptr; }
		~vec() { if (Vec != nullptr) delete[] Vec; }
	};
	vec operator*(double val, const vec& v) { return v * val; }
	vec operator/(double val, const vec& v)
{
	int size = v.Size();
	vec res(size);
	for (register int i = 0; i < size; i++)
	{
		zero_check(v[i]);
		res.set(i, 1 / v[i]);
	}
	return res;
}
	//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
	double sum(const vec& v)
{
	double sum = 0; int size = v.Size();
	for (register int i = 0; i < size; i++)
		sum += v[i];
	return sum;
}
	double Wsum(const vec& v1, const vec& v2)
{
	chk_eq(v1.Size(), v2.Size());
	int size = v1.Size();
	double sum = 0;
	for (register int i = 0; i < size; i++)
		sum += v1[i] * v2[i];
	return sum;
}
	double r_value(const vec& x, const vec& y)
{
	chk_eq(x.Size(), y.Size());
	int n = x.Size();
	double sumxy = Wsum(x, y), sumx = x.sum(), sumy = y.sum(), sumx2 = x.magsqr(), sumy2 = y.magsqr();
	return (n*sumxy - sumx * sumy) / std::sqrt((n*sumx2 - sumx * sumx)*(n*sumy2 - sumy * sumy));
}
	double prod(const vec& v)
{
	double prod = 1; int size = v.Size();
	for (register int i = 0; i < size; i++)
		prod *= v[i];
	return prod;
}
	double mag(const vec& v)
{
	double sum = 0; int size = v.Size();
	for (register int i = 0; i < size; i++)
		sum += v[i] * v[i];
	return std::sqrt(sum);
}
	double min(const vec& v)
{
	double min = v[0]; int size = v.Size();
	for (register int i = 0; i < size; i++)
		if (v[i] < min)
			min = v[i];
	return min;
}
	double max(const vec& v)
{
	double max = v[0]; int size = v.Size();
	for (register int i = 0; i < size; i++)
		if (v[i] > max)
			max = v[i];
	return max;
}
	int mini(const vec& v)
{
	double min = v[0]; int size = v.Size(), ind = 0;
	for (register int i = 0; i < size; i++)
		if (v[i] < min)
		{
			min = v[i];
			ind = i;
		}
	return ind;
}
	int maxi(const vec& v)
{
	double max = v[0]; int size = v.Size(), ind = 0;
	for (register int i = 0; i < size; i++)
		if (v[i] > max)
		{
			max = v[i];
			ind = i;
		}
	return ind;
}
	vec reverse(const vec& v)
{
	int size = v.Size();
	vec res(size);
	for (register int i = 0; i < size; i++)
		res.set(i, v[size - i - 1]);
	return res;
}
	vec hat(const vec& v) { return v / mag(v); }
	vec norm(const vec& v) { return v / max(v); }
	vec _set(const vec& v, double(*F)(double x))
{
	int size = v.Size();
	vec res(size);
	for (register int i = 0; i < size; i++)
		res.set(i, (*F)(v[i]));
	return res;
}
	vec abs(const vec& v) { return _set(v, std::fabs); }
	vec ceil(const vec& v) { return _set(v, std::ceil); }
	vec floor(const vec& v) { return _set(v, std::floor); }
	vec mod(const vec& v, double y)
{
	int size = v.Size();
	vec res(size);
	for (register int i = 0; i < size; i++)
		res.set(i, std::fmod(v[i], y));
	return res;
}

	vec exp(const vec& v) { return _set(v, std::exp); }
	vec pow(const vec& v, double ex)
{
	int size = v.Size();
	vec res(size);
	for (register int i = 0; i < size; i++)
		res.set(i, std::pow(v[i], ex));
	return res;
}
	vec ln(const vec& v) { return _set(v, std::log); }
	vec log(const vec& v, double base = 10)
{
	int size = v.Size();
	vec res(size);
	for (register int i = 0; i < size; i++)
	{
		chk_rng(-v[i], 0);
		res.set(i, std::log(v[i]) / std::log(10));
	}
	return res;
}
	
	vec sin(const vec& v) { return _set(v, std::sin); }
	vec cos(const vec& v) { return _set(v, std::cos); }
	vec tan(const vec& v) { return _set(v, std::tan); }

	vec sec(const vec& v) { return _set(v, [](double val)->double {return 1 / std::cos(val); }); }
	vec csc(const vec& v) { return _set(v, [](double val)->double {return 1 / std::sin(val); }); }
	vec cot(const vec& v) { return _set(v, [](double val)->double {return 1 / std::tan(val); }); }

	vec arcsin(const vec& v) { return _set(v, std::asin); }
	vec arccos(const vec& v) { return _set(v, std::acos); }
	vec arctan(const vec& v) { return _set(v, std::atan); }

	vec arcsec(const vec& v) { return _set(v, [](double val)->double {return std::acos(1 / val); }); }
	vec arccsc(const vec& v) { return _set(v, [](double val)->double {return std::asin(1 / val); }); }
	vec arccot(const vec& v) { return _set(v, [](double val)->double {return std::atan(1 / val); }); }
	
	vec sinh(const vec& v) { return _set(v, std::sinh); }
	vec cosh(const vec& v) { return _set(v, std::cosh); }
	vec tanh(const vec& v) { return _set(v, std::tanh); }

	vec sech(const vec& v) { return _set(v, [](double val)->double {return std::cosh(val); }); }
	vec csch(const vec& v) { return _set(v, [](double val)->double {return std::sinh(val); }); }
	vec coth(const vec& v) { return _set(v, [](double val)->double {return std::tanh(val); }); }

	vec arcsinh(const vec& v) { return _set(v, asinh); }
	vec arccosh(const vec& v) { return _set(v, acosh); }
	vec arctanh(const vec& v) { return _set(v, atanh); }

	vec arcsech(const vec& v) { return _set(v, [](double val)->double {return std::acosh(1 / val); }); }
	vec arccsch(const vec& v) { return _set(v, [](double val)->double {return std::asinh(1 / val); }); }
	vec arccoth(const vec& v) { return _set(v, [](double val)->double {return std::atanh(1 / val); }); }

	std::ostream& operator<<(std::ostream& C, const vec& v) { v.print(); return C; }
	
	vec ones(int s) { return _set(vec(s), [](double val)->double { return 1; }); }
	vec rand(int s) { return _set(vec(s), [](double val)->double { return std::rand() % 10; }); }
	vec ar_seq(double a, double b, int n)
{
	vec res(n);
	res.set(0, a);
	for (register int i = 1; i < n; i++)
		res.set(i, res[i - 1] + b);
	return res;
}
	vec geo_seq(double a, double r, int n)
{
	vec res(n);
	res.set(0, a);
	for (register int i = 1; i < n; i++)
		res.set(i, res[i - 1] * r);
	return res;
}
	vec fib_seq(double f1, double f2, int n)
{
	vec res(n);
	res.set(0, f1); res.set(1, f2);
	for (register int i = 2; i < n; i++)
		res.set(i, res[i - 1] + res[i - 2]);
	return res;
}
	vec linspace(double start, double end, int n) { return ar_seq(start, (end - start) / n, n); }
	vec seq(double(*f)(double*), const vec& args, int n)
{
	vec res(n);
	int num_of_args = args.Size();
	chk_rng(num_of_args, n, "Number of Initial Values Exceeds the Number of the Members of the Vector\nCall from seq");
	for (register int i = 0; i < num_of_args; i++)
		res.set(i, args[i]);
	for (register int i = num_of_args; i < n; i++)
		res.set(i, (*f)(res.begin() + i - num_of_args));
	return res;
}
	vec seq(double(*f)(double), double arg, int n)
{
	vec res(n);
	res.set(0, arg);
	for (register int i = 1; i < n; i++)
		res.set(i, (*f)(res[i - 1]));
	return res;
}
	vec seq(double(*f)(double), int n, double start = 1, double step = 1)
{
	vec res(n);
	for (register int i = 0; i < n; i++)
		res.set(i, (*f)(start + i * step));
	return res;
}
}

namespace com
{
	class vec
	{
		Com* Vec = nullptr; int size;
	public:
		vec() {}
		explicit vec(int s) : size(s) { Vec = new Com[s] {}; }
		template<class T>
		void _assign(T t) { Vec[size - 1] = t; }
		template<class H, class ... T>
		void _assign(H h, T ... t)
		{
			Vec[size - sizeof...(t) - 1] = h;
			_assign(t...);
		}
		template<class H, class ... T>
		vec(H h, T ... t)
		{
			size = sizeof...(T) + 1;
			Vec = new double[size];
			_assign(h, t...);
		}
		vec(const vec& other)
		{
			size = other.size;
			Vec = new Com[size];
			for (register int i = 0; i < size; i++)
				Vec[i] = other.Vec[i];
		}
		vec(vec&& other)
		{
			size = other.size;
			Vec = other.Vec;
			other.Vec = nullptr;
		}
		vec(const vec* other)
		{
			size = other->size;
			Vec = new Com[size];
			for (register int i = 0; i < size; i++)
				Vec[i] = other->Vec[i];
		}
		vec(vec* other)
		{
			size = other->size;
			Vec = new Com[size];
			for (register int i = 0; i < size; i++)
				Vec[i] = other->Vec[i];
		}
		//=============================== << G E T T E R S >> ====================================
		int Size() const { return size; }
		Com*begin() const { return Vec; }
		Com*end() const { return Vec + sizeof(double)*(size - 1); }
		Com operator[](int i) const { chk_rng(i, size); return Vec[i]; }
		Com sum() const
		{
			Com sum = 0;
			for (register int i = 0; i < size; i++)
				sum += Vec[i];
			return sum;
		}
		Com Psum(double p) const
		{
			Com sum = 0;
			for (register int i = 0; i < size; i++)
				sum += Vec[i] ^ p;
			return sum;
		}
		Com Wsum(const vec& other) const
		{
			chk_eq(other.size, size);
			Com sum = 0;
			for (register int i = 0; i < size; i++)
				sum += Vec[i] * other.Vec[i];
			return sum;
		}
		Com min() const
		{
			Com min = Vec[0];
			for (register int i = 0; i < size; i++)
				if (Vec[i] < min)
					min = Vec[i];
			return min;
		}
		Com max() const
		{
			Com max = Vec[0];
			for (register int i = 0; i < size; i++)
				if (Vec[i] > max)
					max = Vec[i];
			return max;
		}
		Com mini() const
		{
			Com min = Vec[0]; int ind = 0;
			for (register int i = 0; i < size; i++)
				if (Vec[i] < min)
				{
					min = Vec[i];
					ind = i;
				}
			return ind;
		}
		Com maxi() const
		{
			Com max = Vec[0]; int ind = 0;
			for (register int i = 0; i < size; i++)
				if (Vec[i] > max)
				{
					max = Vec[i];
					ind = i;
				}
			return ind;
		}
		Com mag() const
		{
			Com sum = 0;
			for (register int i = 0; i < size; i++)
				sum += Vec[i] * Vec[i];
			return com::sqrt(sum);
		}
		Com magsqr() const
		{
			Com sum = 0;
			for (register int i = 0; i < size; i++)
				sum += Vec[i] * Vec[i];
			return sum;
		}
		Com Norm() const
		{
			Com sum = 0;
			for (register int i = 0; i < size; i++)
				sum += com::abs(Vec[i]);
			return sum;
		}
		Com E_Norm() const
		{
			Com sum = 0;
			for (register int i = 0; i < size; i++)
				sum += Vec[i] * Vec[i];
			return com::sqrt(sum);
		}
		Com P_Norm(int p) const
		{
			Com sum = 0;
			for (register int i = 0; i < size; i++)
				sum += com::pow(Vec[i], p);
			return com::pow(sum, 1.0 / p);
		}
		Com inf_Norm() const { return max(); }
		Com range() const { return max() - min(); }
		Com mean() const { return sum() / size; }
		Com median() const
		{
			vec res = sort();
			return size % 2 == 1 ? res.Vec[size / 2] : (res.Vec[size / 2 - 1] + res.Vec[size / 2]) / 2;
		}
		Com variance() const
		{
			Com miu = mean(), sum = 0;
			for (register int i = 0; i < size; i++)
				sum += (Vec[i] - miu)*(Vec[i] - miu);
			return sum / size;
		}
		Com standard_deviation() const { return com::sqrt(variance()); }
		Com Q1() const { return median() - min(); }
		Com Q3() const { return max() - median(); }
		Com iqr()/*InterQuartile Range*/ const { return max() - 2 * median() + min(); }
		bool is_sorted() const
		{
			int c = 0;
			while (Vec[c] == Vec[c + 1] && c < size - 1)
				c++;
			if (Vec[c] > Vec[c + 1])
				for (register int i = c + 2; i < size; i++)
					if (Vec[i - 1] < Vec[i])
						return false;
					else
						continue;
			else
				for (register int i = c + 2; i < size; i++)
					if (Vec[i - 1] > Vec[i])
						return false;
			return true;
		}
		//=============================== << G E T T E R S \E.N.D>> ===============================

		//=============================== << F U N C T I O N S >> =================================
		vec trunc(int start, int end) const
		{
			chk_rng(end, size);
			vec res(end - start);
			for (register int i = start; i < end; i++)
				res.Vec[i - start] = Vec[i];
			return res;
		}
		vec trunc(int end) const
		{
			chk_rng(end, size);
			vec res(end);
			for (register int i = 0; i < end; i++)
				res.Vec[i] = Vec[i];
			return res;
		}
		vec hat() const { return *this / mag(); }
		vec norm() const { return *this / max(); }
		vec deviations() const
		{
			vec res(size); Com miu = mean();
			for (register int i = 0; i < size; i++)
				res.Vec[i] = Vec[i] - miu;
			return res;
		}
		vec adjacent_differences() const
		{
			vec res(size - 1);
			for (register int i = 1; i < size; i++)
				res.Vec[i - 1] = Vec[i] - Vec[i - 1];
			return res;
		}
		vec adjacent_ratios() const
		{
			vec res(size - 1);
			for (register int i = 1; i < size; i++)
			{
				zero_check(Vec[i - 1], "Call from adjacent_ratios()");
				res.Vec[i - 1] = Vec[i] / Vec[i - 1];
			}
			return res;
		}
		vec sort() const
		{
			vec res(this);
			res._sort();
			return res;
		}
		vec reverse() const
		{
			vec res(this);
			res._reverse();
			return res;
		}
		vec shift_right(int n) const
		{
			vec res(size);
			for (register int i = 0; i < size; i++)
				res.Vec[(i + n) % size] = Vec[i];
			return res;
		}
		vec shift_left(int n) const
		{
			vec res(size);
			for (register int i = 0; i < size; i++)
				res.Vec[i] = Vec[(i + n) % size];
			return res;
		}
		int find(double val) const
		{
			for (register int i = 0; i < size; i++)
				if (Vec[i] == val)
					return i;
			return -1;
		}

		/* FIX THE COPYING PROBLEM */
		void _merge(vec& res, int start, int mid, int end)
		{
			Com* temp;
			int i = start, j = mid + 1, k = start;
			while (i <= mid && j <= end)
				if (Vec[i] > Vec[j])
				{
					res.Vec[k] = Vec[j];
					j++; k++;
				}
				else
				{
					res.Vec[k] = Vec[i];
					i++; k++;
				}
			while (i <= mid)
			{
				res.Vec[k] = Vec[i];
				i++; k++;
			}
			while (j <= end)
			{
				res.Vec[k] = Vec[j];
				j++; k++;
			}
			// swap the 2 vectors to avoid unnecessary copy
			temp = res.Vec;
			res.Vec = Vec;
			Vec = temp;
		}
		void _sort(int start = 0, int end = -1)
		{
			vec res(this);
			if (end == -1)
				end = size - 1;
			if (start < end)
			{
				int mid = (start + end) / 2;
				_sort(start, mid);
				_sort(mid + 1, end);
				_merge(res, start, mid, end);
			}
			return;
		}
		void _reverse()
		{
			Com temp;
			for (register int i = 0; i < size / 2; i++)
			{
				temp = Vec[i];
				Vec[i] = Vec[size - i - 1];
				Vec[size - i - 1] = temp;
			}
		}
		void _shift_right(int n) { Vec = shift_right(n).Vec; }
		void _shift_left(int n) { Vec = shift_left(n).Vec; }
		void _hat()
		{
			Com magnitude = mag();
			for (register int i = 0; i < size; i++)
				Vec[i] /= magnitude;
		}
		void _norm()
		{
			Com Max = max();
			zero_check(Max, "Normalizing a Zero Vector");
			for (register int i = 0; i < size; i++)
				Vec[i] /= Max;
		}

		static vec O(int n)
		{
			vec res(n);
			for (register int i = 0; i < n; i++)
				res.Vec[i] = 1;
			return res;
		}
		static vec R(int n)
		{
			vec res(n);
			for (register int i = 0; i < n; i++)
				res.Vec[i] = (std::rand() % 20) - 10;
			return res;
		}
		static vec R(int n, double upperlimit)
		{
			vec res(n);
			for (register int i = 0; i < n; i++)
				res.Vec[i] = std::fmodf(std::rand(), upperlimit);
			return res;
		}
		static vec R(int n, double lowerlimit, double upperlimit)
		{
			vec res(n);
			for (register int i = 0; i < n; i++)
				res.Vec[i] = std::fmodf(std::rand(), upperlimit - lowerlimit) + lowerlimit;
			return res;
		}
		//=============================== << F U N C T I O N S \E.N.D>> ===========================

		//=============================== << O P E R A T O R S >> =================================
		vec operator-() const
		{
			vec res(size);
			for (register int i = 0; i < size; i++)
				res.Vec[i] = -Vec[i];
			return res;
		}
		//================================    VEC X VEC    =====================================
		vec operator+(const vec& other) const
		{
			chk_eq(size, other.size);
			vec res(size);
			for (register int i = 0; i < size; i++)
				res.Vec[i] = Vec[i] + other.Vec[i];
			return res;
		}
		vec operator-(const vec& other) const
		{
			chk_eq(size, other.size);
			vec res(size);
			for (register int i = 0; i < size; i++)
				res.Vec[i] = Vec[i] - other.Vec[i];
			return res;
		}
		vec operator*(const vec& other) const
		{
			chk_eq(size, other.size);
			vec res(size);
			for (register int i = 0; i < size; i++)
				res.Vec[i] = Vec[i] * other.Vec[i];
			return res;
		}
		vec operator/(const vec& other) const
		{
			chk_eq(size, other.size);
			vec res(size);
			for (register int i = 0; i < size; i++)
			{
				zero_check(other.Vec[i]);
				res.Vec[i] = Vec[i] / other.Vec[i];
			}
			return res;
		}
		array<bool> operator==(const vec& other) const
		{
			chk_eq(size, other.size);
			array<bool> res(size);
			for (register int i = 0; i < size; i++)
				res.push_back(Vec[i] == other.Vec[i]);
			return res;
		}
		array<bool> operator!=(const vec& other) const
		{
			chk_eq(size, other.size);
			array<bool> res(size);
			for (register int i = 0; i < size; i++)
				res.push_back(Vec[i] != other.Vec[i]);
			return res;
		}
		array<bool> operator> (const vec& other) const
		{
			chk_eq(size, other.size);
			array<bool> res(size);
			for (register int i = 0; i < size; i++)
				res.push_back(Vec[i] > other.Vec[i]);
			return res;
		}
		array<bool> operator< (const vec& other) const
		{
			chk_eq(size, other.size);
			array<bool> res(size);
			for (register int i = 0; i < size; i++)
				res.push_back(Vec[i] < other.Vec[i]);
			return res;
		}
		array<bool> operator>=(const vec& other) const
		{
			chk_eq(size, other.size);
			array<bool> res(size);
			for (register int i = 0; i < size; i++)
				res.push_back(Vec[i] >= other.Vec[i]);
			return res;
		}
		array<bool> operator<=(const vec& other) const
		{
			chk_eq(size, other.size);
			array<bool> res(size);
			for (register int i = 0; i < size; i++)
				res.push_back(Vec[i] <= other.Vec[i]);
			return res;
		}

		vec& operator=(const vec& other)
		{
			if (Vec == other.Vec)
				return *this;
			this->~vec();
			size = other.size;
			Vec = new Com[size];
			for (register int i = 0; i < size; i++)
				Vec[i] = other.Vec[i];
			return *this;
		}
		vec& operator=(vec&& other)
		{
			size = other.size;
			Vec = other.Vec;
			other.Vec = nullptr;
			return *this;
		}
		vec& operator+=(const vec& other)
		{
			*this = *this + other;
			return *this;
		}
		vec& operator-=(const vec& other)
		{
			*this = *this - other;
			return *this;
		}
		vec& operator*=(const vec& other)
		{
			*this = *this * other;
			return *this;
		}
		vec& operator/=(const vec& other)
		{
			*this = *this / other;
			return *this;
		}
		//===================================    VEC X COM    =====================================
		vec  operator+(Com val) const
		{
			vec res(size);
			for (register int i = 0; i < size; i++)
				res.Vec[i] = Vec[i] + val;
			return res;
		}
		vec  operator-(Com val) const
		{
			vec res(size);
			for (register int i = 0; i < size; i++)
				res.Vec[i] = Vec[i] - val;
			return res;
		}
		vec  operator*(Com val) const
		{
			vec res(size);
			for (register int i = 0; i < size; i++)
				res.Vec[i] = Vec[i] * val;
			return res;
		}
		vec  operator/(Com val) const
		{
			zero_check(val);
			vec res(size);
			for (register int i = 0; i < size; i++)
				res.Vec[i] = Vec[i] / val;
			return res;
		}
		vec  operator^(Com val) const
		{
			vec res(size);
			for (register int i = 0; i < size; i++)
				res.Vec[i] = com::pow(Vec[i], val);
			return res;
		}
		array<bool> operator==(Com val) const
		{
			array<bool> res(size);
			for (register int i = 0; i < size; i++)
				res.push_back(Vec[i] == val);
			return res;
		}
		array<bool> operator!=(Com val) const
		{
			array<bool> res(size);
			for (register int i = 0; i < size; i++)
				res.push_back(Vec[i] != val);
			return res;
		}
		array<bool> operator> (Com val) const
		{
			array<bool> res(size);
			for (register int i = 0; i < size; i++)
				res.push_back(Vec[i] > val);
			return res;
		}
		array<bool> operator< (Com val) const
		{
			array<bool> res(size);
			for (register int i = 0; i < size; i++)
				res.push_back(Vec[i] < val);
			return res;
		}
		array<bool> operator>=(Com val) const
		{
			array<bool> res(size);
			for (register int i = 0; i < size; i++)
				res.push_back(Vec[i] >= val);
			return res;
		}
		array<bool> operator<=(Com val) const
		{
			array<bool> res(size);
			for (register int i = 0; i < size; i++)
				res.push_back(Vec[i] <= val);
			return res;
		}

		vec& operator+=(Com val)
		{
			*this = *this + val;
			return *this;
		}
		vec& operator-=(Com val)
		{
			*this = *this - val;
			return *this;
		}
		vec& operator*=(Com val)
		{
			*this = *this * val;
			return *this;
		}
		vec& operator/=(Com val)
		{
			*this = *this / val;
			return *this;
		}
		vec& operator^=(Com val)
		{
			*this = *this ^ val;
			return *this;
		}
		//===================================    VEC X NUM    =====================================
		vec  operator+(double val) const
		{
			vec res(size);
			for (register int i = 0; i < size; i++)
				res.Vec[i] = Vec[i] + val;
			return res;
		}
		vec  operator-(double val) const
		{
			vec res(size);
			for (register int i = 0; i < size; i++)
				res.Vec[i] = Vec[i] - val;
			return res;
		}
		vec  operator*(double val) const
		{
			vec res(size);
			for (register int i = 0; i < size; i++)
				res.Vec[i] = Vec[i] * val;
			return res;
		}
		vec  operator/(double val) const
		{
			zero_check(val);
			vec res(size);
			for (register int i = 0; i < size; i++)
				res.Vec[i] = Vec[i] / val;
			return res;
		}
		vec  operator^(double val) const
		{
			vec res(size);
			for (register int i = 0; i < size; i++)
				res.Vec[i] = com::pow(Vec[i], val);
			return res;
		}
		array<bool> operator==(double val) const
		{
			array<bool> res(size);
			for (register int i = 0; i < size; i++)
				res.push_back(Vec[i] == val);
			return res;
		}
		array<bool> operator!=(double val) const
		{
			array<bool> res(size);
			for (register int i = 0; i < size; i++)
				res.push_back(Vec[i] != val);
			return res;
		}
		array<bool> operator> (double val) const
		{
			array<bool> res(size);
			for (register int i = 0; i < size; i++)
				res.push_back(Vec[i] > val);
			return res;
		}
		array<bool> operator< (double val) const
		{
			array<bool> res(size);
			for (register int i = 0; i < size; i++)
				res.push_back(Vec[i] < val);
			return res;
		}
		array<bool> operator>=(double val) const
		{
			array<bool> res(size);
			for (register int i = 0; i < size; i++)
				res.push_back(Vec[i] >= val);
			return res;
		}
		array<bool> operator<=(double val) const
		{
			array<bool> res(size);
			for (register int i = 0; i < size; i++)
				res.push_back(Vec[i] <= val);
			return res;
		}

		vec& operator+=(double val)
		{
			*this = *this + val;
			return *this;
		}
		vec& operator-=(double val)
		{
			*this = *this - val;
			return *this;
		}
		vec& operator*=(double val)
		{
			*this = *this * val;
			return *this;
		}
		vec& operator/=(double val)
		{
			*this = *this / val;
			return *this;
		}
		vec& operator^=(double val)
		{
			*this = *this ^ val;
			return *this;
		}
		//=============================== << O P E R A T O R S \E.N.D>> ===================================
		void set(int i, Com val)
		{
			chk_rng(i, size);
			if (i < 0)
				error("Negative Index");
			Vec[i] = val;
		}
		void set(Com val)
		{
			for (register int i = 0; i < size; i++)
				Vec[i] = val;
		}
		void reset(int s)
		{
			if (Vec != nullptr) delete[] Vec;
			size = s;
			Vec = new Com[size] {};
		}
		void reset()
		{
			delete[] Vec;
			Vec = new Com[size] {};
		}
		void print() const
		{
			printf("\n(");
			for (register int i = 0; i < size - 1; i++)
			{
				Vec[i].print();
				printf(",");
			}
			Vec[size - 1].print();
			printf(")\n");
		}
		void erase() { Vec = nullptr; }
		~vec() { if (Vec != nullptr) delete[] Vec; }
	};
	vec operator*(Com val, const vec& v) { return v * val; }
	vec operator/(Com val, const vec& v)
	{
		int size = v.Size();
		vec res(size);
		for (register int i = 0; i < size; i++)
		{
			zero_check(v[i]);
			res.set(i, 1 / v[i]);
		}
		return res;
	}
	vec operator*(double val, const vec& v) { return v * val; }
	vec operator/(double val, const vec& v)
	{
		int size = v.Size();
		vec res(size);
		for (register int i = 0; i < size; i++)
		{
			zero_check(v[i]);
			res.set(i, 1 / v[i]);
		}
		return res;
	}
	//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
	Com sum(const vec& v)
	{
		Com sum = 0; int size = v.Size();
		for (register int i = 0; i < size; i++)
			sum += v[i];
		return sum;
	}
	Com Wsum(const vec& v1, const vec& v2)
	{
		chk_eq(v1.Size(), v2.Size());
		int size = v1.Size();
		Com sum = 0;
		for (register int i = 0; i < size; i++)
			sum += v1[i] * v2[i];
		return sum;
	}
	Com r_value(const vec& x, const vec& y)
	{
		chk_eq(x.Size(), y.Size());
		int n = x.Size();
		Com sumxy = Wsum(x, y), sumx = x.sum(), sumy = y.sum(), sumx2 = x.magsqr(), sumy2 = y.magsqr();
		return (n*sumxy - sumx * sumy) / com::sqrt((n*sumx2 - sumx * sumx)*(n*sumy2 - sumy * sumy));
	}
	Com mag(const vec& v)
	{
		Com sum = 0; int size = v.Size();
		for (register int i = 0; i < size; i++)
			sum += v[i] * v[i];
		return com::sqrt(sum);
	}
	Com min(const vec& v)
	{
		Com min = v[0]; int size = v.Size();
		for (register int i = 0; i < size; i++)
			if (v[i] < min)
				min = v[i];
		return min;
	}
	Com max(const vec& v)
	{
		Com max = v[0]; int size = v.Size();
		for (register int i = 0; i < size; i++)
			if (v[i] > max)
				max = v[i];
		return max;
	}
	int mini(const vec& v)
	{
		Com min = v[0]; int size = v.Size(), ind = 0;
		for (register int i = 0; i < size; i++)
			if (v[i] < min)
			{
				min = v[i];
				ind = i;
			}
		return ind;
	}
	int maxi(const vec& v)
	{
		Com max = v[0]; int size = v.Size(), ind = 0;
		for (register int i = 0; i < size; i++)
			if (v[i] > max)
			{
				max = v[i];
				ind = i;
			}
		return ind;
	}
	vec reverse(const vec& v)
	{
		int size = v.Size();
		vec res(size);
		for (register int i = 0; i < size; i++)
			res.set(i, v[size - i - 1]);
		return res;
	}
	vec hat(const vec& v) { return v / mag(v); }
	vec norm(const vec& v) { return v / max(v); }
	vec _set(const vec& v, Com(*F)(Com x))
	{
		int size = v.Size();
		vec res(size);
		for (register int i = 0; i < size; i++)
			res.set(i, (*F)(v[i]));
		return res;
	}

	vec abs(const vec& v) { return _set(v, com::abs); }
	vec ceil(const vec& v) { return _set(v, com::ceil); }
	vec floor(const vec& v) { return _set(v, com::floor); }
	
	vec exp(const vec& v) { return _set(v, com::exp); }
	vec pow(const vec& v, double ex)
	{
		int size = v.Size();
		vec res(size);
		for (register int i = 0; i < size; i++)
			res.set(i, com::pow(v[i], ex));
		return res;
	}
	vec ln(const vec& v) { return _set(v, com::ln); }
	vec log(const vec& v, double base = 10)
	{
		int size = v.Size();
		vec res(size);
		for (register int i = 0; i < size; i++)
			res.set(i, com::log(v[i]) / std::log(10));
		return res;
	}

	vec sin(const vec& v) { return _set(v, com::sin); }
	vec cos(const vec& v) { return _set(v, com::cos); }
	vec tan(const vec& v) { return _set(v, com::tan); }

	vec sec(const vec& v) { return _set(v, [](Com val)->Com {return 1 / com::cos(val); }); }
	vec csc(const vec& v) { return _set(v, [](Com val)->Com {return 1 / com::sin(val); }); }
	vec cot(const vec& v) { return _set(v, [](Com val)->Com {return 1 / com::tan(val); }); }

	vec arcsin(const vec& v) { return _set(v, com::arcsin); }
	vec arccos(const vec& v) { return _set(v, com::arccos); }
	vec arctan(const vec& v) { return _set(v, com::arctan); }

	vec arcsec(const vec& v) { return _set(v, [](Com val)->Com {return com::arccos(1 / val); }); }
	vec arccsc(const vec& v) { return _set(v, [](Com val)->Com {return com::arcsin(1 / val); }); }
	vec arccot(const vec& v) { return _set(v, [](Com val)->Com {return com::arctan(1 / val); }); }

	vec sinh(const vec& v) { return _set(v, com::sinh); }
	vec cosh(const vec& v) { return _set(v, com::cosh); }
	vec tanh(const vec& v) { return _set(v, com::tanh); }

	vec sech(const vec& v) { return _set(v, [](Com val)->Com {return com::cosh(val); }); }
	vec csch(const vec& v) { return _set(v, [](Com val)->Com {return com::sinh(val); }); }
	vec coth(const vec& v) { return _set(v, [](Com val)->Com {return com::tanh(val); }); }

	vec arcsinh(const vec& v) { return _set(v, arcsinh); }
	vec arccosh(const vec& v) { return _set(v, arccosh); }
	vec arctanh(const vec& v) { return _set(v, arctanh); }

	vec arcsech(const vec& v) { return _set(v, [](Com val)->Com {return com::arccosh(1 / val); }); }
	vec arccsch(const vec& v) { return _set(v, [](Com val)->Com {return com::arcsinh(1 / val); }); }
	vec arccoth(const vec& v) { return _set(v, [](Com val)->Com {return com::arctanh(1 / val); }); }

	std::ostream& operator<<(std::ostream& C, const vec& v) { v.print(); return C; }

	vec ones(int s) { return _set(vec(s), [](Com val)->Com { return 1; }); }
	vec rand(int s) { return _set(vec(s), [](Com val)->Com { return std::rand() % 10; }); }
	vec ar_seq(Com a, Com b, int n)
	{
		vec res(n);
		res.set(0, a);
		for (register int i = 1; i < n; i++)
			res.set(i, res[i - 1] + b);
		return res;
	}
	vec geo_seq(Com a, Com r, int n)
	{
		vec res(n);
		res.set(0, a);
		for (register int i = 1; i < n; i++)
			res.set(i, res[i - 1] * r);
		return res;
	}
	vec fib_seq(Com f1, Com f2, int n)
	{
		vec res(n);
		res.set(0, f1); res.set(1, f2);
		for (register int i = 2; i < n; i++)
			res.set(i, res[i - 1] + res[i - 2]);
		return res;
	}
	vec linspace(Com start, Com end, int n) { return ar_seq(start, (end - start) / n, n); }
	vec seq(Com(*f)(Com*), const vec& args, int n)
	{
		vec res(n);
		int num_of_args = args.Size();
		chk_rng(num_of_args, n, "Number of Initial Values Exceeds the Number of the Members of the Vector\nCall from seq");
		for (register int i = 0; i < num_of_args; i++)
			res.set(i, args[i]);
		for (register int i = num_of_args; i < n; i++)
			res.set(i, (*f)(res.begin() + i - num_of_args));
		return res;
	}
	vec seq(Com(*f)(Com), Com arg, int n)
	{
		vec res(n);
		res.set(0, arg);
		for (register int i = 1; i < n; i++)
			res.set(i, (*f)(res[i - 1]));
		return res;
	}
	vec seq(Com(*f)(Com), int n, Com start = 1, Com step = 1)
	{
		vec res(n);
		for (register int i = 0; i < n; i++)
			res.set(i, (*f)(start + i * step));
		return res;
	}
}

#undef fileName