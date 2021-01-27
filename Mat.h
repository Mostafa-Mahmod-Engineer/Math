#pragma once
#include"Vec.h"
#define fileName "Mat.h"

void construct(int& constructed, int& id, std::string Construction_Statement)
{
	return;
	constructed++;
	id = constructed;
	std::cout << Construction_Statement << "\n";
	printf("Constructed Matrices : %d\nMatrix id : %d\n", constructed, id);
}
void assign(int& constructed, int id, std::string Construction_Statement)
{
	return;
	constructed++;
	std::cout << "Construction Complete.  " << Construction_Statement << "\n";
	printf("Constructed Matrices : %d\nMatrix id : %d\n", constructed, id);
}
//===============================================================================
enum class matstate
{
	non, symmetric, hermitian, positive_definite, semi_definite, orthogonal, orthonormal,
	upper, lower, diagonal, bidiagonal_u, bidiagonal_l, tridiagonal,
	sparce, single
};
namespace real
{
	class mat
	{
		int row, col; double** Mat = nullptr;
		//int id; static int constructed;
	public:
		mat() {}
		explicit mat(int n)
	{
		//=================================================|
		//construct(constructed, id, "Square Matrix");=====|
		//=================================================|
		row = col = n;
		Mat = new double*[n];
		for (register int i = 0; i < n; i++)
			Mat[i] = new double[n] {};
	}
		explicit mat(int n, int m)
	{
		//======================================================|
		//construct(constructed, id, "Rectangular Matrix");=====|
		//======================================================|
		row = n; col = m;
		Mat = new double*[n];
		for (register int i = 0; i < n; i++)
			Mat[i] = new double[m] {};
	}
		mat(const real::vec& v)
	{
		//==========================================|
		//construct(constructed, id, "Vector");=====|
		//==========================================|
		row = v.Size(); col = 1;
		Mat = new double*[row];
		for (register int i = 0; i < row; i++)
			Mat[i] = new double(v[i]);
	}
		mat(int n, const real::vec& v)
	{
		//===================================================|
		//construct(constructed, id, "Vector Reshape");======|
		//===================================================|
		row = n; col = v.Size() / n;
		Mat = new double*[row];
		for (register int i = 0; i < row; i++)
			Mat[i] = new double[col];
		for (register int i = 0; i < row; i++)
			for (register int j = 0; j < col; j++)
				Mat[i][j] = v[i*col + j];
	}
		mat(const mat& other)
	{
		//====================================|
		//id = other.id;                  ====|
		//assign(constructed, id, "Copy");====|
		//====================================|
		row = other.row;
		col = other.col;
		Mat = new double*[row];
		for (register int i = 0; i < row; i++)
			Mat[i] = new double[col];
		for (register int i = 0; i < row; i++)
			for (register int j = 0; j < col; j++)
				Mat[i][j] = other.Mat[i][j];
	}
		mat(mat&& other)
	{
		//id = other.id;
		row = other.row;
		col = other.col;
		Mat = other.Mat;
		other.Mat = nullptr;
	}
		mat(const mat* other)
	{
		//=============================================|
		//id = other->id;                         =====|
		//assign(constructed, id, "Pointer Copy");=====|
		//=============================================|
		row = other->row;
		col = other->col;
		Mat = new double*[row];
		for (register int i = 0; i < row; i++)
			Mat[i] = new double[col];
		for (register int i = 0; i < row; i++)
			for (register int j = 0; j < col; j++)
				Mat[i][j] = other->Mat[i][j];
	}
	
		//=============================== << G E T T E R S >> ====================================
		int Row() const { return row; }
		int Col() const { return col; }
		int Size() const { return row * col; }
		int Rank() const
	{
		mat res(reduce());
		int range = std::fmin(row, col);
		for (register int i = 0; i < range; i++)
			if (res.Mat[i][i] == 0)
				return i;
	}
		bool is_symmetric() const
	{
		if (row != col)
			return false;
		for (register int i = 0; i < row; i++)
			for (register int j = 0; j < i; j++)
				if (Mat[i][j] != Mat[j][i])
					return false;
		return true;
	}
		bool is_positive() const
	{
		if (row != col)
			return false;
		mat res(reduce());
		for (register int i = 0; i < res.row; i++)
			if (res.Mat[i][i] <= 0)
				return false;
		return true;
	}
		bool is_semi_positive() const
	{
		if (row != col)
			return false;
		mat res(reduce());
		for (register int i = 0; i < res.row; i++)
			if (res.Mat[i][i] < 0)
				return false;
		return true;
	}
		bool is_triangular() const
	{
		if (Mat[0][1] == 0)
			for (register int i = 0; i < row; i++)
				for (register int j = i + 1; j < col; j++)
				{
					if (Mat[i][j] != 0)
						return false;
				}
		else if (Mat[1][0] == 0)
			for (register int i = 0; i < col; i++)
				for (register int j = i + 1; j < row; j++)
				{
					if (Mat[j][i] != 0)
						return false;
				}
		else
			return false;
		return true;
	}
		bool is_lowerTriangular() const
	{
		for (register int i = 0; i < row; i++)
			for (register int j = i + 1; j < col; j++)
				if (Mat[i][j] != 0)
					return false;
		return true;
	}
		bool is_upperTriangular() const
	{
		for (register int i = 0; i < col; i++)
			for (register int j = i + 1; j < row; j++)
				if (Mat[j][i] != 0)
					return false;
		return true;
	}
		bool is_square() const { return row == col; }
		double get(int i, int j) const { chk_rng(i, row); chk_rng(j, col); return Mat[i][j]; }
		//=============================== << G E T T E R S \E.N.D>> ===============================
	
		//=============================== << F U N C T I O N S >> =================================
		double tr() const
	{
		chk_eq(row, col, "Not a Square Matrix");
		double sum = 0;
		for (register int i = 0; i < row; i++)
			sum += Mat[i][i];
		return sum;
	}
		mat rsum(matstate state = matstate::non) const
	{
		mat res(row, 1); double sum;
		int dia = row * (row < col) + col * (col <= row);
		switch (state)
		{
		case matstate::upper:
			for (register int i = 0; i < row; i++)
			{
				sum = 0;
				for (register int j = i; j < col; j++)
					sum += Mat[i][j];
				res.Mat[i][0] = sum;
			}
			return res;
		case matstate::lower:
			for (register int i = 0; i < row; i++)
			{
				sum = 0;
				for (register int j = 0; j <= i; j++)
					sum += Mat[i][j];
				res.Mat[i][0] = sum;
			}
			return res;
		case matstate::diagonal:
			for (register int i = 0; i < dia; i++)
				res.Mat[i][0] = Mat[i][i];
			return res;
		default:
			for (register int i = 0; i < row; i++)
			{
				sum = 0;
				for (register int j = 0; j < col; j++)
					sum += Mat[i][j];
				res.Mat[i][0] = sum;
			}
			return res;
		}
	}
		mat csum(matstate state = matstate::non) const
	{
		mat res(1, col); double sum;
		int dia = row * (row < col) + col * (col <= row);
		switch (state)
		{
		case matstate::upper:
			for (register int i = 0; i < col; i++)
			{
				sum = 0;
				for (register int j = 0; j <= i; j++)
					sum += Mat[j][i];
				res.Mat[0][i] = sum;
			}
			return res;
		case matstate::lower:
			for (register int i = 0; i < col; i++)
			{
				sum = 0;
				for (register int j = i; j < row; j++)
					sum += Mat[j][i];
				res.Mat[0][i] = sum;
			}
			return res;
		case matstate::diagonal:
			for (register int i = 0; i < dia; i++)
				res.Mat[0][i] = Mat[i][i];
			return res;
		default:
			for (register int i = 0; i < col; i++)
			{
				sum = 0;
				for (register int j = 0; j < row; j++)
					sum += Mat[j][i];
				res.Mat[0][i] = sum;
			}
			return res;
		}
	}
		double sum(matstate state = matstate::non) const
	{
		double sum = 0;
		int dia = row * (row < col) + col * (col <= row);
		switch (state)
		{
		case matstate::upper:
			for (register int i = 0; i < row; i++)
				for (register int j = i; j < col; j++)
					sum += Mat[i][j];
			return sum;
		case matstate::lower:
			for (register int i = 0; i < row; i++)
				for (register int j = 0; j <= i; j++)
					sum += Mat[i][j];
			return sum;
		case matstate::diagonal:
			for (register int i = 0; i < dia; i++)
				sum += Mat[i][i];
			return sum;
		default:
			for (register int i = 0; i < row; i++)
				for (register int j = 0; j < col; j++)
					sum += Mat[i][j];
			return sum;
		}
	}
		mat rmax() const
	{
		mat res(row, 1);
		for (register int i = 0; i < row; i++)
		{
			double max = Mat[i][0];
			for (register int j = 0; j < col; j++)
				if (Mat[i][j] > max)
					max = Mat[i][j];
			res.Mat[i][0] = max;
		}
		return res;
	}
		mat cmax() const
	{
		mat res(1, col);
		for (register int j = 0; j < row; j++)
		{
			double max = Mat[0][j];
			for (register int i = 0; i < col; i++)
				if (Mat[i][j] > max)
					max = Mat[i][j];
			res.Mat[0][j] = max;
		}
		return res;
	}
		double max() const
	{
		double max = Mat[0][0];
		for (register int i = 0; i < row; i++)
			for (register int j = 1; j < col; j++)
				if (Mat[i][j] > max)
					max = Mat[i][j];
		return max;
	}
		mat rmin() const
	{
		mat res(row, 1);
		for (register int i = 0; i < row; i++)
		{
			double max = Mat[i][0];
			for (register int j = 0; j < col; j++)
				if (Mat[i][j] < max)
					max = Mat[i][j];
			res.Mat[i][0] = max;
		}
		return res;
	}
		mat cmin() const
	{
		mat res(1, col);
		for (register int j = 0; j < row; j++)
		{
			double max = Mat[0][j];
			for (register int i = 0; i < col; i++)
				if (Mat[i][j] < max)
					max = Mat[i][j];
			res.Mat[0][j] = max;
		}
		return res;
	}
		double min() const
	{
		double max = Mat[0][0];
		for (register int i = 0; i < row; i++)
			for (register int j = 1; j < col; j++)
				if (Mat[i][j] < max)
					max = Mat[i][j];
		return max;
	}
		double det(matstate state = matstate::non) const
	{
		chk_eq(row, col, "No Determinant Exists\nNot a Square Matrix");
		mat temp(this);
		double det = 1;
		switch (state)
		{
		case matstate::upper:
			for (register int i = 0; i < row; i++)
				det *= Mat[i][i];
			return det;
		case matstate::lower:
			return 1;
		case matstate::diagonal:
			for (register int i = 0; i < row; i++)
				det *= Mat[i][i];
			return det;
		case matstate::bidiagonal_u:
			for (register int i = 0; i < row; i++)
				det *= Mat[i][i];
			return det;
		case matstate::bidiagonal_l:
			for (register int i = 0; i < row; i++)
				det *= Mat[i][i];
			return det;
		case matstate::tridiagonal:
			for (register int i = 0; i < row - 1; i++)
			{
				temp._add(i + 1, i, -temp.Mat[i + 1][i] / temp.Mat[i][i]);
				det *= temp.Mat[i][i];
			}
			return det * temp.Mat[row - 1][row - 1];
		default:
			short sign = 1;
			for (register int i = 0; i < row; i++)
			{
				if (temp.Mat[i][i] == 0)
					for (register int j = i + 1; j < row; j++)
						if (temp.Mat[j][i] == 0 && j == row - 1)
							return 0;
						else if (temp.Mat[j][i] != 0)
						{
							temp.swap(j, i);
							sign *= -1;
							break;
						}
				for (register int j = i + 1; j < col; j++)
					if (temp.Mat[j][i] != 0)
						temp._add(j, i, -temp.Mat[j][i] / temp.Mat[i][i]);
				det *= temp.Mat[i][i];
			}
			return det * sign;
		}
	}
		mat inv(matstate state = matstate::non) const
	{
		if (row != col)
		{
			mat Inv(~*this);
			return row > col ? !(Inv**this)*Inv : Inv*!(*this*Inv);
		}
		mat res(mat::I(row));
		mat temp(this);
		switch (state)
		{
		case matstate::upper:
			for (register int i = 0; i < row; i++)
				if (res.Mat[i][i] == 0)
					error("Matrix is singular\nInvoked from inv()");
				else
					res._mult(i, 1 / res.Mat[i][i]);
			for (register int i = row - 1; i >= 0; i--)
				for (register int j = i - 1; j >= 0; j--)
					if (temp.Mat[j][i] != 0)
						res._add(j, i, -temp.Mat[j][i] / temp.Mat[i][i]);
			return res;

		case matstate::lower:
			for (register int i = 0; i < col; i++)
				for (register int j = i + 1; j < row; j++)
					if (temp.Mat[j][i] != 0)
						res._add(j, i, -temp.Mat[j][i]);
			return res;
		case matstate::diagonal:
			for (register int i = 0; i < row; i++)
				if (res.Mat[i][i] == 0)
					error("Matrix is singular\nInvoked from inv()");
				else
					res.Mat[i][i] = 1 / res.Mat[i][i];
			return res;

		case matstate::orthogonal:
			double mag;
			for (register int i = 0; i < col; i++)
			{
				mag = 0;
				for (register int j = 0; j < row; j++)
					mag += Mat[j][i] * Mat[j][i];
				mag = std::sqrt(mag);
				for (register int j = 0; j < row; j++)
					res.Mat[j][i] = Mat[j][i] / mag;
			}
			return ~res;

		case matstate::orthonormal:
			return ~*this;

		default:
			for (register int i = 0; i < col; i++)
			{
				if (temp.Mat[i][i] == 0)//_________________swap&zero-check_______________
				{																   ////||
					bool single = true;											   ////||
					for (register int j = i + 1; j < row; j++)					   ////||
						if (temp.Mat[j][i] != 0)								   ////||
						{														   ////||
							res.swap(j, i);									       ////||
							temp.swap(j, i);									   ////||
							single = false;										   ////||
							break;												   ////||
						}														   ////||
					if (single)													   ////||
						error("Matrix is singular\nInvoked from inv()");		   ////||
				}//________________________________________________________________////||
				if (temp.Mat[i][i] != 1)
				{
					res.mult(i, 1 / temp.Mat[i][i]);
					temp._mult(i, 1 / temp.Mat[i][i]);
				}
				for (register int j = i + 1; j < row; j++)
					if (j == i || temp.Mat[j][i] == 0)
						continue;
					else
					{
						res.add(j, i, -temp.Mat[j][i]);
						temp._add(j, i, -temp.Mat[j][i]);
					}
			}
			for (register int i = row - 1; i >= 0; i--)
				for (register int j = i - 1; j >= 0; j--)
					if (temp.Mat[j][i] != 0)
						res.add(j, i, -temp.Mat[j][i]);
			return res;
		}
	}
		//---------------------------------N O R M S--------------------------------|||
		double Norm(matstate state = matstate::non) const { return csum(state).max(); }
		double E_Norm() const
	{
		double sum = 0;
		for (register int i = 0; i < row; i++)
			for (register int j = 0; j < col; j++)
				sum += Mat[i][j] * Mat[i][j];
		return std::sqrt(sum);
	}
		double P_Norm(int p) const
	{
		double sum = 0;
		for (register int i = 0; i < row; i++)
			for (register int j = 0; j < col; j++)
				sum += std::pow(Mat[i][j], p);
		return std::pow(sum, 1 / p);
	}
		double inf_Norm() const
	{
		double* sum = new double[row];
		for (register int i = 0; i < row; i++)
			for (register int j = 0; j < col; j++)
				sum[i] += std::fabs(Mat[i][j]);
		double max = sum[0];
		for (register int i = 0; i < row; i++)
			if (sum[i] > max)
				max = sum[i];
		return max;
	}
		double cond(matstate state = matstate::non) const { return Norm(state) * inv(state).Norm(); }
		double E_cond() const { return E_Norm() * (!*this).E_Norm(); }
		double P_cond(int p) const { return P_Norm(p) * (!*this).P_Norm(p); }
		double inf_cond() const { return inf_Norm() * (!*this).inf_Norm(); }
		//--------------------------------------------------------------------------|||
		mat add(double val) const
	{
		mat res(this);
		for (register int i = 0; i < row; i++)
			for (register int j = 0; j < col; j++)
				res.Mat[i][j] += val;
		return res;
	}
		mat reduce() const
	{
		mat res(this);
		if (row > col)
			res._setrows(col, row, 0);
		int range = std::fmin(row, col);
		for (register int i = 0; i < range; i++)
		{
			bool skip = true;
			if (res.Mat[i][i] == 0)
			{
				for (register int j = i + 1; j < row; j++)
					if (res.Mat[j][i] != 0)
					{
						res.swap(i, j);
						res._mult(j, -1);
						break;
					}
					else if (res.Mat[j][i] == 0 && j == row - 1)
						skip = false;//to skip from the next loop
			}
			for (register int j = i + 1; j < range; j++)
				if (res.Mat[j][i] != 0 && skip)
					res._add(j, i, -res.Mat[j][i] / res.Mat[i][i]);
		} return res;
	}
		mat setrows(int firstrow, int lastrow, double val) const
	{
		chk_rng(lastrow, row);
		mat res(this);
		for (register int i = firstrow; i < lastrow; i++)
			for (register int j = 0; j < col; j++)
				res.Mat[i][j] = val;
		return res;
	}
		mat setrows(int firstrow, int lastrow, const real::vec& v) const
	{
		chk_rng(firstrow, row); chk_rng(lastrow, row);
		chk_eq((lastrow - firstrow)*col, v.Size(), "Vector Size doesn't Match the Rows' Size\nError Call from setrows()");
		mat res(this);
		for (register int i = firstrow; i < lastrow + 1; i++)
			for (register int j = 0; j < col; j++)
				res.Mat[i][j] = v[(i - firstrow)*col + j];
		return res;
	}
		mat setrow(int r, const real::vec& v) const
	{
		mat res(this);
		chk_eq(col, v.Size());
		for (register int i = 0; i < col; i++)
			res.Mat[r][i] = v[i];
		return res;
	}
		mat setcols(int firstcol, int lastcol, double val) const
	{
		chk_rng(lastcol, col);
		mat res(this);
		for (register int i = 0; i < col; i++)
			for (register int j = firstcol; j < lastcol; j++)
				res.Mat[i][j] = val;
		return res;
	}
		mat setcols(int firstcol, int lastcol, const real::vec& v) const
	{
		chk_rng(firstcol, col); chk_rng(lastcol, col);
		chk_eq((lastcol - firstcol)*row, v.Size(), "Vector Size doesn't Match the Columns' Size\nError Call from setcols()");
		mat res(this);
		for (register int i = 0; i < row; i++)
			for (register int j = firstcol; j < lastcol + 1; j++)
				res.Mat[i][j] = v[i*col + j - firstcol];
		return res;
	}
		mat setcol(int c, const real::vec& v) const
	{
		mat res(this);
		chk_eq(row, v.Size());
		for (register int i = 0; i < row; i++)
			res.Mat[i][c] = v[i];
		return res;
	}
		mat augment(const mat& m) const
	{
		chk_eq(m.row, row);
		mat res(row, col + m.col);
		for (register int i = 0; i < res.row; i++)
			for (register int j = 0; j < res.col; j++)
				j < col ? res.Mat[i][j] = Mat[i][j] : res.Mat[i][j] = m.Mat[i][j - col];
		return res;
	}
		mat append(const mat& m) const
	{
		chk_eq(m.col, col);
		mat res(row + m.row, col);
		for (register int i = 0; i < row; i++)
			for (register int j = 0; j < col; j++)
				res.Mat[i][j] = Mat[i][j];
		for (register int i = row; i < res.row; i++)
			for (register int j = 0; j < col; j++)
				res.Mat[i][j] = m.Mat[i - row][j];
		return res;
	}
		mat augment(const real::vec& v) const
	{
		chk_eq(row, v.Size());
		mat res(row, col + 1);
		for (register int i = 0; i < row; i++)
		{
			for (register int j = 0; j < col; j++)
				res.Mat[i][j] = Mat[i][j];
			res.Mat[i][col] = v[i];
		}
		return res;
	}
		mat append(const real::vec& v) const
	{
		chk_eq(col, v.Size());
		mat res(row + 1, col);
		for (register int i = 0; i < row; i++)
		{
			for (register int j = 0; j < col; j++)
				res.Mat[i][j] = Mat[i][j];
			res.Mat[row][i] = v[i];
		}
		return res;
	}
		mat diagonal() const
	{
		mat res(row, col);
		int range = std::fmin(row, col);
		for (register int i = 0; i < range; i++)
			res.Mat[i][i] = Mat[i][i];
		return res;
	}
		mat off_diagonal() const
	{
		mat res(this);
		int range = std::fmin(row, col);
		for (register int i = 0; i < range; i++)
			res.Mat[i][i] = 0;
		return res;
	}
		mat submat(int firstrow, int lastrow, int firstcol, int lastcol) const
	{
		chk_rng(lastrow, row); chk_rng(lastcol, col);
		mat res(lastrow - firstrow + 1, lastcol - firstcol + 1);
		for (register int i = firstrow; i <= lastrow; i++)
			for (register int j = firstcol; j <= lastcol; j++)
				res.Mat[i - firstrow][j - firstcol] = Mat[i][j];
		return res;
	}
		mat submat(int lastrow, int lastcol) const
	{
		chk_rng(lastrow, row); chk_rng(lastcol, col);
		mat res(lastrow + 1, lastcol + 1);
		for (register int i = 0; i <= lastrow; i++)
			for (register int j = 0; j <= lastcol; j++)
				res.Mat[i][j] = Mat[i][j];
		return res;
	}
		mat submat(int pivot) const
	{
		chk_rng(pivot, row); chk_rng(pivot, col);
		mat res(pivot + 1, pivot + 1);
		for (register int i = 0; i <= pivot; i++)
			for (register int j = 0; j <= pivot; j++)
				res.Mat[i][j] = Mat[i][j];
		return res;
	}real::vec unroll() const
	{
		real::vec res(row*col);
		for (register int i = 0; i < row; i++)
			for (register int j = 0; j < col; j++)
				res.set(i*col + j, Mat[i][j]);
		return res;
	}
		real::vec getrow(int n) const
	{
		chk_rng(n, row);
		real::vec res(col);
		for (register int i = 0; i < col; i++)
			res.set(i, Mat[n][i]);
		return res;
	}
		real::vec getcol(int n) const
	{
		chk_rng(n, col);
		real::vec res(row);
		for (register int i = 0; i < row; i++)
			res.set(i, Mat[i][n]);
		return res;
	}
		
		void _add(double val)
	{
		for (register int i = 0; i < row; i++)
			for (register int j = 0; j < col; j++)
				Mat[i][j] += val;
	}
		void _reduce()
	{
		if (row > col)
			_setrows(col - 1, row - 1, 0);
		int range = col * (row > col) + (row - 1)*(row <= col);
		bool sign = false;
		for (register int i = 0; i < range; i++)
		{
			if (Mat[i][i] == 0)
			{
				bool skip = false;
				for (register int k = i; k < col; k++)
				{
					for (register int j = i + 1; j < row; j++)
						if (Mat[j][k] != 0)
						{
							swap(i, j);
							skip = true;
							sign = !sign;
							break;
						}
					if (skip)
						break;
				}
			}
			for (register int j = i + 1; j < range; j++)
				if (Mat[j][i] != 0)
					_add(j, i, -Mat[j][i] / Mat[i][i]);
		}
		if (sign)
			_mult(row - 1, -1);
	}
		void _diagonal()
	{
		_reduce();
		int range = std::fmin(row, col);
		for (register int i = 0; i < range - 1; i++)
			for (register int j = i + 1; j < col; j++)
				Mat[i][j] = 0;
	}
		void _setrows(int firstrow, int lastrow, double val)
	{
		chk_rng(lastrow, row);
		for (register int i = firstrow; i <= lastrow; i++)
			for (register int j = 0; j < col; j++)
				Mat[i][j] = val;
	}
		void _setrows(int firstrow, int lastrow, const real::vec& v)
	{
		chk_rng(firstrow, row); chk_rng(lastrow, row);
		chk_eq((lastrow - firstrow)*col, v.Size(), "Vector Size doesn't Match the Rows' Size\nError Call from setrows()");
		for (register int i = firstrow; i < lastrow + 1; i++)
			for (register int j = 0; j < col; j++)
				Mat[i][j] = v[(i - firstrow)*col + j];
	}
		void _setrow(int r, const real::vec& v)
	{
		chk_eq(col, v.Size()); chk_rng(r, row);
		for (register int i = 0; i < col; i++)
			Mat[r][i] = v[i];
	}
		void _setcols(int firstcol, int lastcol, double val)
	{
		chk_rng(lastcol, col);
		for (register int i = 0; i <= col; i++)
			for (register int j = firstcol; j < lastcol; j++)
				Mat[i][j] = val;
	}
		void _setcols(int firstcol, int lastcol, const real::vec& v)
	{
		chk_rng(firstcol, col); chk_rng(lastcol, col);
		chk_eq((lastcol - firstcol)*row, v.Size(), "Vector Size doesn't Match the Columns' Size\nError Call from setcols()");
		for (register int i = 0; i < row; i++)
			for (register int j = firstcol; j < lastcol + 1; j++)
				Mat[i][j] = v[i*col + j - firstcol];
	}
		void _setcol(int c, const real::vec& v)
	{
		chk_eq(row, v.Size()); chk_rng(c, col);
		for (register int i = 0; i < row; i++)
			Mat[i][c] = v[i];
	}
		void _swapcols(int first, int second)
	{
		chk_rng(first, col); chk_rng(second, col);
		double temp;
		for (register int i = 0; i < row; i++)
		{
			temp = Mat[i][first];
			Mat[i][first] = Mat[i][second];
			Mat[i][second] = temp;
		}
	}
		void _augment(const mat& m)
	{
		chk_eq(row, m.row);
		mat res(row, col + m.col);
		for (register int i = 0; i < row; i++)
		{
			for (register int j = 0; j < col; j++)
				res.Mat[i][j] = Mat[i][j];
			for (register int j = 0; j < m.col; j++)
				res.Mat[i][col + j] = m.Mat[i][j];
		}
		col = res.col;
		Mat = res.Mat;
		res.Mat = nullptr;
	}
		void _augment(const real::vec& v)
	{
		chk_eq(row, v.Size());
		mat res(row, col + 1);
		for (register int i = 0; i < row; i++)
		{
			for (register int j = 0; j < col; j++)
				res.Mat[i][j] = Mat[i][j];
			res.Mat[i][col] = v[i];
		}
		col++;
		Mat = res.Mat;
		res.Mat = nullptr;
	}
		void _append(const mat& m)
	{
		chk_eq(col, m.col);
		mat res(row + m.row, col);
		for (register int i = 0; i < col; i++)
		{
			for (register int j = 0; j < row; j++)
				res.Mat[j][i] = Mat[j][i];
			for (register int j = 0; j < m.row; j++)
				res.Mat[row + j][i] = m.Mat[j][i];
		}
		row = res.row;
		Mat = res.Mat;
		res.Mat = nullptr;
	}
		void _append(const real::vec& v)
	{
		chk_eq(col, v.Size());
		mat res(row + 1, col);
		for (register int i = 0; i < col; i++)
		{
			for (register int j = 0; j < row; j++)
				res.Mat[j][i] = Mat[j][i];
			res.Mat[row][i] = v[i];
		}
		row++;
		Mat = res.Mat;
		res.Mat = nullptr;
	}
		void _submat(int firstrow, int lastrow, int firstcol, int lastcol)
	{
		chk_rng(lastrow, row); chk_rng(lastcol, col);
		mat res(lastrow - firstrow, lastcol - firstcol);
		for (register int i = firstrow; i < lastrow; i++)
			for (register int j = firstcol; j < lastcol; j++)
				res.Mat[i - firstrow][j - firstcol] = Mat[i][j];
		row = res.row;
		col = res.col;
		Mat = res.Mat;
		res.Mat = nullptr;
	}
		void _make_symmetric()
	{
		chk_eq(row, col, "Matrix is not Square\nCan't be Made Symmetric");
		for (register int i = 0; i < row; i++)
			for (register int j = i + 1; j < col; j++)
				Mat[j][i] = Mat[i][j];
	}
		//=============================== << F A C T O R I Z A T I O N >> =================================
		list<mat> lu() const
	{
		mat l(mat::I(row)), u(this), p(l);
		int dia = col * (row > col) + (row - 1)*(row >= col);
		for (register int i = 0; i < dia; i++)
		{
			if (u.Mat[i][i] == 0) // if a pivot is zero:
			{
				bool flag = false;
				for (register int j = i; j < col; j++)
				{
					for (register int k = i + 1; k < row; k++)
						if (u.Mat[k][j] != 0)
						{
							u.swap(i, k);
							p.swap(i, k);
							l._swapcols(i, k);
							flag = true;
							break;
						}
					if (flag)
						break;
				}
			}
			for (register int j = i + 1; j < row; j++)
			{
				if (u.Mat[j][i] == 0)
					continue;
				l.Mat[j][i] = u.Mat[j][i] / u.Mat[i][i];
				u._add(j, i, -u.Mat[j][i] / u.Mat[i][i]);
			}
		}
		list<mat> res;
		res.add(l);
		res.add(u);
		res.add(~p);
		return res;
	}
		void lu(mat& L, mat& U, mat& P) const
	{
		mat l(mat::I(row)), u(this), p(l);
		int dia = col * (row > col) + (row - 1)*(row >= col);
		for (register int i = 0; i < dia; i++)
		{
			if (u.Mat[i][i] == 0) // if a pivot is zero:
			{
				bool flag = false;
				for (register int j = i; j < col; j++)
				{
					for (register int k = i + 1; k < row; k++)
						if (u.Mat[k][j] != 0)
						{
							u.swap(i, k);
							p.swap(i, k);
							l._swapcols(i, k);
							flag = true;
							break;
						}
					if (flag)
						break;
				}
			}
			for (register int j = i + 1; j < row; j++)
			{
				if (u.Mat[j][i] == 0)
					continue;
				l.Mat[j][i] = u.Mat[j][i] / u.Mat[i][i];
				u._add(j, i, -u.Mat[j][i] / u.Mat[i][i]);
			}
		}
		move(l, L);
		move(u, U);
		P = ~p;
	}
		void lu(mat& L, mat& U) const
	{
		mat l(mat::I(row)), u(this), p(l);
		int dia = col * (row > col) + (row - 1)*(row >= col);
		for (register int i = 0; i < dia; i++)
		{
			if (u.Mat[i][i] == 0) // if a pivot is zero:
			{
				bool flag = false;
				for (register int j = i; j < col; j++)
				{
					for (register int k = i + 1; k < row; k++)
						if (u.Mat[k][j] != 0)
						{
							u.swap(i, k);
							p.swap(i, k);
							l._swapcols(i, k);
							flag = true;
							break;
						}
					if (flag)
						break;
				}
			}
			for (register int j = i + 1; j < row; j++)
			{
				if (u.Mat[j][i] == 0)
					continue;
				l.Mat[j][i] = u.Mat[j][i] / u.Mat[i][i];
				u._add(j, i, -u.Mat[j][i] / u.Mat[i][i]);
			}
		}
		move(l, L);
		move(u, U);
	}
	
		list<mat> ldu() const
	{
		mat l(mat::I(row)), d(row), u(this), p(l);
		int dia = col * (row > col) + (row - 1)*(row >= col);
		for (register int i = 0; i < dia; i++)
		{
			if (u.Mat[i][i] == 0) // if a pivot is zero:
			{
				bool flag = false;
				for (register int j = i; j < col; j++)
				{
					for (register int k = i + 1; k < row; k++)
						if (u.Mat[k][j] != 0)
						{
							u.swap(i, k);
							p.swap(i, k);
							l._swapcols(i, k);
							flag = true;
							break;
						}
					if (flag)
						break;
				}
			}
			for (register int j = i + 1; j < row; j++)
			{
				if (u.Mat[j][i] == 0)
					continue;
				l.Mat[j][i] = u.Mat[j][i] / u.Mat[i][i];
				u._add(j, i, -u.Mat[j][i] / u.Mat[i][i]);
			}
		}
		dia += row == col;
		for (register int i = 0; i < dia; i++)
		{
			d.Mat[i][i] = u.Mat[i][i];
			u._mult(i, 1 / u.Mat[i][i]);
		}
		list<mat> res;
		res.add(l);
		res.add(d);
		res.add(u);
		res.add(~p);
		return res;
	}
		void ldu(mat& L, mat& D, mat& U, mat& P) const
	{
		mat l(mat::I(row)), d(row), u(this), p(l);
		int dia = col * (row > col) + (row - 1)*(row >= col);
		for (register int i = 0; i < dia; i++)
		{
			if (u.Mat[i][i] == 0) // if a pivot is zero:
			{
				bool flag = false;
				for (register int j = i; j < col; j++)
				{
					for (register int k = i + 1; k < row; k++)
						if (u.Mat[k][j] != 0)
						{
							u.swap(i, k);
							p.swap(i, k);
							l._swapcols(i, k);
							flag = true;
							break;
						}
					if (flag)
						break;
				}
			}
			for (register int j = i + 1; j < row; j++)
			{
				if (u.Mat[j][i] == 0)
					continue;
				l.Mat[j][i] = u.Mat[j][i] / u.Mat[i][i];
				u._add(j, i, -u.Mat[j][i] / u.Mat[i][i]);
			}
		}
		dia += row == col;
		for (register int i = 0; i < dia; i++)
		{
			d.Mat[i][i] = u.Mat[i][i];
			u._mult(i, 1 / u.Mat[i][i]);
		}
		move(l, L);
		move(d, D);
		move(u, U);
		p = ~p;
	}
		void ldu(mat& L, mat& D, mat& U) const
	{
		mat l(mat::I(row)), d(row), u(this), p(l);
		int dia = col * (row > col) + (row - 1)*(row <= col);
		for (register int i = 0; i < dia; i++)
		{
			if (u.Mat[i][i] == 0) // if a pivot is zero:
			{
				bool flag = false;
				for (register int j = i; j < col; j++)
				{
					for (register int k = i + 1; k < row; k++)
						if (u.Mat[k][j] != 0)
						{
							u.swap(i, k);
							p.swap(i, k);
							l._swapcols(i, k);
							flag = true;
							break;
						}
					if (flag)
						break;
				}
			}
			for (register int j = i + 1; j < row; j++)
			{
				if (u.Mat[j][i] == 0)
					continue;
				l.Mat[j][i] = u.Mat[j][i] / u.Mat[i][i];
				u._add(j, i, -u.Mat[j][i] / u.Mat[i][i]);
			}
		}
		dia += row == col;
		for (register int i = 0; i < dia; i++)
		{
			d.Mat[i][i] = u.Mat[i][i];
			u._mult(i, 1 / u.Mat[i][i]);
		}
		L = (~p)*l;
		move(d, D);
		move(u, U);
	}
	
		pair<mat, mat> qr() const
	{
		mat q(row, col), r(col);
		real::vec* a; a = new real::vec[col];
		for (register int i = 0; i < col; i++)
			a[i] = getcol(i);
		real::vec sum(row);
		for (register int i = 0; i < col; i++)
		{
			for (register int j = 0; j < i; j++)
			{
				r.Mat[j][i] = a[i].Wsum(a[j]);
				sum += r.Mat[j][i] * a[j];
			}
			a[i] -= sum;
			r.Mat[i][i] = a[i].mag();
			a[i] /= r.Mat[i][i];
			q._setcol(i, a[i]);
			sum.reset();
		}
		return pair<mat, mat>(q, r);
	}
		void qr(mat& Q, mat& R) const
	{
		mat q(row, col), r(col);
		real::vec* a; a = new real::vec[col];
		for (register int i = 0; i < col; i++)
			a[i] = getcol(i);
		real::vec sum(row);
		for (register int i = 0; i < col; i++)
		{
			for (register int j = 0; j < i; j++)
			{
				r.Mat[j][i] = a[i].Wsum(a[j]);
				sum += r.Mat[j][i] * a[j];
			}
			a[i] -= sum;
			r.Mat[i][i] = a[i].mag();
			a[i] /= r.Mat[i][i];
			q._setcol(i, a[i]);
			sum.reset();
		}
		Q.row = q.row; Q.col = q.col; Q.Mat = q.Mat; q.Mat = nullptr;
		R.row = r.row; R.col = r.col; R.Mat = r.Mat; r.Mat = nullptr;
	}
	
		//    NOT FINISHED YET
		pair<mat, mat> eig(matstate state = matstate::non) const
	{
		if (row != col)
			error("Matrix is Rectangular\nError Call from eig()");


		mat temp(reduce()), X(row), lambda(row);
		real::vec eig(row);
		int zeroes;
		for (zeroes = row - 1; zeroes >= 0; zeroes--)
			if (temp.Mat[zeroes][zeroes] > 1.0e-7 || temp.Mat[zeroes][zeroes] < -1.0e-7)
				break;
		zeroes = row - zeroes - 1;
		if (zeroes > 0) // zero eigenvalue found:
		{
			for (register int i = row - zeroes - 1; i > 0; i--)
				for (register int j = i - 1; j >= 0; j--)
					temp._add(j, i, -temp.Mat[j][i] / temp.Mat[i][i]);
			real::vec x(row);
			for (register int i = col - 1; i >= col - zeroes; i--)
			{
				for (register int j = 0; j < row - zeroes; j++)
					x.set(j, -temp.Mat[j][i] / temp.Mat[j][j]);
				x.set(i, 1);
				X._setcol(i, x);
				x.set(i, 0);
			}
		}
		//==================================================
		auto power = [](const mat& A, real::vec& x)->double
		{
			mat res = A * A * A;
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
			return x0;
		};
		auto Power = [](const mat& A, real::vec& x)->double
		{
			double x0 = 0, x1, e = 1;
			while (e > 5.0e-7)
			{
				x1 = x0;
				x = A * x;
				x0 = x.abs_max();
				x._norm();
				e = x0 - x1;
			}
			return x0;
		};
		//==================================================
		if (state == matstate::symmetric || is_symmetric())
		{
			auto deflate = [](mat& A, real::vec& X, double eig)
			{
				mat x(X.hat());
				A = A - x * ~x*eig;
			};
			temp = *this; real::vec x(row);
			//___The Main Loop______________________
			for (register int i = 0; i < row - zeroes; i++)
			{
				x.set(1);
				eig.set(i, power(temp, x));
				lambda.Mat[i][i] = eig[i];
				X._setcol(i, x.hat());
				deflate(temp, x, eig[i]);
			}//|||||||||||||||||||||||||||||||||||||
			return pair<mat, mat>(X, lambda);
		}
		else
		{
			double max, min;
			real::vec x = real::vec::O(row);
			eig.set(0, power(temp, x));
			lambda.Mat[0][0] = max = eig[0];
			X._setcol(0, x.hat());
			if (zeroes > 0)
			{
				min = 0;
				zeroes--;
			}
			else
			{
				x = real::vec::O(row);
				eig.set(0, 1 / power(!temp, x));
				lambda.Mat[0][0] = min = eig[0];
				X._setcol(0, x.hat());
			}
			double alpha = (max + min) / 2;
			for (register int i = 1; i <= row - zeroes - 2; i++)
			{
				x = real::vec::O(row);
				eig.set(i, 1 / power(!(temp - alpha) + alpha, x));
				if (eig[i] != max && eig[i] != min)
				{
					lambda.Mat[i][i] = eig[0];
					X._setcol(i, x.hat());
				}
				matstate state = matstate::non;
			}
		}
	}
		void eig(mat& x, mat& e) const;
	
		list<mat> svd() const
	{
		mat res1(*this * ~*this), res2(~*this * *this);

	}
		void svd(mat& u, mat& s, mat& v) const
	{

	}
		//======================================== << S T A T I C >> ========================================
		static mat I(int n)
	{
		mat res(n);
		for (register int i = 0; i < n; i++)
			res.Mat[i][i] = 1;
		return res;
	}
		static mat O(int n)
	{
		mat res(n);
		for (register int i = 0; i < n; i++)
			for (register int j = 0; j < n; j++)
				res.Mat[i][j] = 1;
		return res;
	}
		static mat O(int n, int m)
	{
		mat res(n, m);
		for (register int i = 0; i < n; i++)
			for (register int j = 0; j < m; j++)
				res.Mat[i][j] = 1;
		return res;
	}
		static mat R(int n)
	{
		mat res(n);
		for (register int i = 0; i < n; i++)
			for (register int j = 0; j < n; j++)
				res.Mat[i][j] = (std::rand() % 20) - 10;
		return res;
	}
		static mat R(int n, int m)
	{
		mat res(n, m);
		for (register int i = 0; i < n; i++)
			for (register int j = 0; j < m; j++)
				res.Mat[i][j] = (std::rand() % 20) - 10;
		return res;
	}
		static mat R(int n, double upperlimit)
	{
		mat res(n);
		for (register int i = 0; i < n; i++)
			for (register int j = 0; j < n; j++)
				res.Mat[i][j] = std::fmodf(std::rand(), upperlimit);
		return res;
	}
		static mat R(int n, int m, double upperlimit)
	{
		mat res(n, m);
		for (register int i = 0; i < n; i++)
			for (register int j = 0; j < m; j++)
				res.Mat[i][j] = std::fmodf(std::rand(), upperlimit);
		return res;
	}
		static mat R(int n, double lowerlimit, double upperlimit)
	{
		mat res(n);
		for (register int i = 0; i < n; i++)
			for (register int j = 0; j < n; j++)
				res.Mat[i][j] = std::fmodf(std::rand(), upperlimit - lowerlimit) + lowerlimit;
		return res;
	}
		static mat R(int n, int m, double lowerlimit, double upperlimit)
	{
		mat res(n, m);
		for (register int i = 0; i < n; i++)
			for (register int j = 0; j < m; j++)
				res.Mat[i][j] = std::fmodf(std::rand(), upperlimit - lowerlimit) + lowerlimit;
		return res;
	}
		static mat upper(int n)
	{
		mat res(n);
		for (register int i = 0; i < n; i++)
			for (register int j = i; j < n; j++)
				res.Mat[i][j] = std::rand() % 10;
		return res;
	}
		static mat lower(int n)
	{
		mat res(n);
		for (register int i = 0; i < n; i++)
		{
			for (register int j = 0; j < i; j++)
				res.Mat[i][j] = std::rand() % 10;
			res.Mat[i][i] = 1;
		}
		return res;
	}
		static mat diagonal(int n)
	{
		mat res(n);
		for (register int i = 0; i < n; i++)
			res.Mat[i][i] = std::rand() % 10;
		return res;
	}
		static mat symmetric(int n)
	{
		mat res(n);
		int size = (n + n * n) / 2;
		register int k = 0;
		for (register int i = 0; i < n; i++)
			for (register int j = i; j < n; j++, k++)
				res.Mat[i][j] = res.Mat[j][i] = std::rand() % 10;
		return res;
	}
		static mat grid(int n)
	{
		mat res(n);
		for (register int i = 1; i < n; i++)
			for (register int j = 0; j < n; j++)
				res.Mat[i][j] = i;
		return res;
	}
		static mat magic(int n)
	{
		mat res(n);
		if (n % 2 == 1)
		{
			int i = 1, j = n / 2 + 1, iprev, jprev;
			for (register int k = 1; k <= n * n; k++)
			{
				iprev = i;
				jprev = j;
				res.Mat[i - 1][j - 1] = k;
				i = n - (n - i + 1) % n;
				j = j % n + 1;
				if (res.Mat[i - 1][j - 1] != 0)
				{
					i = iprev % n + 1;
					j = jprev;
				}
			}
		}
		else if (n % 4 == 0)
		{
			printf("doubly even magic square is not implemented yet"); system("pause");
		}
		else
		{
			printf("singly even magic square is not implemented yet"); system("pause");
		}
		return res;
	}
		//=============================== << F U N C T I O N S \E.N.D>> ===========================
	
		//=============================== << O P E R A T O R S >> =================================
		mat operator-() const
	{
		mat res(row, col);
		for (register int i = 0; i < row; i++)
			for (register int j = 0; j < col; j++)
				res.Mat[i][j] = -Mat[i][j];
		return res;
	}
		mat operator~() const
	{
		mat res(col, row);
		for (register int i = 0; i < row; i++)
			for (register int j = 0; j < col; j++)
				res.Mat[j][i] = Mat[i][j];
		return res;
	}
		mat operator!() const
	{
		if (row != col)
		{
			mat inv(~*this);
			return row > col ? !(inv**this)*inv : inv*!(*this*inv);
		}
		mat res(mat::I(row));
		mat temp(this);
		for (register int i = 0; i < col; i++)
		{
			if (temp.Mat[i][i] == 0)//_________________swap&zero-check_______________
			{																   ////||
				bool single = true;											   ////||
				for (register int j = i + 1; j < row; j++)					   ////||
					if (temp.Mat[j][i] != 0)								   ////||
					{														   ////||
						res.swap(j, i);									       ////||
						temp.swap(j, i);									   ////||
						single = false;										   ////||
						break;												   ////||
					}														   ////||
				if (single)													   ////||
					error("Matrix is singular\nInvoked from !-matrix");		   ////||
			}//________________________________________________________________////||
			if (temp.Mat[i][i] != 1)
			{
				res.mult(i, 1 / temp.Mat[i][i]);
				temp._mult(i, 1 / temp.Mat[i][i]);
			}
			for (register int j = i + 1; j < row; j++)
				if (j == i || temp.Mat[j][i] == 0)
					continue;
				else
				{
					res.add(j, i, -temp.Mat[j][i]);
					temp._add(j, i, -temp.Mat[j][i]);
				}
		}
		for (register int i = row - 1; i >= 0; i--)
			for (register int j = i - 1; j >= 0; j--)
				if (temp.Mat[j][i] != 0)
					res.add(j, i, -temp.Mat[j][i]);
		return res;
	}
		//===================================    MAT X MAT    =====================================
		mat operator+(const mat& other) const
	{
		chk_eq(other.row, row); chk_eq(other.col, col);
		mat res(row, col);
		for (register int i = 0; i < row; i++)
			for (register int j = 0; j < col; j++)
				res.Mat[i][j] = Mat[i][j] + other.Mat[i][j];
		return res;
	}
		mat operator-(const mat& other) const
	{
		chk_eq(other.row, row); chk_eq(other.col, col);
		mat res(row, col);
		for (register int i = 0; i < row; i++)
			for (register int j = 0; j < col; j++)
				res.Mat[i][j] = Mat[i][j] - other.Mat[i][j];
		return res;
	}
		mat operator*(const mat& other) const
	{
		chk_eq(col, other.row);
		mat res(row, other.col);
		double sum;
		for (register int i = 0; i < res.row; i++)
			for (register int j = 0; j < res.col; j++)
			{
				sum = 0;
				for (register int k = 0; k < col; k++)
					sum += Mat[i][k] * other.Mat[k][j];
				res.Mat[i][j] = sum;
			}
		return res;
	}
		mat operator/(const mat& other) const { return *this * !other; }
		mat operator|(const mat& other) const { return augment(other); }
		bool operator==(const mat& other) const
	{
		chk_eq(row, other.row); chk_eq(col, other.col);
		for (register int i = 0; i < row; i++)
			for (register int j = 0; j < col; j++)
				if (Mat[i][j] != other.Mat[i][j])
					return false;
		return true;
	}
		bool operator!=(const mat& other) const { return !(*this == other); }
	
		mat& operator =(const mat& other)
	{
		//===========================================|
		//id = other.id;                        =====|
		//assign(constructed, id, "Assignment");=====|
		//===========================================|
		if (Mat == other.Mat)
			return *this;
		this->~mat();
		row = other.row;
		col = other.col;
		Mat = new double*[row];
		for (register int i = 0; i < row; i++)
			Mat[i] = new double[col];
		for (register int i = 0; i < row; i++)
			for (register int j = 0; j < col; j++)
				Mat[i][j] = other.Mat[i][j];
		return *this;
	}
		mat& operator =(mat&& other)
	{
		//id = other.id;
		row = other.row;
		col = other.col;
		Mat = other.Mat;
		other.Mat = nullptr;
		return *this;
	}
		mat& operator+=(const mat& other) { *this = *this + other; return *this; }
		mat& operator-=(const mat& other) { *this = *this - other; return *this; }
		mat& operator*=(const mat& other) { *this = *this * other; return *this; }
		mat& operator/=(const mat& other) { *this = *this / other; return *this; }
		mat& operator|=(const mat& other) { *this = *this | other; return *this; }
		//===================================    MAT X VEC    =====================================
		real::vec operator*(const real::vec& v) const
	{
		chk_eq(col, v.Size());
		real::vec res(row);
		double sum;
		for (register int i = 0; i < row; i++)
		{
			sum = 0;
			for (register int j = 0; j < col; j++)
				sum += Mat[i][j] * v[j];
			res.set(i, sum);
		}
		return res;
	}
		real::vec& operator*=(const real::vec& v)
	{
		*this = *this * v;
		real::vec res(this->unroll());
		return res;
	}
	
		mat& operator =(const real::vec& other)
	{
		//=====================================================|
		//construct(constructed, id, "Vector Assignment");=====|
		//=====================================================|
		row = other.Size();
		col = 1;
		Mat = new double*[row];
		for (register int i = 0; i < row; i++)
			Mat[i] = new double{ other[i] };
		return *this;
	}
		//===================================    MAT X NUM    =====================================
		mat operator*(double val) const
	{
		mat res(row, col);
		for (register int i = 0; i < row; i++)
			for (register int j = 0; j < col; j++)
				res.Mat[i][j] = Mat[i][j] * val;
		return res;
	}
		mat operator/(double val) const
	{
		zero_check(val);
		mat res(row, col);
		for (register int i = 0; i < row; i++)
			for (register int j = 0; j < col; j++)
				res.Mat[i][j] = Mat[i][j] / val;
		return res;
	}
		mat operator%(double val) const
	{
		zero_check(val, "Zero Modulus");
		mat res(row, col);
		for (register int i = 0; i < row; i++)
			for (register int j = 0; j < col; j++)
				res.Mat[i][j] = std::fmodl(Mat[i][j], val);
		return res;
	}
		mat operator+(double val) const { return (*this) + mat::I(std::fmin(row, col))*val; }
		mat operator-(double val) const { return (*this) - mat::I(std::fmin(row, col))*val; }

	
		mat& operator*=(double val) { *this = *this * val; return *this; }
		mat& operator/=(double val) { *this = *this / val; return *this; }
		mat& operator%=(double val) { *this = *this % val; return *this; }
		mat& operator+=(double val) { *this = *this + val; return *this; }
		mat& operator-=(double val) { *this = *this - val; return *this; }
		//=============================== << O P E R A T O R S \E.N.D>> ===================================
	
		//============================= << R O W   O P E R A T I O N S >> =================================
		static mat add(mat res, int row_changed, int row_added, double val = 1)
	{
		chk_rng(row_changed, res.row); chk_rng(row_added, res.row);
		for (register int i = 0; i < res.col; i++)
			res.Mat[row_changed][i] += res.Mat[row_added][i] * val;
		return res;
	}
		static mat mult(mat res, int row, double val)
	{
		chk_rng(row, res.row);
		for (register int i = 0; i < res.col; i++)
			res.Mat[row][i] *= val;
		return res;
	}
		static mat swap(mat res, int row1, int row2)
	{
		chk_rng(row1, res.row); chk_rng(row2, res.row);
		double temp;
		for (register int i = 0; i < res.col; i++)
		{
			temp = res.Mat[row1][i];
			res.Mat[row1][i] = res.Mat[row2][i];
			res.Mat[row2][i] = temp;
		}
		return res;
	}
	
		void add(int row_changed, int row_added, double val = 1)
	{
		chk_rng(row_changed, row); chk_rng(row_added, row);
		for (register int i = 0; i < col; i++)
			Mat[row_changed][i] += Mat[row_added][i] * val;
	}
		void mult(int Row, double val)
	{
		chk_rng(Row, row);
		for (register int i = 0; i < col; i++)
			Mat[Row][i] *= val;
	}
		void swap(int row1, int row2)
	{
		chk_rng(row1, row); chk_rng(row2, row);
		double* temp = Mat[row1];
		Mat[row1] = Mat[row2];
		Mat[row2] = temp;
	}
	
		void _add(int row_changed, int row_added, double val = 1)
	{
		chk_rng(row_changed, row); chk_rng(row_added, row);
		if (row_changed < row_added)
			Mat[row_changed][row_added] = val;
		else
			Mat[row_changed][row_added] = 0;
		for (register int i = row_added + 1; i < col; i++)
			Mat[row_changed][i] += Mat[row_added][i] * val;
	}
		void _mult(int Row, double val)
	{
		chk_rng(Row, row);
		if (Row < col)
		{
			Mat[Row][Row] = 1;
			for (register int i = Row + 1; i < col; i++)
				Mat[Row][i] *= val;
		}
	}
		//========================== << R O W   O P E R A T I O N S \E.N.D >> =====================================
		void setelement(int n, int m, double val) { chk_rng(n, row); chk_rng(m, col); Mat[n][m] = val; }
		void freadf(std::string filename) //   TO-DO
	{
		std::ifstream file(filename);
		file >> row >> col;
		if (Mat != nullptr) this->~mat();
		Mat = new double*[row];
		for (register int i = 0; i < row; i++)
			Mat[i] = new double[col];

		for (register int i = 0; i < row; i++)
			for (register int j = 0; j < col; j++)
				file >> Mat[i][j];
	}
		void print() const
	{
		printf("\n");
		std::string* buff;
		buff = new std::string[row*col];
		int* wid = new int[col] {};
		for (register int i = 0; i < row; i++)
			for (register int j = 0; j < col; j++)
			{
				std::string str = "";
				if (std::fabs(Mat[i][j]) < 0.00005)
					str = "0";
				else if (Mat[i][j] == (int)Mat[i][j])
					str = std::to_string((int)Mat[i][j]);
				else
				{
					str = std::to_string(Mat[i][j]);
					for (register int k = str.length() - 1; k > 0; k--)
						if (str.find('.') != -1 && (str.at(k) == '0' || str.at(k) == '.'))
							str.pop_back();
						else
							break;
				}
				if (str.length() > wid[j])
					wid[j] = str.length();
				buff[i*col + j] = str;
			}
		for (register int i = 0; i < row; i++)
		{
			for (register int j = 0; j < col; j++)
			{
				pad(buff[i*col + j], wid[j] - buff[i*col + j].length() + 3);
				std::cout << buff[i*col + j];
			}
			printf("\n");
		}
		printf("\n");
	}
		void fprintf(std::string filename) const
	{
		std::fstream file(filename + ".txt");
		file << row << " " << col << "\n";
		std::string* buff;
		buff = new std::string[row*col];
		int* wid = new int[col] {};
		//        F*O*R*M*A*T*T*I*N*G
		for (register int i = 0; i < row; i++)
			for (register int j = 0; j < col; j++)
			{
				std::string str = "";
				if (std::fabs(Mat[i][j]) < 0.000000001)
					str = "0";
				else if (Mat[i][j] == (int)Mat[i][j])
					str = std::to_string((int)Mat[i][j]);
				else
				{
					str = std::to_string(Mat[i][j]);
			/* ->*/for (register int k = str.length() - 1; k > 0; k--)
			/* | */		if (str.find('.') != -1 && (str.at(k) == '0' || str.at(k) == '.'))
			/* ^ */			str.pop_back();
			/* | */		else
			/* ---<<---*/	break;
				}
				if (str.length() > wid[j])
					wid[j] = str.length();
				buff[i*col + j] = str;
			}
		//         S*T*R*E*A*M*I*N*G
		for (register int i = 0; i < row; i++)
		{
			for (register int j = 0; j < col; j++)
			{
				pad(buff[i*col + j], wid[j] - buff[i*col + j].length() + 3);
				file << buff[i*col + j];
			}
			file << "\n";
		}
		file.close();
	}
		friend void display(const mat& m1, const mat& m2);
		friend void move(mat& sender, mat& reciever);
		void erase() { Mat = nullptr; }
		~mat()
	{
		if (Mat != nullptr)
		{
			for (register int i = 0; i < row; i++)
				delete[] Mat[i];
			delete[] Mat;
		}
	}
	};
	//int mat::constructed;
	std::ostream& operator<<(std::ostream& C, const mat& m) { m.print(); return C; }
	void display(const mat& m1, const mat& m2)
{
	printf("\n");
	std::string *buff1, *buff2;
	buff1 = new std::string[m1.row*m1.col];
	buff2 = new std::string[m2.row*m2.col];
	int *wid1, *wid2;
	wid1 = new int[m1.col]{};
	wid2 = new int[m2.col]{};
	for (register int i = 0; i < m1.row; i++)
		for (register int j = 0; j < m1.col; j++)
		{
			std::string str = "";
			if (std::fabs(m1.Mat[i][j]) < 0.000000001)
				str = "0";
			else if (m1.Mat[i][j] == (int)m1.Mat[i][j])
				str = std::to_string((int)m1.Mat[i][j]);
			else
			{
				str = std::to_string(m1.Mat[i][j]);
				for (register int k = str.length() - 1; k > 0; k--)
					if (str.find('.') != -1 && (str.at(k) == '0' || str.at(k) == '.'))
						str.pop_back();
					else
						break;
			}
			if (str.length() > wid1[j])
				wid1[j] = str.length();
			buff1[i*m1.col + j] = str;
		}
	for (register int i = 0; i < m2.row; i++)
		for (register int j = 0; j < m2.col; j++)
		{
			std::string str = "";
			if (std::fabs(m2.Mat[i][j]) < 0.000000001)
				str = "0";
			else if (m2.Mat[i][j] == (int)m2.Mat[i][j])
				str = std::to_string((int)m2.Mat[i][j]);
			else
			{
				str = std::to_string(m2.Mat[i][j]);
				for (register int k = str.length() - 1; k > 0; k--)
					if (str.find('.') != -1 && (str.at(k) == '0' || str.at(k) == '.'))
						str.pop_back();
					else
						break;
			}
			if (str.length() > wid2[j])
				wid2[j] = str.length();
			buff2[i*m2.col + j] = str;
		}
	int R, r;
	r = std::fmin(m1.row, m2.row);
	R = std::fmax(m1.row, m2.row);
	for (register int i = 0; i < r; i++)
	{
		for (register int j = 0; j < m1.col; j++)
		{
			pad(buff1[i*m1.col + j], wid1[j] - buff1[i*m1.col + j].length() + 3);
			std::cout << buff1[i*m1.col + j];
		}
		printf("|  ");
		for (register int j = 0; j < m2.col; j++)
		{
			pad(buff2[i*m2.col + j], wid2[j] - buff2[i*m2.col + j].length() + 3);
			std::cout << buff2[i*m2.col + j];
		}
		printf("\n");
	}
	if (m1.row > m2.row)
		for (register int i = r; i < R; i++)
		{
			for (register int j = 0; j < m1.col; j++)
			{
				pad(buff1[i*m1.col + j], wid1[j] - buff1[i*m1.col + j].length() + 3);
				std::cout << buff1[i*m1.col + j];
			}
			printf("|\n");
		}
	else if (m2.row > m1.row)
	{
		for (register int i = r; i < R; i++)
		{
			for (register int j = 0; j < m1.col; j++)
			{
				for (register int k = 0; k < wid1[j] + 3; k++)
					printf(" ");
				printf("|  ");
			}
			for (register int j = 0; j < m2.col; j++)
			{
				pad(buff1[i*m1.col + j], wid1[j] - buff1[i*m1.col + j].length() + 3);
				std::cout << buff1[i*m1.col + j];
			}
			printf("\n");
		}
	}
	printf("\n\n");
}
	void move(mat& sender, mat& reciever)
{
	reciever.row = sender.row;
	reciever.col = sender.col;
	reciever.Mat = sender.Mat;
	sender.Mat = nullptr;
}
	
	mat operator+(double val, const mat& m) { return mat::I(std::fmin(m.Row(), m.Col()))*val + m; }
	mat operator-(double val, const mat& m) { return mat::I(std::fmin(m.Row(), m.Col()))*val - m; }
	mat operator*(double val, const mat& m) { return m * val; }
	mat operator/(double val, const mat& m)
{
	int row = m.Row(), col = m.Col();
	mat res(row, col);
	for (register int i = 0; i < row; i++)
		for (register int j = 0; j < col; j++)
		{
			zero_check(m.get(i, j));
			res.setelement(i, j, val / m.get(i, j));
		}
	return res;
}
	real::vec operator*(const real::vec& v, const mat& m)
{
	chk_eq(v.Size(), m.Row());
	int row = m.Row(), col = m.Col();
	real::vec res(col); double sum;
	for (register int i = 0; i < col; i++)
	{
		sum = 0;
		for (register int j = 0; j < row; j++)
			sum += v[j] * m.get(j, i);
		res.set(i, sum);
	}
	return res;
}
	
	double innerProd(const real::vec& v1, const real::vec& v2)
{
	return v1.Wsum(v2);
}
	mat outerProd(const real::vec& v1, const real::vec& v2)
{
	chk_eq(v1.Size(), v2.Size());
	int n = v1.Size();
	mat res(n);
	for (register int i = 0; i < n; i++)
		for (register int j = 0; j < n; j++)
			res.setelement(i, j, v1[i] * v2[j]);
	return res;
}
	mat _set(const mat& m, double(*F)(double x))
{
	int row = m.Row(), col = m.Col();
	mat res(row, col);
	for (register int i = 0; i < row; i++)
		for (register int j = 0; j < col; j++)
			res.setelement(i, j, (*F)(m.get(i, j)));
	return res;
}
	mat abs(const mat& m) { return _set(m, std::fabs); }
	mat ceil(const mat& m) { return _set(m, std::ceil); }
	mat floor(const mat& m) { return _set(m, std::floor); }
	mat ring(int n, char c = 'm')
{
	mat res(n);
	if (c == 'm')
	{
		for (register int i = 0; i < n; i++)
			for (register int j = 0; j < n; j++)
				res.setelement(i, j, (i * j) % n);
	}
	else if (c == 'a')
	{
		for (register int i = 0; i < n; i++)
			for (register int j = 0; j < n; j++)
				res.setelement(i, j, (i + j) % n);
	}
	else
		error("Please Enter a Valid Character\nError Called from ring()");
	return res;
}
	
	// SOLVING A LINEAR SYSTEM Ax=b WITH TRIANGULAR MATRICES:
	real::vec /*lower triangular*/forward_substitution(const mat& A, const real::vec& b, bool pass = false)
{
	int n = A.Row();
	if (!pass)
	{
		if (!A.is_lowerTriangular() && !A.is_square())
			error("Can't Perform Forward Substitution\nNot a Lower-Triangular Matrix");
		chk_eq(A.Col(), b.Size());
		for (register int i = 0; i < n; i++)
			if (A.get(i, i) == 0)
			{
				printf("\n______________________________________________________________________________________________________\n");
				A.print();
				error("Zero Pivot Found");
			}
	}
	real::vec x(n); double sum;
	for (register int i = 0; i < n; i++)
	{
		sum = 0;
		for (register int j = 0; j < i; j++)
			sum += A.get(i, j) * x[j];
		x.set(i, (b[i] - sum) / A.get(i, i));
	}
	return x;
}
	real::vec /*upper triangular*/backward_substitution(const mat& A, const real::vec& b, bool pass = false)
{
	int n = A.Row();
	if (!pass)
	{
		if (!A.is_upperTriangular() && !A.is_square())
			error("Can't Perform Forward Substitution\nNot an Upper-Triangular Matrix");
		chk_eq(A.Col(), b.Size());
		for (register int i = 0; i < n; i++)
			if (A.get(i, i) == 0)
			{
				printf("\n______________________________________________________________________________________________________\n");
				A.print();
				error("Zero Pivot Found\nError Call from backward_substitution()");
			}
	}
	real::vec x(n); double sum;
	for (register int i = n - 1; i >= 0; i--)
	{
		sum = 0;
		for (register int j = n - 1; j > i; j--)
			sum += A.get(i, j) * x[j];
		x.set(i, (b[i] - sum) / A.get(i, i));
	}
	return x;
}
	real::vec solveLinearSystem(const mat& A, const real::vec& b)
{
	if (!A.is_square())
		error("Not a Square Linear System");
	mat l, u, p; A.lu(l, u, p); l = p * l;
	int n = A.Row();
	for (register int i = 0; i < n; i++)
		if (u.get(i, i) == 0 || l.get(i, i) == 0)
		{
			std::cout << l << u;
			error("A Singular Linear System");
		}
	real::vec d = forward_substitution(l, p*b, 1);
	return backward_substitution(u, d, 1);
}
}

namespace com
{
	class mat
	{
		int row, col; Com** Mat = nullptr;
		//int id; static int constructed;
	public:
		mat() {}
		explicit mat(int n)
		{
			//=================================================|
			//construct(constructed, id, "Square Matrix");=====|
			//=================================================|
			row = col = n;
			Mat = new Com*[n];
			for (register int i = 0; i < n; i++)
				Mat[i] = new Com[n] {};
		}
		explicit mat(int n, int m)
		{
			//======================================================|
			//construct(constructed, id, "Rectangular Matrix");=====|
			//======================================================|
			row = n; col = m;
			Mat = new Com*[n];
			for (register int i = 0; i < n; i++)
				Mat[i] = new Com[m] {};
		}
		mat(const com::vec& v)
		{
			//==========================================|
			//construct(constructed, id, "Vector");=====|
			//==========================================|
			row = v.Size(); col = 1;
			Mat = new Com*[row];
			for (register int i = 0; i < row; i++)
				Mat[i] = new Com(v[i]);
		}
		mat(int n, const com::vec& v)
		{
			//===================================================|
			//construct(constructed, id, "Vector Reshape");======|
			//===================================================|
			row = n; col = v.Size() / n;
			Mat = new Com*[row];
			for (register int i = 0; i < row; i++)
				Mat[i] = new Com[col];
			for (register int i = 0; i < row; i++)
				for (register int j = 0; j < col; j++)
					Mat[i][j] = v[i*col + j];
		}
		mat(const mat& other)
		{
			//====================================|
			//id = other.id;                  ====|
			//assign(constructed, id, "Copy");====|
			//====================================|
			row = other.row;
			col = other.col;
			Mat = new Com*[row];
			for (register int i = 0; i < row; i++)
				Mat[i] = new Com[col];
			for (register int i = 0; i < row; i++)
				for (register int j = 0; j < col; j++)
					Mat[i][j] = other.Mat[i][j];
		}
		mat(mat&& other)
		{
			//id = other.id;
			row = other.row;
			col = other.col;
			Mat = other.Mat;
			other.Mat = nullptr;
		}
		mat(const mat* other)
		{
			//=============================================|
			//id = other->id;                         =====|
			//assign(constructed, id, "Pointer Copy");=====|
			//=============================================|
			row = other->row;
			col = other->col;
			Mat = new Com*[row];
			for (register int i = 0; i < row; i++)
				Mat[i] = new Com[col];
			for (register int i = 0; i < row; i++)
				for (register int j = 0; j < col; j++)
					Mat[i][j] = other->Mat[i][j];
		}

		//=============================== << G E T T E R S >> ====================================
		int Row() const { return row; }
		int Col() const { return col; }
		int Size() const { return row * col; }
		int Rank() const
		{
			mat res(reduce());
			int range = std::fmin(row, col);
			for (register int i = 0; i < range; i++)
				if (res.Mat[i][i] == 0)
					return i;
		}
		bool is_symmetric() const
		{
			if (row != col)
				return false;
			for (register int i = 0; i < row; i++)
				for (register int j = 0; j < i; j++)
					if (Mat[i][j] != Mat[j][i])
						return false;
			return true;
		}
		bool is_positive() const
		{
			if (row != col)
				return false;
			mat res(reduce());
			for (register int i = 0; i < res.row; i++)
				if (res.Mat[i][i] <= 0)
					return false;
			return true;
		}
		bool is_semi_positive() const
		{
			if (row != col)
				return false;
			mat res(reduce());
			for (register int i = 0; i < res.row; i++)
				if (res.Mat[i][i] < 0)
					return false;
			return true;
		}
		bool is_triangular() const
		{
			if (Mat[0][1] == 0)
				for (register int i = 0; i < row; i++)
					for (register int j = i + 1; j < col; j++)
					{
						if (Mat[i][j] != 0)
							return false;
					}
			else if (Mat[1][0] == 0)
				for (register int i = 0; i < col; i++)
					for (register int j = i + 1; j < row; j++)
					{
						if (Mat[j][i] != 0)
							return false;
					}
			else
				return false;
			return true;
		}
		bool is_lowerTriangular() const
		{
			for (register int i = 0; i < row; i++)
				for (register int j = i + 1; j < col; j++)
					if (Mat[i][j] != 0)
						return false;
			return true;
		}
		bool is_upperTriangular() const
		{
			for (register int i = 0; i < col; i++)
				for (register int j = i + 1; j < row; j++)
					if (Mat[j][i] != 0)
						return false;
			return true;
		}
		bool is_square() const { return row == col; }
		Com get(int i, int j) const { chk_rng(i, row); chk_rng(j, col); return Mat[i][j]; }
		//=============================== << G E T T E R S \E.N.D>> ===============================

		//=============================== << F U N C T I O N S >> =================================
		Com tr() const
		{
			chk_eq(row, col, "Not a Square Matrix");
			Com sum = 0;
			for (register int i = 0; i < row; i++)
				sum += Mat[i][i];
			return sum;
		}
		mat rsum(matstate state = matstate::non) const
		{
			mat res(row, 1); Com sum;
			int dia = row * (row < col) + col * (col <= row);
			switch (state)
			{
			case matstate::upper:
				for (register int i = 0; i < row; i++)
				{
					sum = 0;
					for (register int j = i; j < col; j++)
						sum += Mat[i][j];
					res.Mat[i][0] = sum;
				}
				return res;
			case matstate::lower:
				for (register int i = 0; i < row; i++)
				{
					sum = 0;
					for (register int j = 0; j <= i; j++)
						sum += Mat[i][j];
					res.Mat[i][0] = sum;
				}
				return res;
			case matstate::diagonal:
				for (register int i = 0; i < dia; i++)
					res.Mat[i][0] = Mat[i][i];
				return res;
			default:
				for (register int i = 0; i < row; i++)
				{
					sum = 0;
					for (register int j = 0; j < col; j++)
						sum += Mat[i][j];
					res.Mat[i][0] = sum;
				}
				return res;
			}
		}
		mat csum(matstate state = matstate::non) const
		{
			mat res(1, col); Com sum;
			int dia = row * (row < col) + col * (col <= row);
			switch (state)
			{
			case matstate::upper:
				for (register int i = 0; i < col; i++)
				{
					sum = 0;
					for (register int j = 0; j <= i; j++)
						sum += Mat[j][i];
					res.Mat[0][i] = sum;
				}
				return res;
			case matstate::lower:
				for (register int i = 0; i < col; i++)
				{
					sum = 0;
					for (register int j = i; j < row; j++)
						sum += Mat[j][i];
					res.Mat[0][i] = sum;
				}
				return res;
			case matstate::diagonal:
				for (register int i = 0; i < dia; i++)
					res.Mat[0][i] = Mat[i][i];
				return res;
			default:
				for (register int i = 0; i < col; i++)
				{
					sum = 0;
					for (register int j = 0; j < row; j++)
						sum += Mat[j][i];
					res.Mat[0][i] = sum;
				}
				return res;
			}
		}
		Com sum(matstate state = matstate::non) const
		{
			Com sum = 0;
			int dia = row * (row < col) + col * (col <= row);
			switch (state)
			{
			case matstate::upper:
				for (register int i = 0; i < row; i++)
					for (register int j = i; j < col; j++)
						sum += Mat[i][j];
				return sum;
			case matstate::lower:
				for (register int i = 0; i < row; i++)
					for (register int j = 0; j <= i; j++)
						sum += Mat[i][j];
				return sum;
			case matstate::diagonal:
				for (register int i = 0; i < dia; i++)
					sum += Mat[i][i];
				return sum;
			default:
				for (register int i = 0; i < row; i++)
					for (register int j = 0; j < col; j++)
						sum += Mat[i][j];
				return sum;
			}
		}
		mat rmax() const
		{
			mat res(row, 1);
			for (register int i = 0; i < row; i++)
			{
				Com max = Mat[i][0];
				for (register int j = 0; j < col; j++)
					if (Mat[i][j] > max)
						max = Mat[i][j];
				res.Mat[i][0] = max;
			}
			return res;
		}
		mat cmax() const
		{
			mat res(1, col);
			for (register int j = 0; j < row; j++)
			{
				Com max = Mat[0][j];
				for (register int i = 0; i < col; i++)
					if (Mat[i][j] > max)
						max = Mat[i][j];
				res.Mat[0][j] = max;
			}
			return res;
		}
		Com max() const
		{
			Com max = Mat[0][0];
			for (register int i = 0; i < row; i++)
				for (register int j = 1; j < col; j++)
					if (Mat[i][j] > max)
						max = Mat[i][j];
			return max;
		}
		mat rmin() const
		{
			mat res(row, 1);
			for (register int i = 0; i < row; i++)
			{
				Com max = Mat[i][0];
				for (register int j = 0; j < col; j++)
					if (Mat[i][j] < max)
						max = Mat[i][j];
				res.Mat[i][0] = max;
			}
			return res;
		}
		mat cmin() const
		{
			mat res(1, col);
			for (register int j = 0; j < row; j++)
			{
				Com max = Mat[0][j];
				for (register int i = 0; i < col; i++)
					if (Mat[i][j] < max)
						max = Mat[i][j];
				res.Mat[0][j] = max;
			}
			return res;
		}
		Com min() const
		{
			Com max = Mat[0][0];
			for (register int i = 0; i < row; i++)
				for (register int j = 1; j < col; j++)
					if (Mat[i][j] < max)
						max = Mat[i][j];
			return max;
		}
		Com det(matstate state = matstate::non) const
		{
			chk_eq(row, col, "No Determinant Exists\nNot a Square Matrix");
			mat temp(this);
			Com det = 1;
			if (Mat[0][0] != 0) det = Mat[0][0];
			switch (state)
			{
			case matstate::upper:
				for (register int i = 0; i < row; i++)
					det *= Mat[i][i];
				return det;
			case matstate::lower:
				return 1;
			case matstate::diagonal:
				for (register int i = 0; i < row; i++)
					det *= Mat[i][i];
				return det;
			case matstate::bidiagonal_u:
				for (register int i = 0; i < row; i++)
					det *= Mat[i][i];
				return det;
			case matstate::bidiagonal_l:
				for (register int i = 0; i < row; i++)
					det *= Mat[i][i];
				return det;
			case matstate::tridiagonal:
				for (register int i = 0; i < row - 1; i++)
				{
					temp._add(i + 1, i, -temp.Mat[i + 1][i] / temp.Mat[i][i]);
					det *= temp.Mat[i][i];
				}
				return det * temp.Mat[row - 1][row - 1];
//===============================================================================================================================
			default:
				short sign = 1;
				for (register int i = 0; i < row-1; i++)
				{
					if (temp.Mat[i][i] == 0)
						for (register int j = i + 1; j < row; j++)
							if (temp.Mat[j][i] == 0 && j == row - 1)
								return 0;
							else if (temp.Mat[j][i] != 0)
							{
								temp.swap(j, i);
								sign *= -1;
								break;
							}
					for (register int j = i + 1; j < col; j++)
						if (temp.Mat[j][i] != 0)
							temp._add(j, i, -temp.Mat[j][i] / temp.Mat[i][i]);
					det *= temp.Mat[i+1][i+1];
				}
				return det * sign;
//===============================================================================================================================
			}
		}
		mat inv(matstate state = matstate::non) const
		{
			if (row != col)
			{
				mat Inv(~*this);
				row > col ? Inv = !(Inv**this)*Inv : Inv = Inv * !(*this*Inv);
				return Inv;
			}
			mat res(mat::I(row));
			mat temp(this);
			Com mag;
			switch (state)
			{
			case matstate::upper:
				for (register int i = 0; i < row; i++)
					if (res.Mat[i][i] == 0)
						error("Matrix is singular\nInvoked from inv()");
					else
						res._mult(i, 1 / res.Mat[i][i]);
				for (register int i = row - 1; i >= 0; i--)
					for (register int j = i - 1; j >= 0; j--)
						if (temp.Mat[j][i] != 0)
							res._add(j, i, -temp.Mat[j][i] / temp.Mat[i][i]);
				return res;

			case matstate::lower:
				for (register int i = 0; i < col; i++)
					for (register int j = i + 1; j < row; j++)
						if (temp.Mat[j][i] != 0)
							res._add(j, i, -temp.Mat[j][i]);
				return res;
			case matstate::diagonal:
				for (register int i = 0; i < row; i++)
					if (res.Mat[i][i] == 0)
						error("Matrix is singular\nInvoked from inv()");
					else
						res.Mat[i][i] = 1 / res.Mat[i][i];
				return res;

			case matstate::orthogonal:
				for (register int i = 0; i < col; i++)
				{
					mag = 0;
					for (register int j = 0; j < row; j++)
						mag += Mat[j][i] * Mat[j][i];
					mag = com::sqrt(mag);
					for (register int j = 0; j < row; j++)
						res.Mat[j][i] = Mat[j][i] / mag;
				}
				return ~res;

			case matstate::orthonormal:
				return ~*this;

			default:
				for (register int i = 0; i < col; i++)
				{
					if (temp.Mat[i][i] == 0)//_________________swap&zero-check_______________
					{																   ////||
						bool single = true;											   ////||
						for (register int j = i + 1; j < row; j++)					   ////||
							if (temp.Mat[j][i] != 0)								   ////||
							{														   ////||
								res.swap(j, i);									       ////||
								temp.swap(j, i);									   ////||
								single = false;										   ////||
								break;												   ////||
							}														   ////||
						if (single)													   ////||
							error("Matrix is singular\nInvoked from inv()");		   ////||
					}//________________________________________________________________////||
					if (temp.Mat[i][i] != 1)
					{
						res.mult(i, 1 / temp.Mat[i][i]);
						temp._mult(i, 1 / temp.Mat[i][i]);
					}
					for (register int j = i + 1; j < row; j++)
						if (j == i || temp.Mat[j][i] == 0)
							continue;
						else
						{
							res.add(j, i, -temp.Mat[j][i]);
							temp._add(j, i, -temp.Mat[j][i]);
						}
				}
				for (register int i = row - 1; i >= 0; i--)
					for (register int j = i - 1; j >= 0; j--)
						if (temp.Mat[j][i] != 0)
							res.add(j, i, -temp.Mat[j][i]);
				return res;
			}
		}
		//---------------------------------N O R M S--------------------------------|||
		Com Norm(matstate state = matstate::non) const { return csum(state).max(); }
		Com E_Norm() const
		{
			Com sum = 0;
			for (register int i = 0; i < row; i++)
				for (register int j = 0; j < col; j++)
					sum += Mat[i][j] * Mat[i][j];
			return com::sqrt(sum);
		}
		Com P_Norm(int p) const
		{
			Com sum = 0;
			for (register int i = 0; i < row; i++)
				for (register int j = 0; j < col; j++)
					sum += com::pow(Mat[i][j], p);
			return com::pow(sum, 1 / p);
		}
		Com inf_Norm() const
		{
			Com* sum = new Com[row];
			for (register int i = 0; i < row; i++)
				for (register int j = 0; j < col; j++)
					sum[i] += com::abs(Mat[i][j]);
			Com max = sum[0];
			for (register int i = 0; i < row; i++)
				if (sum[i] > max)
					max = sum[i];
			return max;
		}
		Com cond(matstate state = matstate::non) const { return Norm(state) * inv(state).Norm(); }
		Com E_cond() const { return E_Norm() * (!*this).E_Norm(); }
		Com P_cond(int p) const { return P_Norm(p) * (!*this).P_Norm(p); }
		Com inf_cond() const { return inf_Norm() * (!*this).inf_Norm(); }
		//--------------------------------------------------------------------------|||
		mat add(Com val) const
		{
			mat res(this);
			for (register int i = 0; i < row; i++)
				for (register int j = 0; j < col; j++)
					res.Mat[i][j] += val;
			return res;
		}
		mat reduce() const
		{
			mat res(this);
			if (row > col)
				res._setrows(col, row, Com(0));
			int range = std::fmin(row, col);
			for (register int i = 0; i < range; i++)
			{
				bool skip = true;
				if (res.Mat[i][i] == 0)
				{
					for (register int j = i + 1; j < row; j++)
						if (res.Mat[j][i] != 0)
						{
							res.swap(i, j);
							res._mult(j, -1);
							break;
						}
						else if (res.Mat[j][i] == 0 && j == row - 1)
							skip = false;//to skip from the next loop
				}
				for (register int j = i + 1; j < range; j++)
					if (res.Mat[j][i] != 0 && skip)
						res._add(j, i, -res.Mat[j][i] / res.Mat[i][i]);
			} return res;
		}
		mat setrows(int firstrow, int lastrow, Com val) const
		{
			chk_rng(lastrow, row);
			mat res(this);
			for (register int i = firstrow; i < lastrow; i++)
				for (register int j = 0; j < col; j++)
					res.Mat[i][j] = val;
			return res;
		}
		mat setrows(int firstrow, int lastrow, const com::vec& v) const
		{
			chk_rng(firstrow, row); chk_rng(lastrow, row);
			chk_eq((lastrow - firstrow)*col, v.Size(), "Vector Size doesn't Match the Rows' Size\nError Call from setrows()");
			mat res(this);
			for (register int i = firstrow; i < lastrow + 1; i++)
				for (register int j = 0; j < col; j++)
					res.Mat[i][j] = v[(i - firstrow)*col + j];
			return res;
		}
		mat setrow(int r, const com::vec& v) const
		{
			mat res(this);
			chk_eq(col, v.Size());
			for (register int i = 0; i < col; i++)
				res.Mat[r][i] = v[i];
			return res;
		}
		mat setcols(int firstcol, int lastcol, Com val) const
		{
			chk_rng(lastcol, col);
			mat res(this);
			for (register int i = 0; i < col; i++)
				for (register int j = firstcol; j < lastcol; j++)
					res.Mat[i][j] = val;
			return res;
		}
		mat setcols(int firstcol, int lastcol, const com::vec& v) const
		{
			chk_rng(firstcol, col); chk_rng(lastcol, col);
			chk_eq((lastcol - firstcol)*row, v.Size(), "Vector Size doesn't Match the Columns' Size\nError Call from setcols()");
			mat res(this);
			for (register int i = 0; i < row; i++)
				for (register int j = firstcol; j < lastcol + 1; j++)
					res.Mat[i][j] = v[i*col + j - firstcol];
			return res;
		}
		mat setcol(int c, const com::vec& v) const
		{
			mat res(this);
			chk_eq(row, v.Size());
			for (register int i = 0; i < row; i++)
				res.Mat[i][c] = v[i];
			return res;
		}
		mat augment(const mat& m) const
		{
			chk_eq(m.row, row);
			mat res(row, col + m.col);
			for (register int i = 0; i < res.row; i++)
				for (register int j = 0; j < res.col; j++)
					j < col ? res.Mat[i][j] = Mat[i][j] : res.Mat[i][j] = m.Mat[i][j - col];
			return res;
		}
		mat append(const mat& m) const
		{
			chk_eq(m.col, col);
			mat res(row + m.row, col);
			for (register int i = 0; i < row; i++)
				for (register int j = 0; j < col; j++)
					res.Mat[i][j] = Mat[i][j];
			for (register int i = row; i < res.row; i++)
				for (register int j = 0; j < col; j++)
					res.Mat[i][j] = m.Mat[i - row][j];
			return res;
		}
		mat augment(const com::vec& v) const
		{
			chk_eq(row, v.Size());
			mat res(row, col + 1);
			for (register int i = 0; i < row; i++)
			{
				for (register int j = 0; j < col; j++)
					res.Mat[i][j] = Mat[i][j];
				res.Mat[i][col] = v[i];
			}
			return res;
		}
		mat append(const com::vec& v) const
		{
			chk_eq(col, v.Size());
			mat res(row + 1, col);
			for (register int i = 0; i < row; i++)
			{
				for (register int j = 0; j < col; j++)
					res.Mat[i][j] = Mat[i][j];
				res.Mat[row][i] = v[i];
			}
			return res;
		}
		mat diagonal() const
		{
			mat res(row, col);
			int range = std::fmin(row, col);
			for (register int i = 0; i < range; i++)
				res.Mat[i][i] = Mat[i][i];
			return res;
		}
		mat off_diagonal() const
		{
			mat res(this);
			int range = std::fmin(row, col);
			for (register int i = 0; i < range; i++)
				res.Mat[i][i] = 0;
			return res;
		}
		mat submat(int firstrow, int lastrow, int firstcol, int lastcol) const
		{
			chk_rng(lastrow, row); chk_rng(lastcol, col);
			mat res(lastrow - firstrow + 1, lastcol - firstcol + 1);
			for (register int i = firstrow; i <= lastrow; i++)
				for (register int j = firstcol; j <= lastcol; j++)
					res.Mat[i - firstrow][j - firstcol] = Mat[i][j];
			return res;
		}
		mat submat(int lastrow, int lastcol) const
		{
			chk_rng(lastrow, row); chk_rng(lastcol, col);
			mat res(lastrow + 1, lastcol + 1);
			for (register int i = 0; i <= lastrow; i++)
				for (register int j = 0; j <= lastcol; j++)
					res.Mat[i][j] = Mat[i][j];
			return res;
		}
		mat submat(int pivot) const
		{
			chk_rng(pivot, row); chk_rng(pivot, col);
			mat res(pivot + 1, pivot + 1);
			for (register int i = 0; i <= pivot; i++)
				for (register int j = 0; j <= pivot; j++)
					res.Mat[i][j] = Mat[i][j];
			return res;
		}
		com::vec unroll() const
		{
			com::vec res(row*col);
			for (register int i = 0; i < row; i++)
				for (register int j = 0; j < col; j++)
					res.set(i*col + j, Mat[i][j]);
			return res;
		}
		com::vec getrow(int n) const
		{
			chk_rng(n, row);
			com::vec res(col);
			for (register int i = 0; i < col; i++)
				res.set(i, Mat[n][i]);
			return res;
		}
		com::vec getcol(int n) const
		{
			chk_rng(n, col);
			com::vec res(row);
			for (register int i = 0; i < row; i++)
				res.set(i, Mat[i][n]);
			return res;
		}

		void _add(Com val)
		{
			for (register int i = 0; i < row; i++)
				for (register int j = 0; j < col; j++)
					Mat[i][j] += val;
		}
		void _reduce()
		{
			if (row > col)
				_setrows(col - 1, row - 1, Com(0));
			int range = col * (row > col) + (row - 1)*(row <= col);
			bool sign = false;
			for (register int i = 0; i < range; i++)
			{
				if (Mat[i][i] == 0)
				{
					bool skip = false;
					for (register int k = i; k < col; k++)
					{
						for (register int j = i + 1; j < row; j++)
							if (Mat[j][k] != 0)
							{
								swap(i, j);
								skip = true;
								sign = !sign;
								break;
							}
						if (skip)
							break;
					}
				}
				for (register int j = i + 1; j < range; j++)
					if (Mat[j][i] != 0)
						_add(j, i, -Mat[j][i] / Mat[i][i]);
			}
			if (sign)
				_mult(row - 1, -1);
		}
		void _diagonal()
		{
			_reduce();
			int range = std::fmin(row, col);
			for (register int i = 0; i < range - 1; i++)
				for (register int j = i + 1; j < col; j++)
					Mat[i][j] = 0;
		}
		void _setrows(int firstrow, int lastrow, Com val)
		{
			chk_rng(lastrow, row);
			for (register int i = firstrow; i <= lastrow; i++)
				for (register int j = 0; j < col; j++)
					Mat[i][j] = val;
		}
		void _setrows(int firstrow, int lastrow, const com::vec& v)
		{
			chk_rng(firstrow, row); chk_rng(lastrow, row);
			chk_eq((lastrow - firstrow)*col, v.Size(), "Vector Size doesn't Match the Rows' Size\nError Call from setrows()");
			for (register int i = firstrow; i < lastrow + 1; i++)
				for (register int j = 0; j < col; j++)
					Mat[i][j] = v[(i - firstrow)*col + j];
		}
		void _setrow(int r, const com::vec& v)
		{
			chk_eq(col, v.Size()); chk_rng(r, row);
			for (register int i = 0; i < col; i++)
				Mat[r][i] = v[i];
		}
		void _setcols(int firstcol, int lastcol, Com val)
		{
			chk_rng(lastcol, col);
			for (register int i = 0; i <= col; i++)
				for (register int j = firstcol; j < lastcol; j++)
					Mat[i][j] = val;
		}
		void _setcols(int firstcol, int lastcol, const com::vec& v)
		{
			chk_rng(firstcol, col); chk_rng(lastcol, col);
			chk_eq((lastcol - firstcol)*row, v.Size(), "Vector Size doesn't Match the Columns' Size\nError Call from setcols()");
			for (register int i = 0; i < row; i++)
				for (register int j = firstcol; j < lastcol + 1; j++)
					Mat[i][j] = v[i*col + j - firstcol];
		}
		void _setcol(int c, const com::vec& v)
		{
			chk_eq(row, v.Size()); chk_rng(c, col);
			for (register int i = 0; i < row; i++)
				Mat[i][c] = v[i];
		}
		void _swapcols(int first, int second)
		{
			chk_rng(first, col); chk_rng(second, col);
			Com temp;
			for (register int i = 0; i < row; i++)
			{
				temp = Mat[i][first];
				Mat[i][first] = Mat[i][second];
				Mat[i][second] = temp;
			}
		}
		void _augment(const mat& m)
		{
			chk_eq(row, m.row);
			mat res(row, col + m.col);
			for (register int i = 0; i < row; i++)
			{
				for (register int j = 0; j < col; j++)
					res.Mat[i][j] = Mat[i][j];
				for (register int j = 0; j < m.col; j++)
					res.Mat[i][col + j] = m.Mat[i][j];
			}
			col = res.col;
			Mat = res.Mat;
			res.Mat = nullptr;
		}
		void _augment(const com::vec& v)
		{
			chk_eq(row, v.Size());
			mat res(row, col + 1);
			for (register int i = 0; i < row; i++)
			{
				for (register int j = 0; j < col; j++)
					res.Mat[i][j] = Mat[i][j];
				res.Mat[i][col] = v[i];
			}
			col++;
			Mat = res.Mat;
			res.Mat = nullptr;
		}
		void _append(const mat& m)
		{
			chk_eq(col, m.col);
			mat res(row + m.row, col);
			for (register int i = 0; i < col; i++)
			{
				for (register int j = 0; j < row; j++)
					res.Mat[j][i] = Mat[j][i];
				for (register int j = 0; j < m.row; j++)
					res.Mat[row + j][i] = m.Mat[j][i];
			}
			row = res.row;
			Mat = res.Mat;
			res.Mat = nullptr;
		}
		void _append(const com::vec& v)
		{
			chk_eq(col, v.Size());
			mat res(row + 1, col);
			for (register int i = 0; i < col; i++)
			{
				for (register int j = 0; j < row; j++)
					res.Mat[j][i] = Mat[j][i];
				res.Mat[row][i] = v[i];
			}
			row++;
			Mat = res.Mat;
			res.Mat = nullptr;
		}
		void _submat(int firstrow, int lastrow, int firstcol, int lastcol)
		{
			chk_rng(lastrow, row); chk_rng(lastcol, col);
			mat res(lastrow - firstrow, lastcol - firstcol);
			for (register int i = firstrow; i < lastrow; i++)
				for (register int j = firstcol; j < lastcol; j++)
					res.Mat[i - firstrow][j - firstcol] = Mat[i][j];
			row = res.row;
			col = res.col;
			Mat = res.Mat;
			res.Mat = nullptr;
		}
		void _make_symmetric()
		{
			chk_eq(row, col, "Matrix is not Square\nCan't be Made Symmetric");
			for (register int i = 0; i < row; i++)
				for (register int j = i + 1; j < col; j++)
					Mat[j][i] = Mat[i][j];
		}
		//=============================== << F A C T O R I Z A T I O N >> =================================
		list<mat> lu() const
		{
			mat l(mat::I(row)), u(this), p(l);
			int dia = col * (row > col) + (row - 1)*(row >= col);
			for (register int i = 0; i < dia; i++)
			{
				if (u.Mat[i][i] == 0) // if a pivot is zero:
				{
					bool flag = false;
					for (register int j = i; j < col; j++)
					{
						for (register int k = i + 1; k < row; k++)
							if (u.Mat[k][j] != 0)
							{
								u.swap(i, k);
								p.swap(i, k);
								l._swapcols(i, k);
								flag = true;
								break;
							}
						if (flag)
							break;
						else if (j == col - 1)
							goto OVER;
					}
				}
				for (register int j = i + 1; j < row; j++)
				{
					if (u.Mat[j][i] == 0)
						continue;
					l.Mat[j][i] = u.Mat[j][i] / u.Mat[i][i];
					u._add(j, i, -u.Mat[j][i] / u.Mat[i][i]);
				}
			}
			OVER:
			list<mat> res;
			res.add(l);
			res.add(u);
			res.add(~p);
			return res;
		}
		void lu(mat& L, mat& U, mat& P) const
		{
			mat l(mat::I(row)), u(this), p(l);
			int dia = col * (row > col) + (row - 1)*(row >= col);
			for (register int i = 0; i < dia; i++)
			{
				if (u.Mat[i][i] == 0) // if a pivot is zero:
				{
					bool flag = false;
					for (register int j = i; j < col; j++)
					{
						for (register int k = i + 1; k < row; k++)
							if (u.Mat[k][j] != 0)
							{
								u.swap(i, k);
								p.swap(i, k);
								l._swapcols(i, k);
								flag = true;
								break;
							}
						if (flag)
							break;
						else if (j == col - 1)
							goto OVER;
					}
				}
				for (register int j = i + 1; j < row; j++)
				{
					if (u.Mat[j][i] == 0)
						continue;
					l.Mat[j][i] = u.Mat[j][i] / u.Mat[i][i];
					u._add(j, i, -u.Mat[j][i] / u.Mat[i][i]);
				}
			}
			OVER:
			move(l, L);
			move(u, U);
			P = ~p;
		}
		void lu(mat& L, mat& U) const
		{
			mat l(mat::I(row)), u(this), p(l);
			int dia = col * (row > col) + (row - 1)*(row >= col);
			for (register int i = 0; i < dia; i++)
			{
				if (u.Mat[i][i].mag() < 10e-7) // if a pivot is zero:
				{
					bool flag = false;
					for (register int j = i; j < col; j++)
					{
						for (register int k = i + 1; k < row; k++)
							if (u.Mat[k][j] != 0)
							{
								u.swap(i, k);
								p.swap(i, k);
								l._swapcols(i, k);
								flag = true;
								break;
							}
						if (flag)
							break;
						else if (j == col - 1)
							goto OVER;
					}
				}
				for (register int j = i + 1; j < row; j++)
				{
					if (u.Mat[j][i] == 0)
						continue;
					l.Mat[j][i] = u.Mat[j][i] / u.Mat[i][i];
					u._add(j, i, -u.Mat[j][i] / u.Mat[i][i]);
				}
			}
			OVER:
			move(l, L);
			move(u, U);
		}
		list<mat> ldu() const
		{
			mat l(mat::I(row)), d(row), u(this), p(l);
			int dia = col * (row > col) + (row - 1)*(row >= col);
			for (register int i = 0; i < dia; i++)
			{
				if (u.Mat[i][i] == 0) // if a pivot is zero:
				{
					bool flag = false;
					for (register int j = i; j < col; j++)
					{
						for (register int k = i + 1; k < row; k++)
							if (u.Mat[k][j] != 0)
							{
								u.swap(i, k);
								p.swap(i, k);
								l._swapcols(i, k);
								flag = true;
								break;
							}
						if (flag)
							break;
						else if (j == col - 1)
							goto OVER;
					}
				}
				for (register int j = i + 1; j < row; j++)
				{
					if (u.Mat[j][i] == 0)
						continue;
					l.Mat[j][i] = u.Mat[j][i] / u.Mat[i][i];
					u._add(j, i, -u.Mat[j][i] / u.Mat[i][i]);
				}
			}
			OVER:
			dia += row == col;
			for (register int i = 0; i < dia; i++)
			{
				d.Mat[i][i] = u.Mat[i][i];
				u._mult(i, 1 / u.Mat[i][i]);
			}
			list<mat> res;
			res.add(l);
			res.add(d);
			res.add(u);
			res.add(~p);
			return res;
		}
		void ldu(mat& L, mat& D, mat& U, mat& P) const
		{
			mat l(mat::I(row)), d(row), u(this), p(l);
			int dia = col * (row > col) + (row - 1)*(row >= col);
			for (register int i = 0; i < dia; i++)
			{
				if (u.Mat[i][i] == 0) // if a pivot is zero:
				{
					bool flag = false;
					for (register int j = i; j < col; j++)
					{
						for (register int k = i + 1; k < row; k++)
							if (u.Mat[k][j] != 0)
							{
								u.swap(i, k);
								p.swap(i, k);
								l._swapcols(i, k);
								flag = true;
								break;
							}
						if (flag)
							break;
						else if (j == col - 1)
							goto OVER;
					}
				}
				for (register int j = i + 1; j < row; j++)
				{
					if (u.Mat[j][i] == 0)
						continue;
					l.Mat[j][i] = u.Mat[j][i] / u.Mat[i][i];
					u._add(j, i, -u.Mat[j][i] / u.Mat[i][i]);
				}
			}
			OVER:
			dia += row == col;
			for (register int i = 0; i < dia; i++)
			{
				d.Mat[i][i] = u.Mat[i][i];
				u._mult(i, 1 / u.Mat[i][i]);
			}
			move(l, L);
			move(d, D);
			move(u, U);
			p = ~p;
		}
		void ldu(mat& L, mat& D, mat& U) const
		{
			mat l(mat::I(row)), d(row), u(this), p(l);
			int dia = col * (row > col) + (row - 1)*(row <= col);
			for (register int i = 0; i < dia; i++)
			{
				if (u.Mat[i][i] == 0) // if a pivot is zero:
				{
					bool flag = false;
					for (register int j = i; j < col; j++)
					{
						for (register int k = i + 1; k < row; k++)
							if (u.Mat[k][j] != 0)
							{
								u.swap(i, k);
								p.swap(i, k);
								l._swapcols(i, k);
								flag = true;
								break;
							}
						if (flag)
							break;
						else if (j == col - 1)
							goto OVER;
					}
				}
				for (register int j = i + 1; j < row; j++)
				{
					if (u.Mat[j][i] == 0)
						continue;
					l.Mat[j][i] = u.Mat[j][i] / u.Mat[i][i];
					u._add(j, i, -u.Mat[j][i] / u.Mat[i][i]);
				}
			}
			OVER:
			dia += row == col;
			for (register int i = 0; i < dia; i++)
			{
				d.Mat[i][i] = u.Mat[i][i];
				u._mult(i, 1 / u.Mat[i][i]);
			}
			L = (~p)*l;
			move(d, D);
			move(u, U);
		}

		pair<mat, mat> qr() const
		{
			mat q(row, col), r(col);
			com::vec* a; a = new com::vec[col];
			for (register int i = 0; i < col; i++)
				a[i] = getcol(i);
			com::vec sum(row);
			for (register int i = 0; i < col; i++)
			{
				for (register int j = 0; j < i; j++)
				{
					r.Mat[j][i] = a[i].Wsum(a[j]);
					sum += r.Mat[j][i] * a[j];
				}
				a[i] -= sum;
				r.Mat[i][i] = a[i].mag();
				a[i] /= r.Mat[i][i];
				q._setcol(i, a[i]);
				sum.reset();
			}
			return pair<mat, mat>(q, r);
		}
		void qr(mat& Q, mat& R) const
		{
			mat q(row, col), r(col);
			com::vec* a; a = new com::vec[col];
			for (register int i = 0; i < col; i++)
				a[i] = getcol(i);
			com::vec sum(row);
			for (register int i = 0; i < col; i++)
			{
				for (register int j = 0; j < i; j++)
				{
					r.Mat[j][i] = a[i].Wsum(a[j]);
					sum += r.Mat[j][i] * a[j];
				}
				a[i] -= sum;
				r.Mat[i][i] = a[i].mag();
				a[i] /= r.Mat[i][i];
				q._setcol(i, a[i]);
				sum.reset();
			}
			Q.row = q.row; Q.col = q.col; Q.Mat = q.Mat; q.Mat = nullptr;
			R.row = r.row; R.col = r.col; R.Mat = r.Mat; r.Mat = nullptr;
		}

		//    NOT FINISHED YET
		/*pair<mat, mat> eig(matstate state = matstate::non) const
		{
			if (row != col)
				error("Matrix is Rectangular\nError Call from eig()");


			mat temp(reduce()), X(row), lambda(row);
			com::vec eig(row);
			int zeroes;
			for (zeroes = row - 1; zeroes >= 0; zeroes--)
				if (temp.Mat[zeroes][zeroes] > 1.0e-7 || temp.Mat[zeroes][zeroes] < -1.0e-7)
					break;
			zeroes = row - zeroes - 1;
			if (zeroes > 0) // zero eigenvalue found:
			{
				for (register int i = row - zeroes - 1; i > 0; i--)
					for (register int j = i - 1; j >= 0; j--)
						temp._add(j, i, -temp.Mat[j][i] / temp.Mat[i][i]);
				com::vec x(row);
				for (register int i = col - 1; i >= col - zeroes; i--)
				{
					for (register int j = 0; j < row - zeroes; j++)
						x.set(j, -temp.Mat[j][i] / temp.Mat[j][j]);
					x.set(i, 1);
					X._setcol(i, x);
					x.set(i, 0);
				}
			}
			//==================================================
			auto power = [](const mat& A, com::vec& x)->double
			{
				mat res = A * A * A;
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
				return x0;
			};
			auto Power = [](const mat& A, com::vec& x)->double
			{
				double x0 = 0, x1, e = 1;
				while (e > 5.0e-7)
				{
					x1 = x0;
					x = A * x;
					x0 = x.abs_max();
					x._norm();
					e = x0 - x1;
				}
				return x0;
			};
			//==================================================
			if (state == matstate::symmetric || is_symmetric())
			{
				auto deflate = [](mat& A, com::vec& X, double eig)
				{
					mat x(X.hat());
					A = A - x * ~x*eig;
				};
				temp = *this; com::vec x(row);
				//___The Main Loop______________________
				for (register int i = 0; i < row - zeroes; i++)
				{
					x.set(1);
					eig.set(i, power(temp, x));
					lambda.Mat[i][i] = eig[i];
					X._setcol(i, x.hat());
					deflate(temp, x, eig[i]);
				}//|||||||||||||||||||||||||||||||||||||
				return pair<mat, mat>(X, lambda);
			}
			else
			{
				double max, min;
				com::vec x = com::vec::O(row);
				eig.set(0, power(temp, x));
				lambda.Mat[0][0] = max = eig[0];
				X._setcol(0, x.hat());
				if (zeroes > 0)
				{
					min = 0;
					zeroes--;
				}
				else
				{
					x = com::vec::O(row);
					eig.set(0, 1 / power(!temp, x));
					lambda.Mat[0][0] = min = eig[0];
					X._setcol(0, x.hat());
				}
				double alpha = (max + min) / 2;
				for (register int i = 1; i <= row - zeroes - 2; i++)
				{
					x = com::vec::O(row);
					eig.set(i, 1 / power(!(temp - alpha) + alpha, x));
					if (eig[i] != max && eig[i] != min)
					{
						lambda.Mat[i][i] = eig[0];
						X._setcol(i, x.hat());
					}
					matstate state = matstate::non;
				}
			}
		}
		void eig(mat& x, mat& e) const;

		list<mat> svd() const
		{
			mat res1(*this * ~*this), res2(~*this * *this);

		}
		void svd(mat& u, mat& s, mat& v) const
		{

		}
		//*/
		//======================================== << S T A T I C >> ========================================
		static mat I(int n)
		{
			mat res(n);
			for (register int i = 0; i < n; i++)
				res.Mat[i][i] = 1;
			return res;
		}
		static mat O(int n)
		{
			mat res(n);
			for (register int i = 0; i < n; i++)
				for (register int j = 0; j < n; j++)
					res.Mat[i][j] = 1;
			return res;
		}
		static mat O(int n, int m)
		{
			mat res(n, m);
			for (register int i = 0; i < n; i++)
				for (register int j = 0; j < m; j++)
					res.Mat[i][j] = 1;
			return res;
		}
		static mat R(int n)
		{
			mat res(n);
			for (register int i = 0; i < n; i++)
				for (register int j = 0; j < n; j++)
					res.Mat[i][j] = (std::rand() % 20) - 10;//Com((std::rand() % 20) - 10, (std::rand() % 20) - 10);
			return res;
		}
		static mat R(int n, int m)
		{
			mat res(n, m);
			for (register int i = 0; i < n; i++)
				for (register int j = 0; j < m; j++)
					res.Mat[i][j] = Com((std::rand() % 20) - 10, (std::rand() % 20) - 10);
			return res;
		}
		static mat R(int n, Com upperlimit)
		{
			mat res(n);
			for (register int i = 0; i < n; i++)
				for (register int j = 0; j < n; j++)
					res.Mat[i][j] = com::mod(std::rand(), upperlimit);
			return res;
		}
		static mat R(int n, int m, Com upperlimit)
		{
			mat res(n, m);
			for (register int i = 0; i < n; i++)
				for (register int j = 0; j < m; j++)
					res.Mat[i][j] = com::mod(std::rand(), upperlimit);
			return res;
		}
		static mat R(int n, Com lowerlimit, Com upperlimit)
		{
			mat res(n);
			for (register int i = 0; i < n; i++)
				for (register int j = 0; j < n; j++)
					res.Mat[i][j] = com::mod(std::rand(), upperlimit - lowerlimit) + lowerlimit;
			return res;
		}
		static mat R(int n, int m, Com lowerlimit, Com upperlimit)
		{
			mat res(n, m);
			for (register int i = 0; i < n; i++)
				for (register int j = 0; j < m; j++)
					res.Mat[i][j] = com::mod(std::rand(), upperlimit - lowerlimit) + lowerlimit;
			return res;
		}
		static mat upper(int n)
		{
			mat res(n);
			for (register int i = 0; i < n; i++)
				for (register int j = i; j < n; j++)
					res.Mat[i][j] = Com((std::rand() % 20) - 10, (std::rand() % 20) - 10);
			return res;
		}
		static mat lower(int n)
		{
			mat res(n);
			for (register int i = 0; i < n; i++)
			{
				for (register int j = 0; j < i; j++)
					res.Mat[i][j] = Com((std::rand() % 20) - 10, (std::rand() % 20) - 10);
				res.Mat[i][i] = 1;
			}
			return res;
		}
		static mat diagonal(int n)
		{
			mat res(n);
			for (register int i = 0; i < n; i++)
				res.Mat[i][i] = Com((std::rand() % 20) - 10, (std::rand() % 20) - 10);
			return res;
		}
		static mat hermitian(int n)
		{
			mat res(n);
			for (register int i = 0; i < n; i++)
			{
				res.Mat[i][i] = (std::rand() % 20) - 10;
				for (register int j = i+1; j < n; j++)
				{
					res.Mat[i][j] = Com((std::rand() % 20) - 10, (std::rand() % 20) - 10);
					res.Mat[j][i] = ~res.Mat[i][j];
				}
			}
			return res;
		}
		//=============================== << F U N C T I O N S \E.N.D>> ===========================

		//=============================== << O P E R A T O R S >> =================================
		mat operator-() const
		{
			mat res(row, col);
			for (register int i = 0; i < row; i++)
				for (register int j = 0; j < col; j++)
					res.Mat[i][j] = -Mat[i][j];
			return res;
		}
		mat operator~() const
		{
			mat res(col, row);
			for (register int i = 0; i < row; i++)
				for (register int j = 0; j < col; j++)
					res.Mat[j][i] = Mat[i][j];
			return res;
		}
		mat operator!() const
		{
			if (row != col)
			{
				mat inv(~*this);
				return row > col ? !(inv**this)*inv : inv*!(*this*inv);
			}
			mat res(mat::I(row));
			mat temp(this);
			for (register int i = 0; i < col; i++)
			{
				if (temp.Mat[i][i] == 0)//_________________swap&zero-check_______________
				{																   ////||
					bool single = true;											   ////||
					for (register int j = i + 1; j < row; j++)					   ////||
						if (temp.Mat[j][i] != 0)								   ////||
						{														   ////||
							res.swap(j, i);									       ////||
							temp.swap(j, i);									   ////||
							single = false;										   ////||
							break;												   ////||
						}														   ////||
					if (single)													   ////||
						error("Matrix is singular\nInvoked from !-matrix");		   ////||
				}//________________________________________________________________////||
				if (temp.Mat[i][i] != 1)
				{
					res.mult(i, 1 / temp.Mat[i][i]);
					temp._mult(i, 1 / temp.Mat[i][i]);
				}
				for (register int j = i + 1; j < row; j++)
					if (j == i || temp.Mat[j][i] == 0)
						continue;
					else
					{
						res.add(j, i, -temp.Mat[j][i]);
						temp._add(j, i, -temp.Mat[j][i]);
					}
			}
			for (register int i = row - 1; i >= 0; i--)
				for (register int j = i - 1; j >= 0; j--)
					if (temp.Mat[j][i] != 0)
						res.add(j, i, -temp.Mat[j][i]);
			return res;
		}
		//===================================    MAT X MAT    =====================================
		mat operator+(const mat& other) const
		{
			chk_eq(other.row, row); chk_eq(other.col, col);
			mat res(row, col);
			for (register int i = 0; i < row; i++)
				for (register int j = 0; j < col; j++)
					res.Mat[i][j] = Mat[i][j] + other.Mat[i][j];
			return res;
		}
		mat operator-(const mat& other) const
		{
			chk_eq(other.row, row); chk_eq(other.col, col);
			mat res(row, col);
			for (register int i = 0; i < row; i++)
				for (register int j = 0; j < col; j++)
					res.Mat[i][j] = Mat[i][j] - other.Mat[i][j];
			return res;
		}
		mat operator*(const mat& other) const
		{
			chk_eq(col, other.row);
			mat res(row, other.col);
			Com sum;
			for (register int i = 0; i < res.row; i++)
				for (register int j = 0; j < res.col; j++)
				{
					sum = 0;
					for (register int k = 0; k < col; k++)
						sum += Mat[i][k] * other.Mat[k][j];
					res.Mat[i][j] = sum;
				}
			return res;
		}
		mat operator/(const mat& other) const { return *this * !other; }
		mat operator|(const mat& other) const { return augment(other); }
		bool operator==(const mat& other) const
		{
			chk_eq(row, other.row); chk_eq(col, other.col);
			for (register int i = 0; i < row; i++)
				for (register int j = 0; j < col; j++)
					if (Mat[i][j] != other.Mat[i][j])
						return false;
			return true;
		}
		bool operator!=(const mat& other) const { return !(*this == other); }

		mat& operator =(mat&& other)
		{
			//id = other.id;
			row = other.row;
			col = other.col;
			Mat = other.Mat;
			other.Mat = nullptr;
			return *this;
		}
		mat& operator =(const mat& other)
		{
			//===========================================|
			//id = other.id;                        =====|
			//assign(constructed, id, "Assignment");=====|
			//===========================================|
			this->~mat();
			row = other.row;
			col = other.col;
			Mat = new Com*[row];
			for (register int i = 0; i < row; i++)
				Mat[i] = new Com[col];
			for (register int i = 0; i < row; i++)
				for (register int j = 0; j < col; j++)
					Mat[i][j] = other.Mat[i][j];
			return *this;
		}
		mat& operator+=(const mat& other) { *this = *this + other; return *this; }
		mat& operator-=(const mat& other) { *this = *this - other; return *this; }
		mat& operator*=(const mat& other) { *this = *this * other; return *this; }
		mat& operator/=(const mat& other) { *this = *this / other; return *this; }
		mat& operator|=(const mat& other) { *this = *this | other; return *this; }
		//===================================    MAT X VEC    =====================================
		com::vec operator*(const com::vec& v) const
		{
			chk_eq(col, v.Size());
			com::vec res(row);
			Com sum;
			for (register int i = 0; i < row; i++)
			{
				sum = 0;
				for (register int j = 0; j < col; j++)
					sum += Mat[i][j] * v[j];
				res.set(i, sum);
			}
			return res;
		}
		com::vec& operator*=(const com::vec& v)
		{
			*this = *this * v;
			com::vec res(this->unroll());
			return res;
		}

		mat& operator =(const com::vec& other)
		{
			//=====================================================|
			//construct(constructed, id, "Vector Assignment");=====|
			//=====================================================|
			row = other.Size();
			col = 1;
			Mat = new Com*[row];
			for (register int i = 0; i < row; i++)
				Mat[i] = new Com{ other[i] };
			return *this;
		}
		//===================================    MAT X COM    =====================================
		mat operator*(Com val) const
		{
			mat res(row, col);
			for (register int i = 0; i < row; i++)
				for (register int j = 0; j < col; j++)
					res.Mat[i][j] = Mat[i][j] * val;
			return res;
		}
		mat operator/(Com val) const
		{
			zero_check(val);
			mat res(row, col);
			for (register int i = 0; i < row; i++)
				for (register int j = 0; j < col; j++)
					res.Mat[i][j] = Mat[i][j] / val;
			return res;
		}
		mat operator%(Com val) const
		{
			zero_check(val, "Zero Modulus");
			mat res(row, col);
			for (register int i = 0; i < row; i++)
				for (register int j = 0; j < col; j++)
					res.Mat[i][j] = com::mod(Mat[i][j], val);
			return res;
		}
		mat operator+(Com val) const { return (*this) + mat::I(std::fmin(row, col))*val; }
		mat operator-(Com val) const { return (*this) - mat::I(std::fmin(row, col))*val; }

		mat& operator*=(Com val) { *this = *this * val; return *this; }
		mat& operator/=(Com val) { *this = *this / val; return *this; }
		mat& operator%=(Com val) { *this = *this % val; return *this; }
		mat& operator+=(Com val) { *this = *this + val; return *this; }
		mat& operator-=(Com val) { *this = *this - val; return *this; }
		//===================================    MAT X NUM    =====================================
		mat operator*(double val) const
		{
			mat res(row, col);
			for (register int i = 0; i < row; i++)
				for (register int j = 0; j < col; j++)
					res.Mat[i][j] = Mat[i][j] * val;
			return res;
		}
		mat operator/(double val) const
		{
			zero_check(val);
			mat res(row, col);
			for (register int i = 0; i < row; i++)
				for (register int j = 0; j < col; j++)
					res.Mat[i][j] = Mat[i][j] / val;
			return res;
		}
		mat operator%(double val) const
		{
			zero_check(val);
			mat res(row, col);
			for (register int i = 0; i < row; i++)
				for (register int j = 0; j < col; j++)
					res.Mat[i][j] = com::mod(Mat[i][j], val);
			return res;
		}
		mat operator+(double val) const { return (*this) + mat::I(std::fmin(row, col))*val; }
		mat operator-(double val) const { return (*this) - mat::I(std::fmin(row, col))*val; }

		mat& operator*=(double val) { *this = *this * val; return *this; }
		mat& operator/=(double val) { *this = *this / val; return *this; }
		mat& operator%=(double val) { *this = *this % val; return *this; }
		mat& operator+=(double val) { *this = *this + val; return *this; }
		mat& operator-=(double val) { *this = *this - val; return *this; }
		//=============================== << O P E R A T O R S \E.N.D>> ===================================

		//============================= << R O W   O P E R A T I O N S >> =================================
		static mat add(mat res, int row_changed, int row_added, Com val = 1)
		{
			chk_rng(row_changed, res.row); chk_rng(row_added, res.row);
			for (register int i = 0; i < res.col; i++)
				res.Mat[row_changed][i] += res.Mat[row_added][i] * val;
			return res;
		}
		static mat mult(mat res, int row, Com val)
		{
			chk_rng(row, res.row);
			for (register int i = 0; i < res.col; i++)
				res.Mat[row][i] *= val;
			return res;
		}
		static mat swap(mat res, int row1, int row2)
		{
			chk_rng(row1, res.row); chk_rng(row2, res.row);
			Com temp;
			for (register int i = 0; i < res.col; i++)
			{
				temp = res.Mat[row1][i];
				res.Mat[row1][i] = res.Mat[row2][i];
				res.Mat[row2][i] = temp;
			}
			return res;
		}

		void add(int row_changed, int row_added, Com val = 1)
		{
			chk_rng(row_changed, row); chk_rng(row_added, row);
			for (register int i = 0; i < col; i++)
				Mat[row_changed][i] += Mat[row_added][i] * val;
		}
		void mult(int Row, Com val)
		{
			chk_rng(Row, row);
			for (register int i = 0; i < col; i++)
				Mat[Row][i] *= val;
		}
		void swap(int row1, int row2)
		{
			chk_rng(row1, row); chk_rng(row2, row);
			Com* temp = Mat[row1];
			Mat[row1] = Mat[row2];
			Mat[row2] = temp;
		}

		void _add(int row_changed, int row_added, Com val = 1)
		{
			chk_rng(row_changed, row); chk_rng(row_added, row);
			if (row_changed < row_added)
				Mat[row_changed][row_added] = val;
			else
				Mat[row_changed][row_added] = 0;
			for (register int i = row_added + 1; i < col; i++)
				Mat[row_changed][i] += Mat[row_added][i] * val;
		}
		void _mult(int Row, Com val)
		{
			chk_rng(Row, row);
			if (Row < col)
			{
				Mat[Row][Row] = 1;
				for (register int i = Row + 1; i < col; i++)
					Mat[Row][i] *= val;
			}
		}
		//========================== << R O W   O P E R A T I O N S \E.N.D >> =====================================
		void setelement(int n, int m, Com val) { chk_rng(n, row); chk_rng(m, col); Mat[n][m] = val; }
		void freadf(std::string filename) //   TO-DO
		{
			std::ifstream file(filename);
			file >> row >> col;
			if (Mat != nullptr) this->~mat();
			Mat = new Com*[row];
			for (register int i = 0; i < row; i++)
				Mat[i] = new Com[col];

			for (register int i = 0; i < row; i++)
				for (register int j = 0; j < col; j++)
					file >> Mat[i][j].a >> Mat[i][j].b;
		}
		void print() const
		{
			printf("\n");
			std::string* buff;
			buff = new std::string[row*col];
			int* wid = new int[col] {};
			for (register int i = 0; i < row; i++)
				for (register int j = 0; j < col; j++)
				{
					std::string str = com::to_string(Mat[i][j]);
					if (str.length() > wid[j])
						wid[j] = str.length();
					buff[i*col + j] = str;
				}
			for (register int i = 0; i < row; i++)
			{
				for (register int j = 0; j < col; j++)
				{
					pad(buff[i*col + j], wid[j] - buff[i*col + j].length() + 3);
					std::cout << buff[i*col + j];
				}
				printf("\n");
			}
			printf("\n");
		}
		void fprintf(std::string filename) const
		{
			std::fstream file(filename + ".txt");
			file << row << " " << col << "\n";
			std::string* buff;
			buff = new std::string[row*col];
			int* wid = new int[col] {};
			//        F*O*R*M*A*T*T*I*N*G
			for (register int i = 0; i < row; i++)
				for (register int j = 0; j < col; j++)
				{
					std::string str = com::to_string(Mat[i][j]);
					if (str.length() > wid[j])
						wid[j] = str.length();
					buff[i*col + j] = str;
				}
			//         S*T*R*E*A*M*I*N*G
			for (register int i = 0; i < row; i++)
			{
				for (register int j = 0; j < col; j++)
				{
					pad(buff[i*col + j], wid[j] - buff[i*col + j].length() + 3);
					file << buff[i*col + j];
				}
				file << "\n";
			}
			file.close();
		}
		friend void display(const mat& m1, const mat& m2);
		friend void move(mat& sender, mat& reciever);
		void erase() { Mat = nullptr; }
		~mat()
		{
			if (Mat != nullptr)
			{
				for (register int i = 0; i < row; i++)
					delete[] Mat[i];
				delete[] Mat;
			}
		}
	};
	//int mat::constructed;
	std::ostream& operator<<(std::ostream& C, const mat& m) { m.print(); return C; }
	void display(const mat& m1, const mat& m2)
	{
		printf("\n");
		std::string *buff1, *buff2;
		buff1 = new std::string[m1.row*m1.col];
		buff2 = new std::string[m2.row*m2.col];
		int *wid1, *wid2;
		wid1 = new int[m1.col]{};
		wid2 = new int[m2.col]{};
		for (register int i = 0; i < m1.row; i++)
			for (register int j = 0; j < m1.col; j++)
			{
				std::string str = com::to_string(m1.Mat[i][j]);
				if (str.length() > wid1[j])
					wid1[j] = str.length();
				buff1[i*m1.col + j] = str;
			}
		for (register int i = 0; i < m2.row; i++)
			for (register int j = 0; j < m2.col; j++)
			{
				std::string str = com::to_string(m2.Mat[i][j]);
				if (str.length() > wid2[j])
					wid2[j] = str.length();
				buff2[i*m2.col + j] = str;
			}
		int R, r;
		r = std::fmin(m1.row, m2.row);
		R = std::fmax(m1.row, m2.row);
		for (register int i = 0; i < r; i++)
		{
			for (register int j = 0; j < m1.col; j++)
			{
				pad(buff1[i*m1.col + j], wid1[j] - buff1[i*m1.col + j].length() + 3);
				std::cout << buff1[i*m1.col + j];
			}
			printf("|  ");
			for (register int j = 0; j < m2.col; j++)
			{
				pad(buff2[i*m2.col + j], wid2[j] - buff2[i*m2.col + j].length() + 3);
				std::cout << buff2[i*m2.col + j];
			}
			printf("\n");
		}
		if (m1.row > m2.row)
			for (register int i = r; i < R; i++)
			{
				for (register int j = 0; j < m1.col; j++)
				{
					pad(buff1[i*m1.col + j], wid1[j] - buff1[i*m1.col + j].length() + 3);
					std::cout << buff1[i*m1.col + j];
				}
				printf("|\n");
			}
		else if (m2.row > m1.row)
		{
			for (register int i = r; i < R; i++)
			{
				for (register int j = 0; j < m1.col; j++)
				{
					for (register int k = 0; k < wid1[j] + 3; k++)
						printf(" ");
					printf("|  ");
				}
				for (register int j = 0; j < m2.col; j++)
				{
					pad(buff1[i*m1.col + j], wid1[j] - buff1[i*m1.col + j].length() + 3);
					std::cout << buff1[i*m1.col + j];
				}
				printf("\n");
			}
		}
		printf("\n\n");
	}
	void move(mat& sender, mat& reciever)
	{
		reciever.row = sender.row;
		reciever.col = sender.col;
		reciever.Mat = sender.Mat;
		sender.Mat = nullptr;
	}

	mat operator+(double val, const mat& m) { return mat::I(std::fmin(m.Row(), m.Col()))*val + m; }
	mat operator-(double val, const mat& m) { return mat::I(std::fmin(m.Row(), m.Col()))*val - m; }
	mat operator*(double val, const mat& m) { return m * val; }
	mat operator/(double val, const mat& m)
	{
		int row = m.Row(), col = m.Col();
		mat res(row, col);
		for (register int i = 0; i < row; i++)
			for (register int j = 0; j < col; j++)
			{
				zero_check(m.get(i, j));
				res.setelement(i, j, val / m.get(i, j));
			}
		return res;
	}
	com::vec operator*(const com::vec& v, const mat& m)
	{
		chk_eq(v.Size(), m.Row());
		int row = m.Row(), col = m.Col();
		com::vec res(col); Com sum;
		for (register int i = 0; i < col; i++)
		{
			sum = 0;
			for (register int j = 0; j < row; j++)
				sum += v[j] * m.get(j, i);
			res.set(i, sum);
		}
		return res;
	}

	Com innerProd(const com::vec& v1, const com::vec& v2)
	{
		return v1.Wsum(v2);
	}
	mat outerProd(const com::vec& v1, const com::vec& v2)
	{
		chk_eq(v1.Size(), v2.Size());
		int n = v1.Size();
		mat res(n);
		for (register int i = 0; i < n; i++)
			for (register int j = 0; j < n; j++)
				res.setelement(i, j, v1[i] * v2[j]);
		return res;
	}
	mat _set(const mat& m, Com(*F)(Com x))
	{
		int row = m.Row(), col = m.Col();
		mat res(row, col);
		for (register int i = 0; i < row; i++)
			for (register int j = 0; j < col; j++)
				res.setelement(i, j, (*F)(m.get(i, j)));
		return res;
	}
	mat abs(const mat& m) { return _set(m, com::abs); }
	mat ceil(const mat& m) { return _set(m, com::ceil); }
	mat floor(const mat& m) { return _set(m, com::floor); }
	mat ring(int n, char c = 'm')
	{
		mat res(n);
		if (c == 'm')
		{
			for (register int i = 0; i < n; i++)
				for (register int j = 0; j < n; j++)
					res.setelement(i, j, (i * j) % n);
		}
		else if (c == 'a')
		{
			for (register int i = 0; i < n; i++)
				for (register int j = 0; j < n; j++)
					res.setelement(i, j, (i + j) % n);
		}
		else
			error("Please Enter a Valid Character\nError Called from ring()");
		return res;
	}

	// SOLVING A LINEAR SYSTEM Ax=b WITH TRIANGULAR MATRICES:
	com::vec /*lower triangular*/forward_substitution(const mat& A, const com::vec& b, bool pass = false)
	{
		int n = A.Row();
		if (!pass)
		{
			if (!A.is_lowerTriangular() && !A.is_square())
				error("Can't Perform Forward Substitution\nNot a Lower-Triangular Matrix");
			chk_eq(A.Col(), b.Size());
			for (register int i = 0; i < n; i++)
				if (A.get(i, i) == 0)
				{
					printf("\n______________________________________________________________________________________________________\n");
					A.print();
					error("Zero Pivot Found");
				}
		}
		com::vec x(n); Com sum;
		for (register int i = 0; i < n; i++)
		{
			sum = 0;
			for (register int j = 0; j < i; j++)
				sum += A.get(i, j) * x[j];
			x.set(i, (b[i] - sum) / A.get(i, i));
		}
		return x;
	}
	com::vec /*upper triangular*/backward_substitution(const mat& A, const com::vec& b, bool pass = false)
	{
		int n = A.Row();
		if (!pass)
		{
			if (!A.is_upperTriangular() && !A.is_square())
				error("Can't Perform Forward Substitution\nNot an Upper-Triangular Matrix");
			chk_eq(A.Col(), b.Size());
			for (register int i = 0; i < n; i++)
				if (A.get(i, i) == 0)
				{
					printf("\n______________________________________________________________________________________________________\n");
					A.print();
					error("Zero Pivot Found\nError Call from backward_substitution()");
				}
		}
		com::vec x(n); Com sum;
		for (register int i = n - 1; i >= 0; i--)
		{
			sum = 0;
			for (register int j = n - 1; j > i; j--)
				sum += A.get(i, j) * x[j];
			x.set(i, (b[i] - sum) / A.get(i, i));
		}
		return x;
	}
	/*com::vec solveLinearSystem(const mat& A, const com::vec& b)
	{
		if (!A.is_square())
			error("Not a Square Linear System");
		mat l, u, p; A.lu(l, u, p); l = p * l;
		int n = A.Row();
		for (register int i = 0; i < n; i++)
			if (u.get(i, i) == 0 || l.get(i, i) == 0)
			{
				std::cout << l << u;
				error("A Singular Linear System");
			}
		com::vec d = forward_substitution(l, p*b, 1);
		return backward_substitution(u, d, 1);
	}*/
}

#undef fileName