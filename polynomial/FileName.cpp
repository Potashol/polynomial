#include<iostream>
#include<cmath>
#include<vector>
#include <optional>
#include <algorithm>
#include <string>
#include <stack>
#include <random>
#include <functional>
#include <cmath>
#include <iomanip>


//Задание 2

using namespace std;



class polynomial
{
private:
	int n; // cтепень многочлена
	vector<double> c; //коэфициенты
public:
	polynomial();
	polynomial(int _n, vector<double> _c);
	polynomial(int _n);
	void print(); //вывод многочлена на экран
	double D() const;
	double f(int x) const; // подстановка значения в функцию(оказалась ненужной)
	polynomial gorner(double a) const; //схема Горнера, тоже не нужна
	polynomial derivative() const;//производная
	polynomial integral() const; //первообразная
	polynomial operator *(int a) const; //умножение на число 
	polynomial operator*(polynomial other_polynomial) const;
	polynomial operator^(int a) const;
	int multiplicities();//нахождение кратности 
	int get_n();
	double get_c(int n);
	size_t get_c_size();

};

double polynomial::get_c(int n)
{
	if (n >= 0 && n <= c.size())
		return c[n];
	else
		return 0;
}

int polynomial::get_n()
{
	return n;
}

size_t polynomial::get_c_size()
{
	return c.size();
}

polynomial::polynomial(int _n, vector<double> _c)
{
	n = _n;
	c = _c;
}
polynomial::polynomial()
{
	n = 0;
	vector<double>_c(n, 1);
	c = _c;
}
polynomial::polynomial(int _n)
{
	n = _n;
	vector<double>_c(n + 1, 1);
	c = _c;
}
void polynomial::print()
{
	if (n != 0)
	{
		if (c[0] > 0 && c[0] != 1)
			cout << c[0] << "x^" << n;
		else if (c[0] != -1 && c[0] < 0)
			cout << c[0] << "x^" << n;
		else if (c[0] == -1)
			cout << "-x^" << n;
		else
			cout << "x^" << n;
		for (int i = 1; i < c.size() - 1; ++i)
		{
			if (c[i] != 1 && c[i] >= 0)
				cout << " + " << c[i] << "x^" << n - i;
			else if (abs(c[i]) != 1 && c[i] < 0)
				cout << c[i] << "x^" << n - i;
			else if (c[i] == -1)
				cout << " - x^" << n - i;
			else
				cout << " + " << "x^" << n - i;
		}
		if (c[n] >= 0)
			cout << " + " << c[n] << endl;
		else
			cout << c[n] << endl;
	}
	else
		cout << 1 << endl;
}


polynomial polynomial::operator*(int a) const
{
	vector<double> c1;
	for (int i = 0; i < c.size(); ++i)
		c1[i] *= a;
	return(polynomial(n, c1));
}

polynomial polynomial::operator*(polynomial other_polynomial) const
{
	vector<double> c1(c.size() + other_polynomial.c.size() - 1, 0);
	for (int i = 0; i <= n; ++i)
		for (int j = 0; j <= other_polynomial.n; ++j)
			c1[i + j] += c[i] * other_polynomial.c[j];
	return polynomial(n + other_polynomial.n, c1);
}

polynomial polynomial::operator^(int a) const
{
	int power = abs(a);
	polynomial p = *this;
	polynomial r(0);
	while (power > 0)
	{
		if (power % 2 == 1)
			r = r * p;
		p = p * p;
		power /= 2;
	}
	return r;
}

polynomial polynomial::derivative() const
{
	vector<double> c1(n, 1);
	int power = n - 1;
	for (int i = 0; i <= power; ++i)
	{
		c1[i] *= c[i] * (n - i);
	}
	return polynomial(power, c1);
}

polynomial polynomial::integral() const
{
	vector<double>c1(n + 2, 1);
	for (int i = 0; i <= n; ++i)
		c1[i] = c[i] / (n - i + 1);
	return polynomial(n + 1, c1);
}

double polynomial::f(int x) const
{
	double sum = 0;
	for (int i = 0; i < c.size(); ++i)
		sum += pow(x, n - i) * c[i];
	return sum;
}

polynomial polynomial::gorner(double a) const
{
	vector<double> c1(n, 0);
	c1[0] = c[0];
	for (int i = 1; i < n; ++i)
		c1[i] = a * c1[i - 1] + c[i];
	return polynomial(n - 1, c1);

}

vector<vector<double>> silvestr(polynomial p, polynomial px) // составление матрицы Сильвестра
{
	vector<vector<double>> c1(p.get_n() + px.get_n(), vector<double>(p.get_n() + px.get_n())); // размерность матрицы = степень многочлена(n) + степень его производной(m)

	for (int i = 0; i < p.get_n() + px.get_n(); ++i) // вся матрица заполняется нулями
		for (int j = 0 + i; j < p.get_n() + px.get_n(); ++j)
			c1[i][j] = 0;
	int iter = 0;
	for (int i = 0; i < px.get_n(); ++i) //заполнение m строк коэфициентами многочлена
	{
		for (int j = 0 + i; j < p.get_n() + px.get_n(); ++j)
		{
			if (j - i < p.get_n() + 1)
				c1[i][j] = p.get_c(iter);
			++iter;
		}
		iter = 0;
	}
	int it = 0;  iter = 0;
	for (int i = px.get_n(); i < p.get_n() + px.get_n(); ++i) //заполнение n строк коэфициентами производной многочлена
	{
		for (int j = 0 + it; j < p.get_n() + px.get_n(); ++j)
		{
			if (j - it < px.get_n() + 1)
				c1[i][j] = px.get_c(iter);
			++iter;
		}
		++it; iter = 0;
	}
	return c1;
}


double det(vector<vector<double>> s) //определитель матрицы(по алгоритму из лекции)
{
	double result = 1;
	size_t n = s.size();
	for (int i = 0; i < n; ++i)
	{
		auto it = std::max_element(begin(s) + i, end(s), [&i](auto& x, auto& y) {return abs(x[i]) < abs(y[i]); });
		std::iter_swap(begin(s) + i, it);
		double mii = s[i][i];
		if (abs(mii) < 10E-8)
			return 0;
		result *= mii;
		if (it - begin(s) != i)
			result = -result;
		std::transform(s[i].begin() + i, s[i].end(), s[i].begin() + i,
			[&mii](auto el) {return el / mii; });
		for (int j = i + 1; j < n; ++j)
		{
			double mji = s[j][i];
			std::transform(s[i].begin() + i, s[i].end(), s[j].begin() + i,
				s[j].begin() + i, [&mji](auto a, auto b) {return b - mji * a; });
		}
	}
	return result;
}

double polynomial::D() const //нахождение дискриминанта по формуле с результантом
{
	polynomial r = *this;
	polynomial rx = r.derivative();
	return (pow(-1, n * (n - 1) / 2) / c[0]) * det(silvestr(r, rx));
}

int polynomial::multiplicities()
{
	polynomial p = *this;
	if (p.D() != 0)
		return 0;
	else
		if (c[0] != 0)
			return n;
		else
			return n - 1;
}

int main()
{
	polynomial p(2, { 1,2,1 });
	cout << p.multiplicities();

}