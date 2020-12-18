#include "LinearAlgebra.h"
#include <iomanip>
#include <random>
#include <cmath>

typedef long long ll;

struct Data
{
	double a, b;
};

using namespace nm;
using namespace std;

#define f1(x) (1 / (0.35 * pow(cos(0.6 * x), 2) + pow(cos(0.6 * x), 2)))
#define f2(x) ((x * sqrt(0.7 * x * x + 1.2 * x + 2.1))/pow(2.1 + x, 2))
#define EPS 0.00005

double simpson_integration_1(double a, double b, int N)
{
	double h = (b - a) / (2 * N);
	double s = f1(a) + f1(b);
	for (int i = 1; i <= 2 * N - 1; ++i)
	{
		double x = a + h * i;
		s += f1(x) * ((i & 1) ? 4 : 2);
	}
	s *= h / 3;
	return s;
}
double simpson_integration_2(double a, double b, int N)
{
	double h = (b - a) / (2 * N);/// already is dived by two
	double s = f2(a) + f2(b);
	for (int i = 1; i <= 2 * N - 1; ++i)
	{
		double x = a + h * i;
		s += f2(x) * ((i & 1) ? 4 : 2);
	}
	s *= h / 3;
	return s;
}

int main()
{
	try
	{
		/** integral 1 **/
		vector<Data> myData  {{0.1,1.1}, {0.1,2}, {-1,0.2}, {2,2.3}, {0.1,2}};
		for(auto& item: myData)
		{
			ll N = 10;
			double a = item.a, b = item.b;
			double pr_res, cur_res = 0;
			do
			{
				pr_res = cur_res;
				N <<= 1;
				cur_res = simpson_integration_1(a, b, N);
			} while (abs(cur_res - pr_res) > EPS);

			cout << "start = " << a << " end = " << b << '\n';
			cout << "Res Integral = " << cur_res << "\n";
			cout << "For integral 1 specified precision " << EPS << " is reached with N = " << N << '\n';
			cout << '\n';
		}
	}
	catch (std::exception& exception)
	{
		std::cerr << exception.what();
	}
}
/*


 */








