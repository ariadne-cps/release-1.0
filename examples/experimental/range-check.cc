#include "ariadne.h"

using namespace Ariadne;

int main(int argc, char* argv[])
{
	Box test_box(1,1e-14,2.0);
	TaylorModel set_model = TaylorSet(test_box)[0];

	Interval i1(-2.0,1.0);
	cout << "Square of [-2,1] = " << sqr(i1) << endl;

	TaylorModel test_model1 = TaylorModel::scaling(2, 0, Interval(-2.0,1.0));
	cout << "Test model 1: " << test_model1 << endl;
	TaylorModel result = sqr(test_model1);

	cout << "Result: " << result << endl;
	cout << "Result range: " << result.range() << endl;
/*
	cout << "Set model: " << set_model << endl;
	cout << "Set model range: " << set_model.range() << endl;
	TaylorModel squared = sqr(set_model);
	cout << "Squared: " << squared << endl;
	cout << "Squared Range: " << squared.range() << endl;

	TaylorModel sqrt_sqr = sqrt(set_model);

	cout << "Sqrt of squared: " << sqrt_sqr << endl;
	*/
}
