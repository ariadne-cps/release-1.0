#include "ariadne.h"

using namespace Ariadne;

int main(int argc, char* argv[])
{
	Box test_box(1, -0.1,0.1);
	TaylorModel set_model = TaylorSet(test_box)[0];

	cout << "Set model: " << set_model << endl;
	TaylorModel squared = sqr(set_model);
	cout << "Squared: " << squared << endl;

	cout << "Range: " << squared.range() << endl;
}
