#include "ariadne.h"

using namespace Ariadne;

int main(int argc, char* argv[])
{
	Box test_box(1, -0.1,0.1);
	TaylorSet set_model(test_box);

	cout << "Original: " << set_model[0] << endl;
	TaylorModel square = sqr(set_model[0]);
    cout << "Square: " << square << endl;
    cout << "Square range: " << square.range() << endl;
    cout << "Sqrt of square" << sqrt(square) << endl;

}
