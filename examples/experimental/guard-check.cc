#include "ariadne.h"

using namespace Ariadne;

int main(int argc, char* argv[])
{
	Box test_box(1, -0.1,0.1);
	TaylorSet set_model(test_box);

	RealVariable x("x");

	RealExpression expr = Ariadne::sqrt(Ariadne::sqr(x));

	cout << "Expression: " << expr << endl;

    Map<ExtendedRealVariable,Interval> intervalValuation;
    intervalValuation[x]=test_box[0];

    Interval intervalEvaluation = Ariadne::evaluate(expr,intervalValuation);
    cout << "Interval evaluation: " << intervalEvaluation << endl;

    Map<ExtendedRealVariable,TaylorModel> taylorValuation;
    taylorValuation[x]=set_model[0];

	TaylorModel taylorEvaluation = Ariadne::evaluate(expr,taylorValuation);
	cout << "TaylorSet evaluation: " << taylorEvaluation << endl;
}
