#include "ariadne.h"
#include "circle.h"
#include "timer.h"
#include "skin-temperature.h"

using namespace Ariadne;

int main(int argc, char* argv[])
{
    /// Get the automata
    HybridIOAutomaton circle = getCircle();
    HybridIOAutomaton timer = getTimer();
    HybridIOAutomaton skin_temperature = getSkinTemperature();

    HybridIOAutomaton timer_circle = compose("timer-circle",timer,circle,DiscreteLocation("work"),DiscreteLocation("1"));
    HybridIOAutomaton system = compose("laser",timer_circle,skin_temperature,DiscreteLocation("work,1"),DiscreteLocation("resting"));

	RealScalarFunction guard_function = system.guard_function(DiscreteLocation("work,1,resting"),DiscreteEvent("laser_comes"));

	Box test_box(4, 0.0,0.0, 0.00,0.00, 0.161374,0.962411, 2.83028,3.06335);
	//Box test_box(4, 0.0,0.0, 0.00,0.00, 3.0,3.0, 0.0,0.0);
	TaylorSet set_model(test_box);
	cout << "Guard function: " << guard_function << endl;

	Interval boxEvaluation = guard_function.evaluate(test_box);
	cout << "Box evaluation" << boxEvaluation << endl;
	TaylorModel taylorEvaluation = guard_function.evaluate(set_model.models());
	cout << "TaylorSet evaluation: " << taylorEvaluation << endl;
}
