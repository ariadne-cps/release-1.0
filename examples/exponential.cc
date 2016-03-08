/*****************************************************************************************************
 *            exponential.cc
 *
 *  Copyright  2016  Luca Geretti
 *
 * Provides the behavior of a CMOS inverter fed by a sinusoidal (thus analog) input.
 *
 *****************************************************************************************************/

#include <cstdarg>
#include "ariadne.h"

using namespace Ariadne;

/// Variables:
/*
 t: absolute time;
 Vi: input voltage, consequently gate-source voltage for the nMOS, while Vdd-Vi is the gate-source voltage for the pMOS;
 Vo: output voltage, consequently drain-source voltage for the nMOS, while Vdd-Vo is the drain-source voltage for the

*/

/// Function for the behavior of the system in the nMOS linear mode, pMOS subthreshold mode (Vi >= Vth, Vo <= Vi-Vth)
/// (t' = 1; Vi' = 2*pi*f*Vdd*cos(2*pi*f*t); Vo' = -beta_n*Sn/Cl*((Vi-Vth)*Vo-Vo^2/2) + Id0/Cl*e^((-Vi-Vth+Vdd)/(nVT)) )
struct nl_pt_df : VectorFunctionData<3,3,11> {
    template<class R, class A, class P> static void
    compute(R& r, const A& x, const P& p) {
	r[0] = 1.0;
        r[1] = p[4]*2.0*pi<Real>()*p[10]*Ariadne::cos(2.0*pi<Real>()*p[10]*x[0]);
        r[2] = -p[5]*p[7]/p[3] * ((x[1]-p[1])*x[2] - x[2]*x[2]/2.0) + p[0]/p[3] * Ariadne::exp((-x[1]-p[1]+p[4])/p[2]);
    }
};

/// Function for plotting the orbit and reachability set
template<class SET> void plot(const char* filename, const int& xaxis, const int& yaxis, const int& numVariables, const Box& bbox, const Colour& fc, const SET& set, const int& MAX_GRID_DEPTH) {
    // Assigns local variables
    Figure fig;
    array<uint> xy(2,xaxis,yaxis);

    fig.set_projection_map(ProjectionFunction(xy,numVariables));
    fig.set_bounding_box(bbox);

    // If the grid must be shown
    if (MAX_GRID_DEPTH >= 0)
    {
	// The rectangle to be drawn
	Box rect = Box(numVariables);
	// Chooses the fill colour
        fig << fill_colour(Colour(1.0,1.0,1.0));

	// Gets the number of times each variable interval would be divided by 2
        int numDivisions = MAX_GRID_DEPTH / numVariables;
	// Gets the step in the x direction, by 1/2^(numDivisions+h), where h is 1 if the step is to be further divided by 2, 0 otherwise
	double step_x = 1.0/(1 << (numDivisions + ((MAX_GRID_DEPTH - numDivisions*numVariables > xaxis) ? 1 : 0)));
	// Initiates the x position to the bounding box left bound
        double pos_x = bbox[0].lower();
        // Sets the rectangle 2-nd interval to the corresponding bounding box interval (while the >2 intervals are kept at [0,0])
	rect[yaxis] = bbox[1];
        // While between the interval
        while (pos_x < bbox[0].upper())
        {
	    rect[xaxis] = Interval(pos_x,pos_x+step_x); // Sets the rectangle x coordinate
	    pos_x += step_x; // Shifts the x position
	    fig << rect; // Appends the rectangle
        }

	// Repeats for the rectangles in the y direction
	double step_y = 1.0/(1 << (numDivisions + ((MAX_GRID_DEPTH - numDivisions*numVariables > yaxis) ? 1 : 0)));
        double pos_y = bbox[1].lower();
	rect[xaxis] = bbox[0];
        while (pos_y < bbox[1].upper())
        {
	    rect[yaxis] = Interval(pos_y,pos_y+step_y);
   	    fig << rect;
	    pos_y += step_y;
        }
    }
    // Draws and creates file
    fig.set_fill_colour(fc);
    fig << set;
    fig.write(filename);
}


int main(int argc, char* argv[])
{
    /// Dynamics parameters
    Vector<Float> dp(11);

    dp[0] = 1e-6; /// Subthreshold current, Id0
    dp[1] = 0.15; /// Threshold voltage, Vth
    dp[2] = 0.035; /// n-Thermal voltage, n*VT
    dp[3] = 1e-3; /// Load capacitance, Cl
    dp[4] = 1.0; /// Operating voltage, Vdd
    dp[5] = 1.0e-2; /// Beta of the nMOS, beta_n = mu_n * Cox
    dp[6] = 0.5e-2; /// Beta of the pMOS, beta_p = mu_p * Cox
    dp[7] = 1.0; /// Shape factor of the nMOS, Sn
    dp[8] = 2.0; /// Shape factor of the pMOS, Sp
    dp[9] = 0.01; /// Early effect constant, lambda
    dp[10] = 0.25; /// Frequency, f

    /// Constants
    float EVOL_TIME = 2.0/dp[10];   /// Evolution time
    int EVOL_TRANS = 22;            /// Evolution transitions
    float MAX_ENCL_WIDTH = 0.1;   /// Maximum enclosure width
    float MAX_STEP_SIZE = 1e-3;     /// Maximum step size
    int VERBOSITY = 1;              /// Verbosity of the HybridEvolver
	if (argc > 1)
		VERBOSITY = atoi(argv[1]);

    /// Build the Hybrid System

    /// Create a HybridAutomaton object
    HybridIOAutomaton inverter;

    /// Create the discrete states
    DiscreteLocation work("work");

    RealVariable t("t");
    RealVariable vi("vi");
    RealVariable vo("vo");

    inverter.add_output_var(t);
    inverter.add_output_var(vi);
    inverter.add_output_var(vo);

	inverter.new_mode(work);

	RealExpression dyn_t = 1.0;
	inverter.set_dynamics(work, t, dyn_t);
	RealExpression dyn_vi = dp[4]*2.0*pi<Real>()*dp[10]*Ariadne::cos(2.0*pi<Real>()*dp[10]*t);
	inverter.set_dynamics(work, vi, dyn_vi);
	RealExpression dyn_vo = -dp[5]*dp[7]/dp[3] * ((vi-dp[1])*vo - vo*vo/2.0) + dp[0]/dp[3] * Ariadne::exp((-vi-dp[1]+dp[4])/dp[2]);
	inverter.set_dynamics(work, vo, dyn_vo);

    /// Compute the system evolution

    /// Create a HybridEvolver object
    HybridEvolver evolver(inverter);
    evolver.verbosity = VERBOSITY;

    evolver.settings().hybrid_maximum_step_size[work] = MAX_STEP_SIZE;
    evolver.settings().minimum_discretised_enclosure_widths[work] = Vector<Float>(3,MAX_ENCL_WIDTH);

    // Declare the type to be used for the system evolution
    typedef HybridEvolver::EnclosureType HybridEnclosureType;
    typedef HybridEvolver::OrbitType OrbitType;

    Box initial_box(3, 0.161699,0.161699, 0.850000,0.850000, 0.246241,0.246241);
    HybridEnclosureType initial_enclosure(work,initial_box);

    HybridTime evolution_time(EVOL_TIME,EVOL_TRANS);

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    Box graphic_box_vo(2, 0.0, EVOL_TIME, -2*dp[4], 2*dp[4]);
    Box graphic_box_vi(2, 0.0, EVOL_TIME, -2*dp[4], 2*dp[4]);

    plot("exponential_t_vi", 0, 1, 3, graphic_box_vi, Colour(0.0,0.5,1.0), orbit, -1);
    plot("exponential_t_vo", 0, 2, 3, graphic_box_vo, Colour(0.0,0.5,1.0), orbit, -1);
}
