/***************************************************************************
 *            springs.cc
 ****************************************************************************/

#include <cstdarg>
#include "ariadne.h"

using namespace Ariadne;

/// Control variables:
/// x1: position of the first mass
/// x2: position of the second mass
/// v1: speed of the first mass
/// v2: speed of the second mass

int main() 
{    

    /// Constants
    float m1 = 4.0; // Mass of the first ball
    float m2 = 0.75; // Mass of the second ball
    float k1 = 2.0; // Elastic constant of the first spring
    float k2 = 1.0; // Elastic constant of the second spring
    float p1 = 1.0; // Neutral position for the first spring
    float p2 = 2.0; // Neutral position for the second spring 
    float x1_0 = 0.0; // Initial position for the first spring
    float x2_0 = 3.0; // Initial position for the second spring
    float st = 1.9; // Stickyness
    float EVOL_TIME = 25.0; // Evolution time
    int   EVOL_TRANS = 4; // Evolution transitions
    float MAX_ENCLOSURE_WIDTH = 0.02; // Maximum enclosure width
    float MAX_STEP_SIZE = 0.05; // Maximum integration step size

    /// Build the Hybrid System
  
    /// Create a HybridAutomaton object
    HybridAutomaton system("springs");
  
    /// Create the discrete states
    DiscreteLocation free(1);
    DiscreteLocation stuck(2);

    /// Create the discrete events
    DiscreteEvent sticking(1);
    DiscreteEvent unsticking(2);
    
    /// Create the dynamics

    /// Free oscillation (x1'=vx1; x2'=vx2; vx1'=k1*(p1-x1)/m1; vx2'=k2*(p2-x2)/m2; t'=1 )
    Matrix<Float> A1 = Matrix<Float>(5,5);
    Vector<Float> b1 = Vector<Float>(5);
    A1[0][2] = 1.0;
    A1[1][3] = 1.0;
    A1[2][0] = -k1/m1;
    b1[2] = k1*p1/m1;
    A1[3][1] = -k2/m2;
    b1[3] = k2*p2/m2;
    b1[4] = 1.0;
    VectorAffineFunction free_d(A1,b1);
    
    /// Stuck oscillation (x1'=vx1; x2'=vx2; vx1'=vx2'=(k1*p1+k2*p2-(k1+k2)*x1)/(m1+m2); t'=1 )
    Matrix<Float> A2 = Matrix<Float>(5,5);
    Vector<Float> b2 = Vector<Float>(5);
    A2[0][2] = 1.0;
    A2[1][3] = 1.0;
    A2[2][0] = -(k1+k2)/(m1+m2);
    b2[2] = (k1*p1+k2*p2)/(m1+m2);
    A2[3][1] = -(k1+k2)/(m1+m2);
    b2[3] = (k1*p1+k2*p2)/(m1+m2);
    b2[4] = 1.0;
    VectorAffineFunction stuck_d(A2,b2);

    /// Create the resets

    /// Stick (x1^=x1; x2^=x2; vx1^=vx2^=(m1*vx1+m2*vx2)/(m1+m2) ) 
    VectorAffineFunction stick_r(Matrix<Float>(5,5, 1.0, 0.0, 0.0,        0.0,        0.0,
                                              0.0, 1.0, 0.0,        0.0,        0.0,
                                              0.0, 0.0, m1/(m1+m2), m2/(m1+m2), 0.0,
                                              0.0, 0.0, m1/(m1+m2), m2/(m1+m2), 0.0,
                                              0.0, 0.0, 0.0,        0.0,        1.0),
                           Vector<Float>(5));
    /// Unstick (do nothing)
    IdentityFunction unstick_r(5);

    /// Create the guards
    
    /// Guard for the transition between free and stuck (x1 >= x2)
    VectorAffineFunction free2stuck_g(Matrix<Float>(1,5,1.0,-1.0,0.0,0.0,0.0),Vector<Float>(1));
    /// Guard for the transition between stuck and just free ((k1-k2)*x1 + k2*p2 - k1*p1 >= st)
    VectorAffineFunction stuck2free_g(Matrix<Float>(1,5,k1-k2,0.0,0.0,0.0,0.0),Vector<Float>(1,k2*p2 -k1*p1 -st));
    //VectorAffineFunction stuck2free_g(Matrix<Float>(1,4,0.0,0.0,-1.0,1.0),Vector<Float>(1,-st)); 


    /// Build the automaton
    
    /// Locations
    system.new_mode(free,free_d);
    system.new_mode(stuck,stuck_d);
    /// Invariants
    //springs.new_invariant(free,free2stuck_g);
    /// Events
    system.new_forced_transition(sticking,free,stuck,stick_r,free2stuck_g);
    system.new_forced_transition(unsticking,stuck,free,unstick_r,stuck2free_g);

    /// Finished building the automaton

    cout << "Automaton = " << system << endl << endl;

    /// Compute the system evolution

    /// Create a HybridEvolver object
    HybridEvolver evolver(system);
    evolver.verbosity = 1;

    /// Set the evolution parameters
    evolver.settings().minimum_discretised_enclosure_widths[1] = Vector<Float>(5,MAX_ENCLOSURE_WIDTH);
    evolver.settings().minimum_discretised_enclosure_widths[2] = Vector<Float>(5,MAX_ENCLOSURE_WIDTH);
    evolver.settings().hybrid_maximum_step_size[1] = MAX_STEP_SIZE;
    evolver.settings().hybrid_maximum_step_size[2] = MAX_STEP_SIZE;
    std::cout <<  evolver.settings() << std::endl;

    // Declare the type to be used for the system evolution
    typedef HybridEvolver::EnclosureType HybridEnclosureType;
    typedef HybridEvolver::OrbitType OrbitType;

    std::cout << "Computing evolution..." << std::endl;

    Box initial_box(5, x1_0,x1_0, x2_0,x2_0, 0.0,0.0, 0.0,0.0, 0.0,0.0);

    HybridEnclosureType initial_enclosure(free,initial_box);
    
    std::cout << "Initial set=" << initial_enclosure << std::endl;
  
    HybridTime evolution_time(EVOL_TIME,EVOL_TRANS);

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    PlotHelper plotter(system.name());
    plotter.plot(orbit.reach(),"reach");
}
