#include "behavior.h"
#include <stdlib.h>
#include <stdio.h>
#include <vector>

class Grid ;

class Community{
    public:

        Community()
        {
        };
        /* Construct the community class.
         * Inputs:
         *      nags = number of agents in
         *      the community.
         *      L = size of box. It is fixed
         *      to be the same in all directions
         *      for facilitating the use of a Grid.
         *      ags = pointer to the array storing
         *      the *nags* Agent instances.
         *      p = pointer to the array to store
         *      the positions of agents
         *      (size depens on DIM).
         *      v = pointer to the array to store
         *      the velocities of agents
         *      (size depens on DIM).
         */
        Community(int nags, double L, Agent* ags, double* p, double* v, double* vN,
                  double* agentSepInfos);
        /* Return the pointer to the position array.
         * The position of agent *i* corresponds to
         * the values
         *      get_pos()[i*DIM : i*DIM + DIM]
         */
        double* get_pos() ;
        /* Return the pointer to the velocity array.
         * The velocity of agent *i* corresponds to
         * the values
         *      get_pos()[i*DIM : i*DIM + DIM]
         */
        double* get_vel() ;
        /* Return the pointer to the array of agents.*/
        Agent* get_agents() ;
        /* Return how many agents are in the community.*/
        int get_num_agents() ;
        /* Return average number of neighbors (is function of opening angle) */
        double get_av_num_neighbors();
        /* Return pointer to agent sep info matrix */
        double* get_AgentSepInfo();
        /* Return box size.*/
        double get_box_size() ;
        /* Store random values in the position array *pos*.
         * The values for each dimension are in the
         * [0:box_size] range.
         */
        void randomize_positions(float width = 2.0f);
        /* Store random values in the velocity array *vel*.
         * Each random velocity has fixed norm *v0*.
         */
        void randomize_directions(double v0) ;
        /* Move all the agents during *dt* time, meaning
         * pos += dt * vel
         */
        void move(double dt) ;

        // updates the position separation vector (and its magnitude) of all agent pairs
        void updateAgentSepInfo();
        /* Sense the consensus velocities using
         * each agent's behavior and store the 
         * result in *vel_sensed*.
         * Calls Agent->behavior->sense_velocity
         * and uses a Grid if setup. If using
         * Grid, each call to this also fills the grid.
         * Note: it is not safe to use the class'
         * own *vel* as *vel_sensed*.
         */
        void sense_velocities(double* vel_sensed);
        double sense_velocities_and_velsq(double* vel_sensed);

        void update_velocities(double* vel_sensed) ;
        /* Print the position and velocity of each
         * agent to stdin. The format used is
         * x    y   vz  vz                  (2D)
         * x    y   z   vz  vy  vz          (3D)
         * Prints two blank lines at the end.
         */
        void mean_position(double* meanpos) ;
        /* Store the mean position (center of mass)
         * of all the agents in *meanpos*. Takes
         * into account the periodic boundary
         * conditions to compute a "mean angular position".
         * See https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions
         */
        double mean_velocity(double* meanvel);
        /* Store the mean velocity direction
         * of all agents in *polar*. return norm of *polar*
         */
        double polarization(double* polar);
        /* Return the order parameter of velocity alignment,
         * defined as
         *      op = sqrt{|(1/N)\sum_{i=0}^N v_i|^2} / v0 
         *  (the norm of the mean velocity divided by *v0*.)
         *  If the norm of each velocity is =< v0, then the
         *  order parameter is a number between 0 and 1.
         *  The higher the order parameter, the more aligned
         *  the agents are.
         */
        double order_parameter(double v0) ;
        /*
        * Computes the correlations in speed fluctuations in the system.
        *
        * The function computes the total correlation *totalcorr* and the
        * count of pairs of agents at a certain distance *count*.
        */
        void correlation_histo(double maxRadius, int n_bins, double v0, 
                               std::vector<double>& totalcorr,
                               std::vector<int>& count);

        /* Return the distance between the two farthest points in
        * the computation box with periodic boundary conditions.
        */
        double max_distance() ;
        /* Start using a Grid to compute
         * the neighbors of each agent.
         * The Grid instance *g* has to
         * be initialized by the user before
         * calling this function.
         */
        void setup_grid(Grid* g) ;
        /* Fill the grid calling its own
         * fill_grid method. This method
         * does NOT check that if grid has
         * been setup or not.
         */
        void fill_grid() ;
    protected:
        /* Number of agents. */
        int num_agents ;

        double av_num_neighbors;
        /* positions of the agents 
         *      Size: num_agents * DIM
         */
        double* pos ;
        /* velocities of the agents
         *      Size: num_agents * DIM
         */
        double* vel ;
        /* normal velocities of the agents
        *      Size: num_agents * DIM
        */
        double* velNorm;

        /* Array of agents
         *      Size: num_agents
         */
        Agent* agents ;
        /* Box size, same in all directions. */
        double box_size ;
        /* False by default. Turns to True
         * when a Grid instance is inseted 
         * in Community through setup_grid().
         */
        bool use_grid ;
        /* Grid instance to store the coarse
         * location of each agent. See Grid
         * documentation for more info.
         * Is it initialized to NULL and can
         * be set with setup_grid().
         */
        Grid* grid ;

        // matrix of separtion info between all agents
        double* agentSepInfo;
} ;

// Utils for automatization of the setup of a Community.

/* Allocate space for storing a vector per agent,
 * i.e. num_agents * DIM doubles, and return the
 * pointer to the array.
 */
double* spp_community_alloc_space(int num_agents) ;
/* Allocate space for a vector of Agent instances
 * and return the pointer to the array.
 */
Agent* spp_community_alloc_agents(int num_agents) ;
/* Allocate space for a vector of Agent POINTERS
 * and return the pointer to the array.
 */
Agent** spp_community_alloc_neighbors(int num_agents) ;
/* Construct an array of *num_agents* Agents where
 * each agent uses consequetive slots of pos (vel) to
 * store its position (velocity). All the agents
 * have the same behavior *behavior* and the same
 * list to store pointers to neighbors *neis*.
 *
 * WARNING: All the agents sharing the same *neis*
 * can yield to problems if paralelization is used.
 *
 * TODO add a spp_community_build_agents_parallel(...)
 */
Agent* spp_community_build_agents(int num_agents, double* pos, double* vel, Agent** neis, Behavior* behavior) ;
/* Return a Community instance i"ready to use" from scratch.
 * Allocate the required space, initialize the agents
 * and randomize their positions and velocities.
 * The speed *speed* is only used to randomize velocities
 * with a constant norm.
 *
 * Uses spp_community_alloc_neighbors and the returned Community
 * is therefore not thread-safe for parallelization.
 */
Community spp_community_autostart(int num_agents, double speed, double box_size, Behavior* behavior) ;


/* inlines */
inline int modulo(int a, int b) {
    const int result = a % b;
    return result < 0 ? result+b: result ;
}

inline double fmodulo(double a, double b) {
    const double result = fmod(a,b);
    return result < 0. ? result+b: result ;
}

