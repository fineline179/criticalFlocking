#include <random>
#include <chrono>

class Behavior ; 

/*
 * Class describing one self-propagating agent
 * perfoming multi-agent consensus.
 * Mostly a placeholder for ease of use, the
 * algorithms for the consesus protocol are
 * defined by the Behavior class.
 *
 * This class can be used alone to simulate
 * multi-agent consensus, but for most uses
 * the class Community is a useful wrapper
 * to control a large amount of agents.
 *
 *      pos =   position of agent.
 *      vel =   velocity of agent.
 *      neis =  space to allocate a list
 *              of neighbors of the agent.
 *      beh =   Behavior instance containing
 *              the consensus protocol.
 *
 * WARNING: Expanding on Agent by creating a
 * class that inherits from it can lead to
 * errors if the size of the new class is 
 * different from Agent. The reason is that 
 * Agent uses a **Agent to store the neighbors,
 * and an array of pointers don't work very
 * well with inheriting classes.
 * For that reason the behavior of Agent is
 * almost exclusively determined by its
 * Behavior instance. Creating new behaviors
 * that inherit from Behavior is perfectly
 * fine.
 *
 */
class Agent {
    public:
        Agent() {} ;
        Agent(double* pos, double* vel, double* velNorm, Agent** nn, Behavior* bb) ;
        /* Update the agent position according to its velocity,
         *  pos += dt * vel
         */
        void move( double dt) ;
        /* Set the velocity of the agent to *new_vel*
         * (copy the values, not the pointer)
         */
        void update_vel(double* new_vel) ;
        /* Return position */
        double* get_pos() ;
        /* Return velocity */
        double* get_vel();
        double* get_velNorm();
        /* Return array of neighbors */
        Agent** get_neis() ;
        /* Turn the agent into a copy
         * of *ag*. The agent and *ag*
         * are still independen
         * The behavior instance is copied.
         * TODO: Check if that is always safe.
         * If both agent share the same behavior
         * that copies to and from the same
         * memory space.
         */
        void copy(Agent* ag) ;
        /* Get the distance from agent to *point*
         */
        double distance2(double* point) ;
        /* Return 1 if *nei* is a neighbor
         * of the agent.
         */
        int is_neighbor(Agent* nei) ;
        /* Call the beh->inter->look_around()
         * method for non-local interactions.
         */
        void look_around(int n_agents, Agent* ags) ;
        /* Copy A POINTER to all the
         * agents in *ags* that are neighbors
         * of the agent into *neis*.
         */
        int get_neighbors(int n_agents, Agent* ags);

        // returns pointer to list of neighbors. Note only the first
        Agent** get_neighbor_list();
        // returns current number of neighbors of agent (used in balanced topological)
        int get_num_neighs();
        /* Store the sensed velocity in *new_vel*.
         * Defined by the Behavior *beh*.
         */
        void sense_velocity(int num_agents, Agent* ags, double* new_vel);
        void sense_velocity(int num_agents, Agent* ags, double* new_vel, double* neis_vel_sq);
        void sense_velocity(int num_agents, Agent* ags, double* new_vel, double* neis_vel_sq,
                            int* num_neighbors, double* posPairs, bool updateNeighbors);
        /* Same as sense_velocity but using
         * the sense_danger function from *beh*.
         */
        int sense_danger(int num_threats, Agent* threats, double* new_vel) ;
        /* Call the sense_victims function in *beh*. */
        int sense_victims(int num_agents, Agent* ags) ;
        /* Call the hunt function in *beh*. */
        int hunt(Agent* prey, double deltat) ;

        void setNormalParams(double mean, double std);
        inline double getNormalVariate() { return nd(rng); };
    protected:
        /* position of the agent */
        double* pos ;
        /* velocity of the agent */
        double* vel ;
        /* direction of velocity of the agent */
        double* velNorm;
        /* neighbors of the agent */
        Agent** neis ;
        int     num_neighs;
        /* behavior (consensus protocol)
         * of the agent. */
        Behavior* beh ;

        std::mt19937 rng;
        std::normal_distribution<double> nd;
} ;

// Interface to generate random numbers in the spp library.

/* Initialize the seed for the random number generator. */
void spp_set_seed(long int seed) ;

/* Return a double random number in [0,1] interval.
 * Statistical properties of the randoms not 
 * throroughly tested.
 */
double spp_frandom() ; 
 
/* Store a random vector with norm <= 1.
 * Return the norm2 of that vector.
 * Calls spp_frandom.
 */
double spp_random_vector(double* v) ;
