#include "math.h"
#include "agent.h"
#include "behavior.h"
#include <stdlib.h>
#include <time.h>


#define DIM 3

void spp_set_seed(long int seed){
    //srandom(seed);
    // WINDOWS substitution
    srand(seed);
}

double spp_frandom(){
    //return random()*1. / RAND_MAX;
    // WINDOWS substitution
    return rand()*1. / RAND_MAX;
}

double spp_random_vector(double* v){
    /* Note that the higher the dimension the more amount
     * of 'trial' vectors will have to be sampled, as the
     * proportion between the volume of a DIM-sphere
     * and a DIM-cube decreases with increasing DIM.
     */
    double norm2 = 2.0 ; // or any value >1 to enter the loop
    int i ;
    while(norm2 > 1.0){
        norm2 = 0.0 ;
        for(i=0; i<DIM ; i++){
            v[i] = (spp_frandom()-0.5)*2.0 ;
            norm2 += v[i]*v[i] ;
        }
    }
    return norm2 ;
}

Agent::Agent(double* p , double* v, double* vN, Agent** ns, Behavior* bb){
    pos = p ;
    vel = v ;
    velNorm = vN;
    neis = ns ;
    num_neighs = 0;
    beh = bb ;

    spp_set_seed(time(NULL));

    rng = std::mt19937(std::chrono::system_clock::now().time_since_epoch().count());
    setNormalParams(0.0, sqrt(beh->getdt()) );
}

// Kinetic stuff

void Agent::move(double dt){
    for(int i=0; i<DIM ; i++)
        pos[i] += vel[i]*dt ;
}

void Agent::update_vel(double* new_vel){
    /* copy values, different from vel = new_vel */
    for(int i=0; i<DIM ; i++)
        vel[i] = new_vel[i] ;
}

double* Agent::get_pos(){
    return pos ;
}

double* Agent::get_vel(){
    return vel ;
}

double* Agent::get_velNorm(){
    return velNorm;
}

Agent** Agent::get_neis()
{
    return neis ;
}

void Agent::copy(Agent* ag){
    /* Deep copy the values of
     * ag into self.
     * TODO is *beh = *(ag->beh) always safe?
     * What if it copies from and to the
     * same memory space?
     */
    for(int i=0; i<DIM ; i++){
        pos[i] = ag->get_pos()[i] ;
        vel[i] = ag->get_vel()[i] ;
    }
    neis = ag->get_neis() ;
    *beh = *(ag->beh) ;
}

// Consensus protocol

double Agent::distance2(double *point){
    return beh->inter->g->distance2( this->pos , point) ;
}
    
int Agent::is_neighbor(Agent* nei){
    return beh->inter->is_neighbor( this, nei) ;
}

int Agent::get_neighbors(int n_agents, Agent* ags){
    return beh->inter->get_neighbors(this, n_agents, ags, neis);
}

Agent** Agent::get_neighbor_list()
{
    return neis;
}

int Agent::get_num_neighs()
{
    return num_neighs;
}

void Agent::look_around(int n_agents, Agent* ags){
    beh->inter->look_around( this, n_agents, ags);
}

void Agent::sense_velocity(int num_agents, Agent* ags, double* new_vel){
    beh->sense_velocity(this, num_agents, ags, new_vel) ; 
}

void Agent::sense_velocity(int num_agents, Agent* ags, double* new_vel, double* neis_vel_sq){
    beh->sense_velocity(this, num_agents, ags, new_vel, neis_vel_sq);
}

void Agent::sense_velocity(int num_agents, Agent* ags, double* new_vel, double* neis_vel_sq,
                           int* num_neighbors, double* posPairs)
{
    beh->sense_velocity(this, num_agents, ags, new_vel, neis_vel_sq, num_neighbors, posPairs);
    num_neighs = *num_neighbors;
}

int Agent::sense_danger(int num_threats, Agent* threats, double* new_vel){
    return beh->sense_danger(this, num_threats, threats, new_vel) ;
}
        
int Agent::sense_victims(int num_agents, Agent* ags){
    return beh->sense_victims(this, num_agents, ags) ;
}

int Agent::hunt(Agent* prey, double deltat){
    return beh->hunt(this, prey, deltat) ;
}

void Agent::setNormalParams(double mean, double std)
{
    nd = std::normal_distribution<double>(mean, std);
}

