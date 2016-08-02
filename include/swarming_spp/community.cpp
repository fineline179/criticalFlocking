#include "community.h"
#include "grid.h"
#include "cinder/Rand.h"
#include <math.h>
#include <vector>
#include <algorithm>

// this is normally supplied in the command line when compiling the library.
//  since I am just going from the source here, I've defined it manually
#define DIM 3
#define M_PI 3.14159265358979

/*----------------------- Community class --------------------------*/

Community::Community(int nags , double L, Agent* ags , double* p, double* v, double* vN,
                     double* agSepInfo){
    num_agents = nags ;
    agents = ags ;
    pos = p ;
    vel = v ;
    velNorm = vN;
    agentSepInfo = agSepInfo;
    box_size = L ;
    use_grid = false ;
    grid = NULL ;
    av_num_neighbors = 0.0;
}

double* Community::get_pos(){ return pos ; }

double* Community::get_vel(){ return vel ; }

Agent* Community::get_agents(){ return agents ;} 

int Community::get_num_agents(){ return num_agents ; }

double Community::get_av_num_neighbors() { return av_num_neighbors ; }

double* Community::get_AgentSepInfo() { return agentSepInfo; }

double Community::get_box_size(){ return box_size ;} 

// Initialization

void Community::randomize_positions(float width)
{
    for (int i = 0; i < DIM*num_agents; i++)
        pos[i] = (spp_frandom() - 0.5f) * width + box_size / 2;

}

void Community::randomize_directions(double v0){
    /*
     * Fill the vel list with random
     * velocities of norm v0.
     * Uses spp_random_vector from agent.cpp
     * (imported through grid.h -> agent.h)
     */
    //int i, j ;
    //double v2 ;
    //for(i=0; i<num_agents; i++){
    //    v2 = spp_random_vector(vel + i*DIM);
    //    for (j = 0; j < DIM; j++)
    //        vel[i*DIM + j] *= v0 / sqrt(v2);
    //}

    // Using cinder's random unit Vec3f instead
    int i, j;
    ci::Vec3f unitVec;
    for (i = 0; i < num_agents; i++)
    {
        unitVec = ci::Rand::randVec3f();
        for (j = 0; j < DIM; j++)
        {// random velocity
            //vel[i*DIM + j] = v0 * unitVec[j];
            // constant velocity (positive z axis)
            //vel[i*DIM + j] = v0 * ci::Vec3f(1.0, 0.0, 0.0)[j] +   (v0 / 3)*unitVec[j];
            vel[i*DIM + j] = v0 * ci::Vec3f(1.0, 0.0, 0.0)[j];
        }
    }

}

// Kinematic

void Community::move(double dt){
    for(int i=0; i<num_agents*DIM; i++)
        pos[i] += dt * vel[i] ;
}


void Community::updateAgentSepInfo()
{
    int i, j, k;
    double temp;

    for (i = 0; i < num_agents - 1; i++){
        for (j = i + 1; j < num_agents; j++){
            // make sure to accumulate sep squared on 0
            agentSepInfo[i * 4 * num_agents + j * 4 + 3] = 0.0;
            agentSepInfo[j * 4 * num_agents + i * 4 + 3] = 0.0;
            for (k = 0; k < 3; k++){
                temp = pos[DIM*j + k] - pos[DIM*i + k];
                agentSepInfo[i * 4 * num_agents + j * 4 + k] = temp;
                // symmetric entry in matrix
                agentSepInfo[j * 4 * num_agents + i * 4 + k] = -temp;
                // ij separation squared
                agentSepInfo[i * 4 * num_agents + j * 4 + 3] += temp*temp;
                agentSepInfo[j * 4 * num_agents + i * 4 + 3] += temp*temp;
            }
        }
    }
}

// Consensus protocol

void Community::sense_velocities(double* vel_sensed){
    /* 
     * If using grid, this fills the grid from scratch
     * at every iteration.
     */
    int num_neis ;
    Agent* neis ;
    if(use_grid){
        fill_grid() ;
        for(int i=0; i<num_agents; i++){
            neis = grid->get_neighborhood(agents+i , &num_neis ) ; 
            agents[i].sense_velocity(num_neis , neis , vel_sensed + i*DIM) ;
        }
    }else{
        for(int i=0; i<num_agents; i++)
            agents[i].sense_velocity(num_agents , agents , vel_sensed + i*DIM) ;
    }
}

double Community::sense_velocities_and_velsq(double* vel_sensed)
{
    double sum_of_vel_sq = 0;
    double* vel_sq = new double;
    *vel_sq = 0;
    av_num_neighbors = 0.0;
    double* num_neighbors = new double;
    /*
    * If using grid, this fills the grid from scratch
    * at every iteration.
    */
    int num_neis;
    Agent* neis;
    if (use_grid)
    {
        fill_grid();
        for (int i = 0; i<num_agents; i++)
        {
            neis = grid->get_neighborhood(agents + i, &num_neis);
            agents[i].sense_velocity(num_neis, neis, vel_sensed + i*DIM, vel_sq,
                                     num_neighbors,agentSepInfo);
            sum_of_vel_sq += *vel_sq;
            av_num_neighbors += *num_neighbors;
        }
        av_num_neighbors /= num_agents;
    }
    else
    {
        for (int i = 0; i < num_agents; i++)
        {
            agents[i].sense_velocity(num_agents, agents, vel_sensed + i*DIM, vel_sq,
                                     num_neighbors, agentSepInfo + i * 4 * num_agents);
            sum_of_vel_sq += *vel_sq;
            av_num_neighbors += *num_neighbors;
        }
        av_num_neighbors /= num_agents;
    }

    return sum_of_vel_sq;
}

void Community::update_velocities(double* vel_sensed)
{
    for (int i = 0; i < num_agents; i++)
    {
        double speed = sqrt(vel_sensed[i*DIM] * vel_sensed[i*DIM] +
                            vel_sensed[i*DIM + 1] * vel_sensed[i*DIM + 1] + 
                            vel_sensed[i*DIM + 2] * vel_sensed[i*DIM + 2]);
        for (int j = 0; j < DIM; j++)
        {
            vel[i*DIM + j] = vel_sensed[i*DIM + j];
            velNorm[i*DIM + j] = vel_sensed[i*DIM + j] / speed;
        }
    }
}

// Statistical properties

void Community::mean_position(double* meanpos){
    int i, ia ;
    for(i=0; i<DIM; i++){
        meanpos[i] = 0.0 ;
        for(ia=0; ia<num_agents; ia++){
            meanpos[i] += pos[ia*DIM + i] ;
        }
        meanpos[i] /= num_agents ;
    }
}

double Community::mean_velocity(double* meanvel){
    int i, ia ;
    double speed2 = 0. ;

    for (i = 0; i<DIM; i++)
        meanvel[i] = 0.0 ;
    for(ia=0; ia<num_agents; ia++){
        for(i=0; i<DIM; i++)
            meanvel[i] += vel[ia*DIM + i ] ;
    }
    for(i=0; i<DIM; i++){
        meanvel[i] /= num_agents ;
        speed2 += meanvel[i] * meanvel[i] ;
    }
    return sqrt(speed2);
}

double Community::polarization(double* polar)
{
    int i, ia;
    double polarizationMag = 0.;
    for (i = 0; i < DIM; i++)
        polar[i] = 0.0;
    for (ia = 0; ia < num_agents; ia++){
        for (i = 0; i < DIM; i++)
            polar[i] += velNorm[ia*DIM + i];
    }
    for (i = 0; i < DIM; i++)
    {
        polar[i] /= num_agents;
        polarizationMag += polar[i] * polar[i];
    }
    return polarizationMag;
}

double Community::order_parameter(double v0)
{
    double meanvel[DIM] ;
    double mv2 = this->mean_velocity(meanvel) ;
    return sqrt(mv2)/v0 ;
}


// TODO: fix this to not require all the square root calls (ie the distances are calculated 
//       somewhere else already)
void Community::correlation_histo(double maxRadius, int n_bins, double v0, 
                                  std::vector<double>& totalcorr, 
                                  std::vector<int>& count 
                                  /*double* totalcorr, int* count*/)
{
    /*
     * Compute the correlations in speed fluctuations in the system.
     * Mathematical formulation based on the work by Attanasi et al. in
     *      PLoS Comput Biol 10, e1003697 (2014)
     *
     * The function computes the total correlation *totalcorr* and the
     * count of pairs of agents at a certain distance *count* separate
     * instead of returning just *totalcorr/count* (Eq 2) to be able
     * to compute correctly the cumulative correlation (Eq 3) and
     * the susceptibility.
     *
     * WARNING: The normalizing factor *norm* assumes that all the agents
     * have velocity with modulus *v0*.
     *
     */
    int i, ia, ja, bin;
    double mv[DIM] ;
    double *v1, *v2 ;
    double dist, meanSpeed, speed1, speed2;
    double bindist = n_bins / maxRadius;

    for (i = 0; i < n_bins; i++)
    {
        totalcorr[i] = 0.; 
        count[i] = 0;
    }

    meanSpeed = this->mean_velocity(mv) ;
    
    for(ia=0; ia<num_agents; ia++){
        v1 = agents[ia].get_vel();
        speed1 = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
        for(ja=ia+1; ja<num_agents; ja++){
            v2 = agents[ja].get_vel() ;
            speed2 = sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
            dist = sqrt( agents[ia].distance2( agents[ja].get_pos() ) ) ;
            bin = int(dist * bindist);
            // make sure we don't write outside of memory bounds when an agent gets outside our max range
            if (bin < n_bins)
            {
                count[bin] += 1;
                // note we are using the CURRENT measured mean speed of the flock in mv, rather than 
                //  the set mean speed parameter
                totalcorr[bin] += (speed1 - meanSpeed) * (speed2 - meanSpeed);
            }
        }
    }

    for (i = 0; i < n_bins; i++)
    {
        //totalcorr[i] *= norm;
        if (count[i] != 0)
            totalcorr[i] /= count[i];
    }
}

// Other
double Community::max_distance(){
    /* 
     * The furthest distance in an N-dimensional
     * periodic hypercube is that between the center
     * of the cube and any of its vertices.
     * For a cube with sizes L_i, that is:
     *      d^2 = sum_i (L_i/2.)^2
     *  for L_i=L forall i,
     *      d^2 = L * N / 2^2
     */
    return box_size * sqrt(DIM / 4.) ;
}

// Optimization related

void Community::setup_grid(Grid *g){
    grid = g ;
    use_grid = true ;
}

void Community::fill_grid(){
    grid->fill_grid( num_agents, agents ) ;
}

/*------------------- End Community class --------------------------*/

double* spp_community_alloc_space(int num_agents){
    return new double[ num_agents * DIM ] ;
}

Agent* spp_community_alloc_agents(int num_agents){
    return new Agent[num_agents] ;
}

Agent** spp_community_alloc_neighbors(int num_agents){
    return new Agent*[num_agents] ;
    }

Agent* spp_community_build_agents(int num_agents, double* pos, double* vel, double* velNorm, Agent** neis, Behavior* behavior){
    /* 
     * WARNING: all the agents share the same
     * *behavior* and the same *neis* pointer.
     * Sharing the same *neis* pointer may screw with
     * paralellization!
     */
    Agent* ags = spp_community_alloc_agents(num_agents) ;
    for(int ia=0; ia<num_agents; ia++)
        ags[ia] = Agent(pos + ia*DIM, vel + ia*DIM, velNorm + ia*DIM, neis, behavior);
    return ags ;
}

Community spp_community_autostart(int num_agents, double speed, double box_size, Behavior* behavior){
    double* pos     = spp_community_alloc_space(num_agents) ;
    double* vel     = spp_community_alloc_space(num_agents);
    double* velNorm = spp_community_alloc_space(num_agents);
    double* agentSepInfos = new double[num_agents * num_agents * (DIM + 1)];
    std::fill(agentSepInfos, agentSepInfos + num_agents * num_agents * (DIM + 1), 0.0);
    Agent** neis    = spp_community_alloc_neighbors(num_agents) ;
    Agent* ags      = spp_community_build_agents(num_agents, pos, vel, velNorm, neis, behavior) ;
    Community com = Community(num_agents, box_size, ags, pos, vel, velNorm, agentSepInfos);

    /* Starting positions and velocities */
    com.randomize_positions(5.0) ;
    com.randomize_directions(speed) ;
    return com ;
}
