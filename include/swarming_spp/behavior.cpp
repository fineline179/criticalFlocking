#include "behavior.h"
#include "agent.h"
#include "cinder/Rand.h"
#include <random>

#define DIM 3
#define M_PI 3.14159265358979
 
/*
* Bialek fixed velocity Consensus - 1107.0604
*/
Bialek_fixedVel_consensus::Bialek_fixedVel_consensus(Interaction* ii, double vzero, 
                                                     double alph, double bet, double nois,
                                                     double r_a, double r_b, double r_e, 
                                                     double r_0)
{
    inter = ii;
    v0    = vzero;
    alpha = alph; beta = bet;
    noise = nois;
    ra = r_a; rb = r_b; re = r_e; r0 = r_0;
    ra_minus_rb = ra - rb;
    dt = 1.0;
}

void Bialek_fixedVel_consensus::sense_velocity(Agent* ag, int num_agents, Agent* ags, 
                                               double* new_vel)
{
    int i, j;
    double v2 = 0.0;
    int num_neis;
    double neiMatch_vel[DIM];
    double neiAttract_force[DIM];
    double noise_force[DIM];

    for (i = 0; i < DIM; i++)
    {
        new_vel[i] = 0.;
        neiMatch_vel[i] = 0.;
        neiAttract_force[i] = 0.;
        noise_force[i] = 0.;
    }

    // gets _location_ of neighbor list for source agent "ag"
    Agent** neis = ag->get_neis();

    // fill neighborlist for *ag* with all relevant topological neighbors.
    // there should be num_neis - 1 actual neighbors (ie we don't include the source agent) 
    num_neis = inter->get_neighbors(ag, num_agents, ags, neis);

    // calc av velocity of all neighbors
    for (j = 0; j < num_neis; j++)
    {
        for (i = 0; i < DIM; i++)
            neiMatch_vel[i] += neis[j]->get_vel()[i];
    }
    // subtract off velocity of source agent from sum
    for (i = 0; i < DIM; i++)
        neiMatch_vel[i] -= ag->get_vel()[i];
    for (i = 0; i < DIM; i++)
        neiMatch_vel[i] *= alpha;

    // calc sum of attraction of all neighbors
    Geometry* ge = inter->g;
    // cache source agent position
    double* ag_pos = ag->get_pos();
    double* nei_pos;
    double ag_nei_dist;
    double ag_nei_dir[DIM];
    double softcoreScaleFac;
    for (j = 0; j < num_neis; j++)
    {
        // cache neighbor j's position
        nei_pos = neis[j]->get_pos();
        // distance between source and neighbor j
        ag_nei_dist = sqrt( ge->distance2(ag_pos, nei_pos) );
        // skip self interaction
        if (ag_nei_dist == 0.0f)
            continue;
        // unit dir from source to neighbor j
        ge->displacement(ag_pos, nei_pos, ag_nei_dir);
        for (i = 0; i < DIM; i++)
            ag_nei_dir[i] /= ag_nei_dist;
        // scale factor as in eq. (12)
        softcoreScaleFac = .25 * (ag_nei_dist - re) / ra_minus_rb;
        // add attraction force from neighbor j to sum
        if (ag_nei_dist < rb) // hardcore repel
        {
            // ignore all other terms and return direction away from neighbor j, normalized to v0
            for (i = 0; i < DIM; i++)
                new_vel[i] = -v0 * ag_nei_dir[i];
            return;
        }
        else if (rb <= ag_nei_dist && ag_nei_dist < ra) // softcore attract
        {
            for (i = 0; i < DIM; i++)
                neiAttract_force[i] += ag_nei_dir[i] * softcoreScaleFac;
        }
        else if (ra <= ag_nei_dist && ag_nei_dist < r0) // distant attract
        {
            for (i = 0; i < DIM; i++)
                neiAttract_force[i] += ag_nei_dir[i];
        }
    }
    for (i = 0; i < DIM; i++)
        neiAttract_force[i] *= beta;

    // calc noise term
    auto unitNoise = ci::Rand::randVec3();
    for (i = 0; i < DIM; i++)
        noise_force[i] = (num_neis - 1) * noise * unitNoise[i];

    // sum of all three force terms
    for (i = 0; i < DIM; i++)
    {
        new_vel[i] = neiMatch_vel[i] + neiAttract_force[i] + noise_force[i];
        v2 += new_vel[i] * new_vel[i];
    }
    // normalize result to length v0
    for (i = 0; i < DIM; i++) new_vel[i] *= v0 / sqrt(v2);
}

void Bialek_fixedVel_consensus::rotate(double* v)
{
}


// Bialek (full) Consensus - 1307.5563
Bialek_consensus::Bialek_consensus(Interaction* ii, double vzero, double gamma, 
                                   double deltat, double j, double g, double temp,
                                   double alph, double r_a, double r_b, double r_e, double r_0,
                                   double distAttScale)
{
    inter = ii;
    v0 = vzero; gam = gamma; dt = deltat;
    J = j; G = g; T = temp;
    ra = r_a; rb = r_b; re = r_e; r0 = r_0;
    alpha = alph;
    ra_minus_rb = ra - rb;
    distAttractScale = distAttScale;
}

void Bialek_consensus::sense_velocity(Agent* ag, int num_agents, Agent* ags,
                                      double* new_vel, double* neis_vel_sq, 
                                      int* num_neighbors, double* posPairs,
                                      bool updateNeighbors)
{
    int i, j;
    int num_neis;
    double neiMatch_vel[DIM];
    double match_vel_Av[DIM];
    double neiAttract_force[DIM];
    double noise_force[DIM];

    *neis_vel_sq = 0.;
    for (i = 0; i < DIM; i++)
    {
        new_vel[i] = 0.;
        neiMatch_vel[i] = 0.;
        match_vel_Av[i] = 0.;
        neiAttract_force[i] = 0.;
        noise_force[i] = 0.;
    }

    // gets _location_ of neighbor list for source agent "ag"
    Agent** neis = ag->get_neis();

    // fill neighborlist for *ag* with all relevant balanced topological neighbors. Or use
    //  neighbors from previous call to inter->get_neighbors
    if (updateNeighbors)
        // NB: this is a hack and will break for grid use
        num_neis = inter->get_neighbors(ag, ag - ags, num_agents, ags, neis, posPairs);
    else
        num_neis = ag->get_num_neighs();

    *num_neighbors = num_neis;
    // 1) calc sum of velocity differences between agent and all neighbors
    // NB: agent is included in its neighbor list, but the term for which the neighbor is the agent
    //     itself contributes zero to the sum.

    // cache source agent velocity
    double* ag_vel = ag->get_vel();

    for (j = 0; j < num_neis; j++)
    {
        for (i = 0; i < DIM; i++)
        {
            neiMatch_vel[i] += (ag_vel[i] - neis[j]->get_vel()[i]);
            *neis_vel_sq += (ag_vel[i] - neis[j]->get_vel()[i])*(ag_vel[i] - neis[j]->get_vel()[i]);
        }
    }
    for (i = 0; i < DIM; i++)
        neiMatch_vel[i] *= -( dt * J / (/*2**/v0 * v0) );

    Geometry* ge = inter->g;
    // 2) calc difference between agent's speed and average flock speed
    double ag_vel_dir[DIM];
    double ag_speed = sqrt(ge->length2(ag_vel));
    for (i = 0; i < DIM; i++)
    {
        ag_vel_dir[i] = ag_vel[i] / ag_speed;
        match_vel_Av[i] = -( dt * G * ag_vel_dir[i] * (ag_speed - v0) / (v0*v0) );
    }

    // 3) calc sum of attraction of all neighbors
    // cache source agent position
    double* ag_pos = ag->get_pos();
    double* nei_pos;
    double ag_nei_dist;
    double ag_nei_dir[DIM];
    double softcoreScaleFac;
    for (j = 0; j < num_neis; j++)
    {
        // cache neighbor j's position
        nei_pos = neis[j]->get_pos();
        // distance between source and neighbor j
        ag_nei_dist = sqrt(ge->distance2(ag_pos, nei_pos));
        // TODO: see if we can remove this since balanced topological neighbors DON'T include 
        //       the source
        // skip self interaction
        if (ag_nei_dist == 0.0f)
            continue;
        // unit dir from source to neighbor j
        ge->displacement(ag_pos, nei_pos, ag_nei_dir);
        for (i = 0; i < DIM; i++)
            ag_nei_dir[i] /= ag_nei_dist;
        // scale factor as in eq. (12)
        softcoreScaleFac = .25 * (ag_nei_dist - re) / ra_minus_rb;
        // add attraction force from neighbor j to sum

        // NOTE: infinite repulsion force is not present in 1307.5563!
        //if (ag_nei_dist < rb) // hardcore repel
        //{
        //    // ignore all other terms and return direction away from neighbor j, normalized to v0
        //    for (i = 0; i < DIM; i++)
        //        new_vel[i] = -v0 * ag_nei_dir[i];
        //    return;
        //}
        if (/*rb <= ag_nei_dist &&*/ ag_nei_dist < ra) // softcore attract
        {
            for (i = 0; i < DIM; i++)
                neiAttract_force[i] += alpha * ag_nei_dir[i] * softcoreScaleFac;
        }
        else if (ra <= ag_nei_dist && ag_nei_dist < r0) // distant attract
        {
            for (i = 0; i < DIM; i++)
                neiAttract_force[i] += alpha * distAttractScale * ag_nei_dir[i];
        }
    }
    for (i = 0; i < DIM; i++)
        neiAttract_force[i] *= (dt / num_neis);

    // 4) calc noise term
    for (i = 0; i < DIM; i++)
        //noise_force[i] = 2 * T * ag->getNormalVariate();
        //noise_force[i] = sqrt(2 * T * dt) * ag->getNormalVariate();
        noise_force[i] = sqrt(24 * T * dt) *(ci::Rand::randFloat() - 0.5f);

    // sum of all FOUR terms
    double delta_vel[DIM];
    for (i = 0; i < DIM; i++)
    {
        new_vel[i] = ag_vel[i] + (neiMatch_vel[i] + match_vel_Av[i] + neiAttract_force[i] + noise_force[i]);
        //if (new_vel[i] > 100000.)
        //{
        //    int breakit = 1;
        //}
    }

}