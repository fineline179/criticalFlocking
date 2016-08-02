#include "interaction.h"
#include "agent.h"
#include <list>
#include <tuple>

using namespace std;

#define DIM 3

/*
 * Geometry
 */
double Geometry::length2(double* vect){
    double l2 = 0. ;
    int i ;
    for(i=0 ; i<DIM ; i++)
        l2 += vect[i]*vect[i] ;
    return l2 ;
}


/*
 * Cartesian
 */
Cartesian::Cartesian(double l ){
    L = l ;
}

void Cartesian::displacement(double* x0, double* x1, double* dis){
    for(int i=0 ; i<DIM ; i++)
       dis[i] = x1[i] - x0[i] ;
}

double Cartesian::distance2(double* x0, double* x1){
    double dis = 0. ;
    for(int i=0 ; i<DIM ; i++)
       dis += (x1[i] - x0[i])*(x1[i] - x0[i]) ;
   return dis ;
}


/*
 * Cartesian Periodic
 */
CartesianPeriodic::CartesianPeriodic(double l ){
    L = l ;
}

void CartesianPeriodic::displacement(double* x0, double* x1, double* dis){
    for(int i=0 ; i<DIM ; i++)
        dis[i] = (x1[i] - x0[i]) - rint( (x1[i] - x0[i])/L ) * L ; 
}

double CartesianPeriodic::distance2(double* x0, double* x1){
    double dis = 0. ;
    double tmp ;
    for(int i=0 ; i<DIM ; i++){
        tmp = (x1[i] - x0[i]) - rint( (x1[i] - x0[i])/L ) * L ; 
        dis += tmp * tmp ;
    }
   return dis ;
}


/*
 * Interactions
 */

/*
 * Metric
 */
Metric::Metric(double r , Geometry* gg){
    rad2 = r*r ;
    g = gg ;
}

int Metric::is_neighbor(Agent* a0 , Agent* a1){
    if(g->distance2( a0->get_pos() , a1->get_pos()) <= rad2)
        return 1 ;
    return 0 ;
}

int Metric::get_neighbors(Agent* a0, int n_agents, Agent* ags, Agent** neis){
    int ia ;
    int n_neis = 0 ;
    double* pos = a0->get_pos() ;
    for(ia=0; ia < n_agents ; ia++){
        if(g->distance2( pos , (ags+ia)->get_pos()) <= rad2){
            neis[n_neis] = ags + ia ;
            n_neis += 1 ;
        }
    }
    return n_neis ;
}


/*
 * Topologic
 */
Topologic::Topologic(int kk, Geometry* gg, double* dd){
    k = kk;
    g = gg ;
    rad2 = 0.0 ;
    dists2 = dd ;
}

int Topologic::is_neighbor(Agent* a0 , Agent* a1){
    if(g->distance2( a0->get_pos() , a1->get_pos()) <= rad2)
        return 1 ;
    return 0 ;
}

int Topologic::get_neighbors(Agent* a0, int n_agents, Agent* ags, Agent** neis){
    int ia;
    int n_neis = 0;
    double* pos = a0->get_pos();
    /* determine the effective radius */
    for (ia = 0; ia < n_agents; ia++)
        dists2[ia] = g->distance2( pos , (ags+ia)->get_pos()) ;
    rad2 = quickselect(dists2, n_agents, k ) ;
    // TODO: add break once we have gotten to k neighbors (at least when not using grid)
    for(ia=0; ia < n_agents ; ia++){
        if(g->distance2( pos , (ags+ia)->get_pos()) <= rad2){
            neis[n_neis] = ags + ia;
            n_neis += 1 ;
        }
    }
    return n_neis;
}

void Topologic::look_around(Agent* a0, int n_agents, Agent* ags){
    int ia;
    double* pos = a0->get_pos() ;
    /* determine the effective radius */
    for(ia=0; ia < n_agents ; ia++)
        dists2[ia] = g->distance2( pos , (ags+ia)->get_pos()) ;
    rad2 = quickselect(dists2, n_agents, k ) ;
}

/*
* Balanced Topologic
*/
TopoBalanced::TopoBalanced(int kk, double mumu, Geometry* gg, double* dd)
{
    k = kk;
    mu = mumu;
    g = gg;
    dists = dd;
}


int TopoBalanced::get_neighbors(Agent* a0, int a0_num, int n_agents, Agent* ags, Agent** neis,
                                double* pos_pairs)
{
    typedef tuple<Agent*, double*> agAndPos;
    int ia;
    double pairDotProb, denom;

    list<agAndPos> nei_pairs;
    // add pointers into this agent's slice of pos_pair struct, for all its neighbor pairs.
    //  skip agent's position from itself.
    for (ia = 0; ia < a0_num; ia++)
        nei_pairs.push_back(agAndPos(ags + ia, pos_pairs + ia * 4));
    for (ia = a0_num + 1; ia < n_agents; ia++)
        nei_pairs.push_back(agAndPos(ags + ia, pos_pairs + ia * 4));

    // cull neighbors by angles from source agent
    list<agAndPos>::iterator it1, it2;
    for (it1 = nei_pairs.begin(); it1 != nei_pairs.end();){
        for (it2 = next(it1); it2 != nei_pairs.end();){
            pairDotProb = get<1>(*it1)[0]*get<1>(*it2)[0] + get<1>(*it1)[1]*get<1>(*it2)[1] + get<1>(*it1)[2]*get<1>(*it2)[2];

            // if neighbors opening angle >= 90 degrees, both are fine
            if (pairDotProb <= 0.){
                it2++;
                continue;
            }
            // otherwise compute cos^2 between them
            pairDotProb *= pairDotProb;
            denom = 1 / (get<1>(*it1)[3] * get<1>(*it2)[3]);
            pairDotProb *= denom;
            // if we're above critical angle, both are fine
            if (pairDotProb < .75 /* cos(30 degrees)^2, hardcoded*/)
            {
                it2++;
                continue;
            }
            // otherwise we're below critical angle, so delete whichever neighbor is further away
            else
            {
                if (get<1>(*it1)[3] >= get<1>(*it2)[3]){
                    it1 = nei_pairs.erase(it1);
                    // start over, comparing the next it1 with the it2 immediately after it.
                    //it2 = nei_pairs.end();
                    it2 = next(it1);
                }
                else it2 = nei_pairs.erase(it2);
            }
        }
        it1++;
    }

    for (ia = 0, it1 = nei_pairs.begin(); it1 != nei_pairs.end(); ia++, it1++)
    {
        neis[ia] = get<0>(*it1);
    }
    return nei_pairs.size();
}

int TopoBalanced::is_neighbor(Agent* a0, Agent* a1)
{
    return 0;
}

void TopoBalanced::look_around(Agent* a0, int n_agents, Agent* ags)
{

}


/*
 * No interaction
 */
NoInteraction::NoInteraction(Geometry* gg){
    g = gg ;
}

int NoInteraction::is_neighbor(Agent* a0 , Agent* a1){
    /* Beware of incosistency, a0 is not neighbor of itself
     * but it is returned in get_neighbors
     */
    return 0 ;
}

int NoInteraction::get_neighbors(Agent* a0, int n_agents, Agent* ags, Agent** neis){
    neis[0] = a0 ;
    return 1 ;
}

/*
 * Other stuff
 */

#define SWAP(a,b) { temp=(a);(a)=(b);(b)=temp; }

double quickselect(double *arr, int n, int k) {
    int i,ir,j,l,mid;
    double a,temp;

    l=0;
    ir=n-1;
    for(;;) {
        if (ir <= l+1) { 
            if (ir == l+1 && arr[ir] < arr[l]) {
	            SWAP(arr[l],arr[ir]);
            }
            return arr[k];
        } else {
            mid=(l+ir) >> 1; 
            SWAP(arr[mid],arr[l+1]);
            if (arr[l] > arr[ir])
                SWAP(arr[l],arr[ir]) ; 
            if (arr[l+1] > arr[ir])
                SWAP(arr[l+1],arr[ir]) ; 
            if (arr[l] > arr[l+1])
                SWAP(arr[l],arr[l+1]) ; 
            i=l+1; 
            j=ir;
            a=arr[l+1]; 
            for (;;) { 
	            do i++; while (arr[i] < a); 
	            do j--; while (arr[j] > a); 
	            if (j < i) break; 
	            SWAP(arr[i],arr[j]);
            } 
            arr[l+1]=arr[j]; 
            arr[j]=a;
            if (j >= k) ir=j-1; 
            if (j <= k) l=i;
        }
    }
}
