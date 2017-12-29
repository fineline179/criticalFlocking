#include "cinder/app/App.h"
#include "cinder/app/RendererGl.h"
#include "cinder/Vector.h"
#include "cinder/Utilities.h"
#include "cinder/params/Params.h"
#include "cinder/Camera.h"
#include "cinder/gl/gl.h"
#include "ParticleController.h"
#include "swarming_spp/community.h"
#include "cinder/Rand.h"
#include "swarming_spp/grid.h"
#include <math.h>

using namespace ci;
using namespace ci::app;

class FlockingApp : public App {
 public:
	void keyDown( KeyEvent event );
	void setup();
	void update();
	void draw();
    void drawGrid(float boxSize, float cellSpacing, float gridRadius);
    void drawC_sp_Graph();
    void drawBirdToBirdArrow(int fromIndex, int toIndex);
    void drawBirdToNeighborsArrows(int fromIndex);

    bool     mPaused;
    bool     mShowGrid;
    float    birdFrom;
    float    birdTo;
    float    birdDimFactor;

	// PARAMS
	params::InterfaceGlRef	mParams;
	
	// CAMERA
	CameraPersp			mCam;
	quat				mSceneRotation;
	vec3				mEye, mCenter, mUp;
	float				mCameraDistance;
	
	ParticleController	mParticleController;
	float				mZoneRadius;
	bool				mCentralGravity;
	bool				mFlatten;

    //// SWARMING_SPP
    Cartesian*          g;
    TopoBalanced*       interaction;
    Bialek_consensus*   behavior;
    Grid*               grid;
    Community           com;

    // Size of bounding box of simulation (note: possibly not currently used)
    double mBoxSize;
    // Number of agents
    int mN;
    // Time step for iteration of Langevin stochastic diff eq
    double mDt;
    // Average speed of flock
    double mV0;
    // Parameters of MaxEnt model
    double mJ, mG;
    // Temperature of noise in Langevin eq
    double mTemp; 
    // Number of nearest neighbors
    float  m_nc;
    // Minimum angle of separation of nearest neighbors
    double  mBalanceAngle;
    // Characteristic radii of attraction force
    double ra, rb, re, r0;

    // space to store temporary data for SWARMING_SPP
    double* dist2;
    double* v2;

    std::vector<Particle>  spp_particles;

    int mRadialCorMaxRange;
    int mNumRadialBins;
    // update neighbor list every this number of frames
    int mNeighborUpdateFrequency;
    int mFrameCounter;
    // Correlation of speed of bird pairs as func of separation
    std::vector<double> mC_sp_r;
    // Counts of bird pairs as a function of separation (used in calc of previous quantity
    std::vector<int> mC_sp_r_count;
    double* mFlockCenter;
    double* mFlockAvVel;
    double  mFlockAvSpeed;
    double  mAvNumNeighbors;
    double* mFlockPolarization;
    double  mFlockPolarization_mag;

    // similarity of velocity of nearest neighbors in flock. Eq. (1) in 1307.5563
    double mQ_int;
};


void FlockingApp::setup()
{	
    // SETUP STATISTICAL VARIABLES
    mNumRadialBins = 20;
    mRadialCorMaxRange = 25;
    mNeighborUpdateFrequency = 1;
    // set frameCounter so neighbors update on first frame
    mFrameCounter = mNeighborUpdateFrequency;
    mC_sp_r = std::vector<double>(20);
    mC_sp_r_count = std::vector<int>(20);
    mFlockCenter = new double[3];
    mFlockAvVel = new double[3];
    mAvNumNeighbors = 0.0;
    mFlockPolarization = new double[3];

    birdFrom        = -1.0;
    birdTo          = 2.;
    mPaused         = true;
    mShowGrid       = true;
    birdDimFactor   = 0.4f;

    // SETUP SIMULATION PARAMETERS
    mBoxSize = 500.;
    mN = 300;
    //mJ = 10.0; mG = 0.08;
    mJ = 19.0; mG = 0.2;
    mDt = .03;
    mV0 = 12.57;
    mTemp = 0.55;
    m_nc = 8.;
    mBalanceAngle = 51.0;
    ra = 0.8; rb = 0.2; re = 0.5; 
    // make this large for unlimited small attraction range on non-balanced topological interactions
    r0 = 1000.;

	// SETUP CAMERA
    mCameraDistance = 20.0f;
    mEye =    vec3(0.0, 0.0, mCameraDistance);
    mCenter = vec3(0.0, 0.0, 0.0);
    mUp =     vec3(0.0, 1.0, 0.0);
	mCam.setPerspective( 75.0f, getWindowAspectRatio(), 5.0f, 2000.0f );

	// SETUP UI PARAMS
    mParams = params::InterfaceGl::create("Flocking", ivec2(200, 400));
    mParams->addParam("Scene Rotation", &mSceneRotation, "opened=1");
    mParams->addSeparator();
    mParams->addParam("J", &mJ, "min=1.0 max=60.0 step=1.0 keyIncr=f keyDecr=d");
    mParams->addParam("g",  &mG,  "min=0.01 max=0.5 step=0.01 keyIncr=v keyDecr=c");
    mParams->addParam("Temp", &mTemp, "min=0.025 max=2.0 step=0.025 keyIncr=y keyDecr=t");
    mParams->addParam("Balance Angle", &mBalanceAngle, 
                      "min=5.0 max=180.0 step=1.0 keyIncr=o keyDecr=i");
    mParams->addSeparator();
    mParams->addParam("Eye Distance", &mCameraDistance,
                      "min=5.0 max=1500.0 step=1.0 keyIncr=s keyDecr=w" );
    // which bird to highlight (1-indexed). If 0, no highlighting.
    // TODO: unhardcode max and replace with NUM_INITIAL_PARTICLES
    mParams->addParam("Bird From", &birdFrom,
                      "min=-1.0 max=9.0 step=1.0 keyIncr=. keyDecr=,");
    mParams->addParam("Bird To", &birdTo,
                      "min=0.0 max=9.0 step=1.0 keyIncr=l keyDecr=k");
    mParams->addParam("Paused", &mPaused, "keyIncr=space");
    mParams->addParam("Draw grid", &mShowGrid, "keyIncr=g");
    mParams->addSeparator();
    mParams->addParam("Mean speed", &mFlockAvSpeed, true);
    mParams->addParam("Av Num Neighs", &mAvNumNeighbors, true);
    mParams->addParam("Polarization", &mFlockPolarization_mag, true);
    mParams->addParam("Q_int", &mQ_int, true);

    // CREATE SWARMING_SPP containers
    dist2 = spp_community_alloc_space(mN);
    v2    = spp_community_alloc_space(mN);

    g = new Cartesian(mBoxSize);
    //interaction = new Topologic((int) m_nc, g, dist2);
    interaction = new TopoBalanced((int) m_nc, mBalanceAngle, g, dist2);                       
    behavior = new Bialek_consensus(interaction, mV0, 1.0,
                                    mDt, mJ, mG, mTemp,
                                    0.95, ra, rb, re, r0,
                                    1.0);
    com = spp_community_autostart(mN, mV0, mBoxSize, behavior);

    // initialize agent separation matrix
    com.updateAgentSepInfo();

    // setup SWARMING_SPP grid for faster execution
    //int nSlots = (int) sqrt( (0.5*NUM_INITIAL_PARTICLES) / mNc);
    //grid        = new Grid(nSlots, BOX_SIZE, NUM_INITIAL_PARTICLES);
    //com.setup_grid(grid);
    
    // create Cinder particles
    for (int i = 0; i < mN; i++)
        spp_particles.push_back(Particle(vec3(), vec3()));
}


void FlockingApp::keyDown( KeyEvent event )
{
}


void FlockingApp::update()
{
    gl::rotate( mSceneRotation );

    com.mean_position(mFlockCenter);
    // average speed of flock
    mFlockAvSpeed = com.mean_velocity(mFlockAvVel);
    // magnitude of flock polarization
    mFlockPolarization_mag = com.polarization(mFlockPolarization);

    // update camera
    mCenter = vec3(mFlockCenter[0], mFlockCenter[1], mFlockCenter[2]);
    vec3 cam_offset = mSceneRotation * vec3(0, 0, mCameraDistance);
    mEye = mCenter - cam_offset;
    mUp = mSceneRotation * vec3(0, 1, 0);

    mCam.lookAt(mEye, mCenter, mUp);
    gl::setMatrices( mCam );
	
    behavior->setJandG(mJ, mG);
    behavior->setTemp(mTemp);
    interaction->setCriticalAngle(mBalanceAngle);

    if (!mPaused)
    {
        // compute C_sp(r): speed correlation of agents as function of agent separation
        com.correlation_histo(mRadialCorMaxRange, mNumRadialBins, mV0, mC_sp_r, mC_sp_r_count);

        double sum_of_velsq = com.sense_velocities_and_velsq(v2);
        mAvNumNeighbors = com.get_av_num_neighbors();

        // update Q_int
        mQ_int = sum_of_velsq / (2 * mV0*mV0*mN*mAvNumNeighbors);
        com.update_velocities(v2);
        com.move(mDt);

        // update matrix of agent separations
        com.updateAgentSepInfo();

        double* cp = com.get_pos();
        double* cv = com.get_vel();

        // update Cinder particle objects from spp data
        for (int i = 0; i < mN; i++)
        {
            spp_particles[i].updatePos_spp(cp[3 * i], cp[3 * i + 1], cp[3 * i + 2]);
            spp_particles[i].updateVel_spp(cv[3 * i], cv[3 * i + 1], cv[3 * i + 2]);
            spp_particles[i].mVelNormal = normalize(spp_particles[i].mVel);
            spp_particles[i].mTailPos = spp_particles[i].mPos -
                                        spp_particles[i].mVelNormal * spp_particles[i].mLength;
        }
    }
}

void FlockingApp::draw()
{	
    gl::clear(Color(0, 0, 0), true);
	gl::enableDepthRead();
	gl::enableDepthWrite();
	
    if (mShowGrid)
        drawGrid(mBoxSize, 4.0f, 16.0f);
	
    gl::color(ColorA(1.0f, 1.0f, 1.0f, 1.0f));
    gl::drawStrokedCube(vec3(mBoxSize / 2, mBoxSize / 2, mBoxSize / 2), 
                        vec3(mBoxSize, mBoxSize, mBoxSize));

    //// Draw particles
    int birdFromInd = (int) (floorf(birdFrom + .0001));
    int birdToInd = (int) (floorf(birdTo + .0001));

    // no highlighting
    if (birdFromInd == -1)
    {
        gl::color(ColorA(1.0f, 1.0f, 1.0f, 1.0f));
        for (int i = 0; i < mN; i++)
            spp_particles[i].draw();
        
        gl::begin(GL_LINES);
        for (int i = 0; i < mN; i++)
            spp_particles[i].drawTail();
        gl::end();
    }
    // highlighting
    else
    {
        //// draw all but from and to birds dimmed
        // particles
        gl::color(ColorA(birdDimFactor, birdDimFactor, birdDimFactor, 1.0f));
        for (int i = 0; i < birdFromInd; i++)
            spp_particles[i].draw();
        for (int i = birdFromInd + 1; i < mN; i++)
            spp_particles[i].draw();

        // tails
        gl::begin(GL_LINES);
        for (int i = 0; i < birdFromInd; i++)
            spp_particles[i].drawTail(birdDimFactor);
        for (int i = birdFromInd + 1; i < mN; i++)
            spp_particles[i].drawTail(birdDimFactor);
        gl::end();

        // draw from bird in standard color, slightly bigger
        gl::color(ColorA(1.0f, 1.0f, 1.0f, 1.0f));
        spp_particles[birdFromInd].draw(1.6);
      
        // draw to bird in blue, slightly bigger
        gl::color(ColorA(0.0f, 0.0f, 1.0f, 1.0f));
        spp_particles[birdToInd].draw(1.6);
        
        // draw from and to bird's tails
        gl::begin(GL_LINES);
        spp_particles[birdFromInd].drawTail();
        spp_particles[birdToInd].drawTail();
        gl::end();

        // draw arrows from 'from' bird to all its neighbors
        drawBirdToNeighborsArrows(birdFromInd);
    }

    // Draw average velocity of flock
    auto avPos = vec3(mFlockCenter[0], mFlockCenter[1], mFlockCenter[2]);
    auto avVel = vec3((float) mFlockAvVel[0], (float) mFlockAvVel[1], (float) mFlockAvVel[2]);
    gl::color(ColorA(0.0f, 1.0f, 0.0f, 1.0f));
    gl::drawVector(avPos, avPos + avVel / 3.0f, 0.3f, .1f);

    // Draw pair velocity correlation as function of separation graph
    drawC_sp_Graph();
	
	// Draw Params window
	mParams->draw();
}


//////////////////////////////////////////////////////////////////////////
// Helper graphing methods
//////////////////////////////////////////////////////////////////////////

void FlockingApp::drawGrid(float boxSize, float cellSpacing, float gridRadius)
{
    vec3 begin, end;
    gl::color(ColorA(0.25f, 0.25f, 0.25f, 1.0f));

    float minX = int((mFlockCenter[0] - gridRadius) / cellSpacing) * cellSpacing;
    float maxX = (int((mFlockCenter[0] + gridRadius) / cellSpacing) + 1) * cellSpacing;
    float minY = int((mFlockCenter[1] - gridRadius) / cellSpacing) * cellSpacing;
    float maxY = (int((mFlockCenter[1] + gridRadius) / cellSpacing) + 1) * cellSpacing;
    float minZ = int((mFlockCenter[2] - gridRadius) / cellSpacing) * cellSpacing;
    float maxZ = (int((mFlockCenter[2] + gridRadius) / cellSpacing) + 1) * cellSpacing;
    
    // x lines
    for (float i = minY; i <= maxY; i = i + cellSpacing)
        for (float j = minZ; j <= maxZ; j = j + cellSpacing)
        {
            begin = vec3(minX, i, j); end = vec3(maxX, i, j);
            gl::drawLine(begin, end);
        }

    // y lines
    for (float i = minX; i <= maxX; i = i + cellSpacing)
        for (float j = minZ; j <= maxZ; j = j + cellSpacing)
        {            
            begin = vec3(i, minY, j); end = vec3(i, maxY, j);
            gl::drawLine(begin, end);
        }

    //z lines
    for (float i = minX; i <= maxX; i = i + cellSpacing)
        for (float j = minY; j <= maxY; j = j + cellSpacing)
        {            
            begin = vec3(i, j, minZ); end = vec3(i, j, maxZ);
            gl::drawLine(begin, end);
        }
}

void FlockingApp::drawC_sp_Graph()
{
    // set graph height scale to max correlation magnitude
    double scale = 0.;
    for (size_t i = 0; i < mC_sp_r.size(); i++)
    {
        double val = fabs(mC_sp_r[i]);
        if (val > scale)
            scale = val;
    }
    // dont run until initialized
    if (scale == 0.) return;

    double binW = 200.0 / mNumRadialBins;
    double barW = 0.9 * binW;
    gl::disableDepthRead();
    gl::disableDepthWrite();
    gl::setMatricesWindow(getWindowSize());
    gl::pushModelView();
      gl::translate(vec3(10.0f, getWindowHeight() - 10.0f - 600.0 / 2, 0.0f));
      gl::color(ColorA(1.0f, 1.0f, 1.0f, 1.0f));
      // draw graph axes
      gl::drawLine(vec2(0.0f, 0.0f), vec2(200.0f, 0.0f)); // x axis
      gl::drawLine(vec2(0.0f, 100.0f), vec2(0.0f, -100.0f)); // y axis
      // draw graph bars
      gl::color(ColorA(0.25f, 0.25f, 1.0f, 1.0f));
      // width of graph will be 200 pixels
      // draw graph bars
      for (int i = 1; i <= mNumRadialBins; i++)
      {
          gl::drawSolidRect(
              Rectf(i*binW - barW / 2, 0.,
                    i*binW + barW / 2, -(100 * mC_sp_r[i - 1] / scale)));
      }
    gl::popModelView();
}

void FlockingApp::drawBirdToBirdArrow(int fromIndex, int toIndex)
{
    if (toIndex >= 0 && toIndex != fromIndex)
    {
        double* cp = com.get_pos();
        double* agSepInfo = com.get_AgentSepInfo();

        auto fromPos = vec3(cp[3 * fromIndex], cp[3 * fromIndex + 1], cp[3 * fromIndex + 2]);
        auto toDisplacement = vec3(agSepInfo[fromIndex * 4 * mN + toIndex * 4],
                                   agSepInfo[fromIndex * 4 * mN + toIndex * 4 + 1],
                                   agSepInfo[fromIndex * 4 * mN + toIndex * 4 + 2]);
        gl::color(ColorA(0.0f, 1.0f, 1.0f, 1.0f));
        gl::drawVector(fromPos, fromPos + toDisplacement / 2.0f, 0.2f, .06f);
        gl::drawLine(fromPos + toDisplacement / 2.0f, fromPos + toDisplacement);
    }
}

void FlockingApp::drawBirdToNeighborsArrows(int fromIndex)
{
    if (fromIndex >= 0)
    {
        double* cp = com.get_pos();
        double* agSepInfo = com.get_AgentSepInfo();

        Agent* ags = com.get_agents();
        Agent** neis = ags[fromIndex].get_neighbor_list();
        int num_neis = ags[fromIndex].get_num_neighs();

        for (int i = 0; i < num_neis; i++)
        {
            // neighbor index is address of neighbor - address of beginning of agent array
            int toIndex = neis[i] - ags;
            auto fromPos = vec3(cp[3 * fromIndex], cp[3 * fromIndex + 1], cp[3 * fromIndex + 2]);
            auto toDisplacement = vec3(agSepInfo[fromIndex * 4 * mN + toIndex * 4],
                                       agSepInfo[fromIndex * 4 * mN + toIndex * 4 + 1],
                                       agSepInfo[fromIndex * 4 * mN + toIndex * 4 + 2]);
            gl::color(ColorA(0.0f, 1.0f, 1.0f, 1.0f));
            gl::drawVector(fromPos, fromPos + toDisplacement / 2.0f, 0.2f, .06f);
            gl::drawLine(fromPos + toDisplacement / 2.0f, fromPos + toDisplacement);
        }
    }
}

void prepareSettings(App::Settings *settings)
{
    settings->setWindowSize(1920, 1080);
    settings->setFrameRate(60.0f);
}

CINDER_APP(FlockingApp, RendererGl(), prepareSettings)
