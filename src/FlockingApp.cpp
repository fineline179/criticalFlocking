#include "cinder/app/AppBasic.h"
#include "cinder/Vector.h"
#include "cinder/Utilities.h"
#include "cinder/params/Params.h"
#include "cinder/Camera.h"
#include "ParticleController.h"
#include "swarming_spp/community.h"
#include "cinder/Rand.h"
#include "swarming_spp/grid.h"
#include <math.h>

using namespace ci;
using namespace ci::app;

class FlockingApp : public AppBasic {
 public:
	void prepareSettings( Settings *settings );
	void keyDown( KeyEvent event );
	void setup();
	void update();
	void draw();
    void drawGrid(float boxSize, float cellSpacing, float gridRadius);
    void drawC_sp_Graph();
    void drawBirdToBirdArrow(int fromIndex, int toIndex);

    bool     mPaused;
    bool     mShowGrid;
    float    birdFrom;
    float    birdTo;
    float    birdDimFactor;

	// PARAMS
	params::InterfaceGlRef	mParams;
	
	// CAMERA
	CameraPersp			mCam;
	Quatf				mSceneRotation;
	Vec3f				mEye, mCenter, mUp;
	float				mCameraDistance;
	
	ParticleController	mParticleController;
	float				mZoneRadius;
	bool				mCentralGravity;
	bool				mFlatten;

    //// SWARMING_SPP
    Cartesian*           g;
    //Topologic*           interaction;
    TopoBalanced*          interaction;
    Bialek_consensus*    behavior;
    Grid*                grid;
    Community            com;

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
    double*              dist2;
    double*              v2;

    std::vector<Particle>  spp_particles;

    int mRadialCorMaxRange;
    int mNumRadialBins;
    // Correlation of speed of bird pairs as func of separation
    //double* mC_sp_r;
    std::vector<double> mC_sp_r;
    // Counts of bird pairs as a function of separation (used in calc of previous quantity
    //int*    mC_sp_r_count;
    std::vector<int>    mC_sp_r_count;
    double* mFlockCenter;
    double* mFlockAvVel;
    double  mFlockAvSpeed;
    double  mAvNumNeighbors;
    double* mFlockPolarization;
    double  mFlockPolarization_mag;

    // similarity of velocity of nearest neighbors in flock. Eq. (1) in 1307.5563
    double mQ_int;
};


void FlockingApp::prepareSettings( Settings *settings )
{
    //settings->setWindowSize(1280, 720);
    settings->setWindowSize(1920, 1080);
	settings->setFrameRate( 60.0f );
}


void FlockingApp::setup()
{	
    // SETUP STATISTICAL VARIABLES
    mNumRadialBins = 20;
    mRadialCorMaxRange = 25.0;
    //mC_sp_r = new double[20];
    //mC_sp_r_count = new int[20];
    mC_sp_r = std::vector<double>(20);
    mC_sp_r_count = std::vector<int>(20);
    mFlockCenter = new double[3];
    mFlockAvVel = new double[3];
    mAvNumNeighbors = 0.0;
    mFlockPolarization = new double[3];

    birdFrom        = 0;
    birdTo          = 0;
    mPaused         = true;
    mShowGrid       = true;
    birdDimFactor   = 0.25f;

    // SETUP SIMULATION PARAMETERS
    mBoxSize = 500.;
    mN = 442;
    mJ = 10.0; mG = 0.08;
    mDt = .03;
    mV0 = 10.13;
    //mTemp = 0.93;
    mTemp = 0.5;
    m_nc = 8.;
    mBalanceAngle = 32.0;
    ra = 0.8; rb = 0.2; re = 0.5; 
    // make this large for unlimited small attraction range on non-balanced topological interactions
    r0 = 1000.;

	// SETUP CAMERA
    mCameraDistance = 20.0f;
    mEye = Vec3f(0.0, 0.0, mCameraDistance);
    mCenter = Vec3f(0.0, 0.0, 0.0);
    mUp = Vec3f::yAxis();
	mCam.setPerspective( 75.0f, getWindowAspectRatio(), 5.0f, 2000.0f );

	// SETUP UI PARAMS
    mParams = params::InterfaceGl::create("Flocking", Vec2i(200, 400));
    mParams->addParam("Scene Rotation", &mSceneRotation, "opened=1");
    mParams->addSeparator();
    mParams->addParam("J", &mJ, "min=1.0 max=60.0 step=1.0 keyIncr=f keyDecr=d");
    mParams->addParam("g",  &mG,  "min=0.1 max=20.0 step=0.1 keyIncr=v keyDecr=c");
    mParams->addParam("Temp", &mTemp, "min=0.025 max=2.0 step=0.025 keyIncr=y keyDecr=t");
    mParams->addParam("nc", &m_nc, "min=2.0 max=25.0 step=1.0 keyIncr=r keyDecr=e");
    mParams->addParam("Balance Angle", &mBalanceAngle, 
                      "min=5.0 max=180.0 step=1.0 keyIncr=o keyDecr=i");
    mParams->addSeparator();
    mParams->addParam("Eye Distance", &mCameraDistance,
                      "min=5.0 max=1500.0 step=1.0 keyIncr=s keyDecr=w" );
    // which bird to highlight (1-indexed). If 0, no highlighting.
    // TODO: unhardcode max and replace with NUM_INITIAL_PARTICLES
    mParams->addParam("Bird From", &birdFrom,
                      "min=0.0 max=10.0 step=1.0 keyIncr=. keyDecr=,");
    mParams->addParam("Bird To", &birdTo,
                      "min=0.0 max=10.0 step=1.0 keyIncr=l keyDecr=k");
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
    
    // create particles
    for (int i = 0; i < mN; i++)
        spp_particles.push_back(Particle(Vec3f(), Vec3f()));
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

    // UPDATE CAMERA
    mCenter = Vec3f(mFlockCenter[0], mFlockCenter[1], mFlockCenter[2]);
    Vec3f cam_offset = mSceneRotation * Vec3f(0, 0, mCameraDistance);
    mEye = mCenter - cam_offset;
    mUp = mSceneRotation * Vec3f(0, 1, 0);

    mCam.lookAt(mEye, mCenter, mUp);
    gl::setMatrices( mCam );
	
    behavior->setJandG(mJ, mG);
    behavior->setTemp(mTemp);
    interaction->setOutdegree((int) m_nc);
    interaction->setCriticalAngle(mBalanceAngle);

    if (!mPaused)
    {
        // compute C_sp(r): speed correlation of agents as function of agent separation
        com.correlation_histo(mRadialCorMaxRange, mNumRadialBins, mV0, mC_sp_r, mC_sp_r_count);

        // update spp particle velocities
        double sum_of_velsq = com.sense_velocities_and_velsq(v2);
        mAvNumNeighbors = com.get_av_num_neighbors();

        // update Q_int
        mQ_int = sum_of_velsq / (2 * mV0*mV0*mN*m_nc);
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
            spp_particles[i].mVelNormal = spp_particles[i].mVel.normalized();
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
        drawGrid(mBoxSize, 2.0f, 5.0f);
	
    gl::color(ColorA(1.0f, 1.0f, 1.0f, 1.0f));
    gl::drawStrokedCube(Vec3f(mBoxSize / 2, mBoxSize / 2, mBoxSize / 2), 
                        Vec3f(mBoxSize, mBoxSize, mBoxSize));

    // DRAW SWARMING_SPP PARTICLES
    // no highlighting
    if ((int) birdFrom == 0)
    {
        gl::color(ColorA(1.0f, 1.0f, 1.0f, 1.0f));
        for (int i = 0; i < mN; i++)
            spp_particles[i].draw();
        
        glBegin(GL_LINES);
        for (int i = 0; i < mN; i++)
            spp_particles[i].drawTail();
        glEnd();
    }

    // highlighting
    else
    {
        //// draw all but from and to birds dimmed
        // particles
        gl::color(ColorA(birdDimFactor, birdDimFactor, birdDimFactor, 1.0f));
        for (int i = 0; i < (int) birdFrom - 1; i++)
            spp_particles[i].draw();
        for (int i = (int) birdFrom; i < mN; i++)
            spp_particles[i].draw();

        // tails
        glBegin(GL_LINES);
        for (int i = 0; i < (int) birdFrom - 1; i++)
            spp_particles[i].drawTail(birdDimFactor);
        for (int i = (int) birdFrom; i < mN; i++)
            spp_particles[i].drawTail(birdDimFactor);
        glEnd();

        int birdFromIndex = (int) birdFrom - 1;
        int birdToIndex   = (int) birdTo - 1;

        // draw from bird in standard color, slightly bigger
        gl::color(ColorA(1.0f, 1.0f, 1.0f, 1.0f));
        spp_particles[birdFromIndex].draw(1.6);
        // draw to bird in blue, slightly bigger
        gl::color(ColorA(0.0f, 0.0f, 1.0f, 1.0f));
        spp_particles[birdToIndex].draw(1.6);
        // draw from and to bird's tails
        glBegin(GL_LINES);
        spp_particles[birdFromIndex].drawTail();
        spp_particles[birdToIndex].drawTail();
        glEnd();

        // draw arrow from 'from' bird to 'to' bird
        drawBirdToBirdArrow(birdFromIndex, birdToIndex);
    }

    // Draw average velocity of flock
    auto avPos = Vec3f(mFlockCenter[0], mFlockCenter[1], mFlockCenter[2]);
    auto avVel = Vec3f((float) mFlockAvVel[0], (float) mFlockAvVel[1], (float) mFlockAvVel[2]);
    gl::color(ColorA(0.0f, 1.0f, 0.0f, 1.0f));
    gl::drawVector(avPos, avPos + avVel / 3, 0.3, .1);

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
    Vec3f begin, end;
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
            begin.set(minX, i, j); end.set(maxX, i, j);
            gl::drawLine(begin, end);
        }

    // y lines
    for (float i = minX; i <= maxX; i = i + cellSpacing)
        for (float j = minZ; j <= maxZ; j = j + cellSpacing)
        {
            begin.set(i, minY, j); end.set(i, maxY, j);
            gl::drawLine(begin, end);
        }

    //z lines
    for (float i = minX; i <= maxX; i = i + cellSpacing)
        for (float j = minY; j <= maxY; j = j + cellSpacing)
        {
            begin.set(i, j, minZ); end.set(i, j, maxZ);
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
      gl::translate(Vec3f(10.0f, getWindowHeight() - 10.0f - 600.0 / 2, 0.0f));
      gl::color(ColorA(1.0f, 1.0f, 1.0f, 1.0f));
      // draw graph axes
      gl::drawLine(Vec2f(0.0f, 0.0f), Vec2f(200.0f, 0.0f)); // x axis
      gl::drawLine(Vec2f(0.0f, 100.0f), Vec2f(0.0f, -100.0f)); // y axis
      // draw graph bars
      gl::color(ColorA(0.25f, 0.25f, 1.0f, 1.0f));
      // width of graph will be 200 pixels
      // draw graph bars
      for (int i = 1; i <= mNumRadialBins; i++)
      {
          gl::drawSolidRect(
              Rectf(i*binW - barW / 2, 0.,
                    i*binW + barW / 2, -(100 * mC_sp_r[i - 1] / scale))
          );
      }
    gl::popModelView();
}

void FlockingApp::drawBirdToBirdArrow(int fromIndex, int toIndex)
{
    if (toIndex >= 0 && toIndex != fromIndex)
    {
        double* cp = com.get_pos();
        double* agSepInfo = com.get_AgentSepInfo();

        auto fromPos = Vec3f(cp[3 * fromIndex], cp[3 * fromIndex + 1], cp[3 * fromIndex + 2]);
        auto toDisplacement = Vec3f(agSepInfo[fromIndex * 4 * mN + toIndex * 4],
                                    agSepInfo[fromIndex * 4 * mN + toIndex * 4 + 1],
                                    agSepInfo[fromIndex * 4 * mN + toIndex * 4 + 2]);
        gl::color(ColorA(0.0f, 1.0f, 1.0f, 1.0f));
        gl::drawVector(fromPos, fromPos + toDisplacement / 2, 0.2, .06);
        gl::drawLine(fromPos + toDisplacement / 2, fromPos + toDisplacement);
    }
}


CINDER_APP_BASIC(FlockingApp, RendererGl)
