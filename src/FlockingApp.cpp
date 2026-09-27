#include "cinder/app/App.h"
#include "cinder/app/RendererGl.h"
#include "cinder/Vector.h"
#include "cinder/Utilities.h"
#include "cinder/CinderImGui.h"
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
	void mouseDown( MouseEvent event );
	void mouseDrag( MouseEvent event );
	void mouseWheel( MouseEvent event );
	void resize();
	void setup();
	void update();
	void draw();
    void drawParamsWindow();
    void drawGrid(float boxSize, float cellSpacing, float gridRadius);
    void drawC_sp_Graph();
    void drawBirdToBirdArrow(int fromIndex, int toIndex);
    void drawBirdToNeighborsArrows(int fromIndex);

    bool     mPaused;
    bool     mShowGrid;
    // which birds to highlight (0-indexed). If birdFrom is -1, no highlighting.
    int      birdFrom;
    int      birdTo;
    float    birdDimFactor;

	// CAMERA
	CameraPersp			mCam;
	quat				mSceneRotation;
	vec3				mEye, mCenter, mUp;
	float				mCameraDistance;
	ivec2				mLastMousePos;
	
	ParticleController	mParticleController;
	gl::BatchRef		mCubeBatch;
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
    // Number of agents (birds)
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
    mQ_int = 0.0;

    birdFrom        = -1;
    birdTo          = 2;
    mPaused         = true;
    mShowGrid       = true;
    birdDimFactor   = 0.4f;

    // SETUP SIMULATION PARAMETERS
    mBoxSize = 500.;
    mN = 300;
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

	// SETUP UI (Dear ImGui; the params window itself is built each frame in drawParamsWindow())
    ImGui::Initialize();
    // ImGui works in framebuffer pixels, so scale it up on high-density (Retina) displays
    float uiScale = getWindowContentScale();
    ImGui::GetStyle().ScaleAllSizes(uiScale);
    ImGui::GetStyle().FontScaleDpi = uiScale;

    // CREATE SWARMING_SPP containers
    dist2 = spp_community_alloc_space(mN);
    v2    = spp_community_alloc_space(mN);

    g = new Cartesian(mBoxSize);
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
    mCubeBatch = gl::Batch::create(geom::Cube(), gl::getStockShader(gl::ShaderDef().color()));
    for (int i = 0; i < mN; i++)
        spp_particles.push_back(Particle(vec3(), vec3()));
}


void FlockingApp::keyDown( KeyEvent event )
{
    // keyboard shortcuts (formerly the AntTweakBar keyIncr/keyDecr bindings)
    if (event.getCode() == KeyEvent::KEY_SPACE)
    {
        mPaused = !mPaused;
        return;
    }

    switch (event.getChar())
    {
        case 'f': mJ = glm::clamp(mJ + 1.0, 1.0, 60.0); break;
        case 'd': mJ = glm::clamp(mJ - 1.0, 1.0, 60.0); break;
        case 'v': mG = glm::clamp(mG + 0.01, 0.01, 0.5); break;
        case 'c': mG = glm::clamp(mG - 0.01, 0.01, 0.5); break;
        case 'y': mTemp = glm::clamp(mTemp + 0.025, 0.025, 2.0); break;
        case 't': mTemp = glm::clamp(mTemp - 0.025, 0.025, 2.0); break;
        case 'o': mBalanceAngle = glm::clamp(mBalanceAngle + 1.0, 5.0, 180.0); break;
        case 'i': mBalanceAngle = glm::clamp(mBalanceAngle - 1.0, 5.0, 180.0); break;
        case 's': mCameraDistance = glm::clamp(mCameraDistance + 1.0f, 5.0f, 1500.0f); break;
        case 'w': mCameraDistance = glm::clamp(mCameraDistance - 1.0f, 5.0f, 1500.0f); break;
        case '.': birdFrom = glm::clamp(birdFrom + 1, -1, mN - 1); break;
        case ',': birdFrom = glm::clamp(birdFrom - 1, -1, mN - 1); break;
        case 'l': birdTo = glm::clamp(birdTo + 1, 0, mN - 1); break;
        case 'k': birdTo = glm::clamp(birdTo - 1, 0, mN - 1); break;
        case 'g': mShowGrid = !mShowGrid; break;
        case 'r': mSceneRotation = quat(); break;
    }
}

// Left-drag orbits the camera around the flock center (replaces the AntTweakBar
// "Scene Rotation" widget). Events over the ImGui window are consumed by ImGui.
void FlockingApp::mouseDown( MouseEvent event )
{
    mLastMousePos = event.getPos();
}

void FlockingApp::mouseDrag( MouseEvent event )
{
    ivec2 delta = event.getPos() - mLastMousePos;
    mLastMousePos = event.getPos();

    const float radiansPerPixel = 0.005f;
    mSceneRotation = normalize(mSceneRotation *
                               angleAxis(-delta.x * radiansPerPixel, vec3(0, 1, 0)) *
                               angleAxis( delta.y * radiansPerPixel, vec3(1, 0, 0)));
}

void FlockingApp::mouseWheel( MouseEvent event )
{
    mCameraDistance = glm::clamp(mCameraDistance * powf(0.9f, event.getWheelIncrement()),
                                 5.0f, 1500.0f);
}

void FlockingApp::resize()
{
    mCam.setAspectRatio(getWindowAspectRatio());
}

void FlockingApp::drawParamsWindow()
{
    auto sliderDouble = [](const char* label, double* value, double min, double max,
                           const char* format)
    {
        ImGui::SliderScalar(label, ImGuiDataType_Double, value, &min, &max, format,
                            ImGuiSliderFlags_AlwaysClamp);
    };

    float uiScale = getWindowContentScale();
    ImGui::SetNextWindowPos(ImVec2(10 * uiScale, 10 * uiScale), ImGuiCond_FirstUseEver);
    ImGui::Begin("Flocking", nullptr, ImGuiWindowFlags_AlwaysAutoResize);

    if (ImGui::Button("Reset view"))
        mSceneRotation = quat();
    ImGui::SameLine();
    ImGui::TextDisabled("drag: rotate, wheel: zoom");
    ImGui::Separator();

    sliderDouble("J", &mJ, 1.0, 60.0, "%.1f");
    sliderDouble("g", &mG, 0.01, 0.5, "%.2f");
    sliderDouble("Temp", &mTemp, 0.025, 2.0, "%.3f");
    sliderDouble("Balance Angle", &mBalanceAngle, 5.0, 180.0, "%.1f");
    ImGui::Separator();

    ImGui::SliderFloat("Eye Distance", &mCameraDistance, 5.0f, 1500.0f, "%.1f",
                       ImGuiSliderFlags_Logarithmic | ImGuiSliderFlags_AlwaysClamp);
    ImGui::SliderInt("Bird From", &birdFrom, -1, mN - 1, "%d", ImGuiSliderFlags_AlwaysClamp);
    ImGui::SliderInt("Bird To", &birdTo, 0, mN - 1, "%d", ImGuiSliderFlags_AlwaysClamp);
    ImGui::Checkbox("Paused", &mPaused);
    ImGui::Checkbox("Draw grid", &mShowGrid);
    ImGui::Separator();

    ImGui::Text("Mean speed:     %.4f", mFlockAvSpeed);
    ImGui::Text("Av Num Neighs:  %.4f", mAvNumNeighbors);
    ImGui::Text("Polarization:   %.4f", mFlockPolarization_mag);
    ImGui::Text("Q_int:          %.4f", mQ_int);

    if (ImGui::CollapsingHeader("Keyboard shortcuts"))
    {
        ImGui::TextUnformatted(
            "space  pause / run\n"
            "f / d  J +/-\n"
            "v / c  g +/-\n"
            "y / t  Temp +/-\n"
            "o / i  Balance Angle +/-\n"
            "s / w  Eye Distance +/-\n"
            ". / ,  Bird From +/-\n"
            "l / k  Bird To +/-\n"
            "g      toggle grid\n"
            "r      reset view");
    }

    ImGui::End();
}


void FlockingApp::update()
{
    drawParamsWindow();

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
    int birdFromInd = birdFrom;
    int birdToInd = birdTo;

    // no highlighting
    if (birdFromInd == -1)
    {
        gl::color(ColorA(1.0f, 1.0f, 1.0f, 1.0f));
        for (int i = 0; i < mN; i++)
            spp_particles[i].draw(mCubeBatch);
        
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
            spp_particles[i].draw(mCubeBatch);
        for (int i = birdFromInd + 1; i < mN; i++)
            spp_particles[i].draw(mCubeBatch);

        // tails
        gl::begin(GL_LINES);
        for (int i = 0; i < birdFromInd; i++)
            spp_particles[i].drawTail(birdDimFactor);
        for (int i = birdFromInd + 1; i < mN; i++)
            spp_particles[i].drawTail(birdDimFactor);
        gl::end();

        // draw from bird in standard color, slightly bigger
        gl::color(ColorA(1.0f, 1.0f, 1.0f, 1.0f));
        spp_particles[birdFromInd].draw(mCubeBatch, 1.6);
      
        // draw to bird in blue, slightly bigger
        gl::color(ColorA(0.0f, 0.0f, 1.0f, 1.0f));
        spp_particles[birdToInd].draw(mCubeBatch, 1.6);
        
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

	// (the ImGui params window is rendered automatically after draw())
}


//////////////////////////////////////////////////////////////////////////
// Helper graphing methods
//////////////////////////////////////////////////////////////////////////

void FlockingApp::drawGrid(float boxSize, float cellSpacing, float gridRadius)
{
    gl::color(ColorA(0.25f, 0.25f, 0.25f, 1.0f));

    float minX = int((mFlockCenter[0] - gridRadius) / cellSpacing) * cellSpacing;
    float maxX = (int((mFlockCenter[0] + gridRadius) / cellSpacing) + 1) * cellSpacing;
    float minY = int((mFlockCenter[1] - gridRadius) / cellSpacing) * cellSpacing;
    float maxY = (int((mFlockCenter[1] + gridRadius) / cellSpacing) + 1) * cellSpacing;
    float minZ = int((mFlockCenter[2] - gridRadius) / cellSpacing) * cellSpacing;
    float maxZ = (int((mFlockCenter[2] + gridRadius) / cellSpacing) + 1) * cellSpacing;
    
    // all lines go in a single batch; one gl::drawLine() per segment is very slow on macOS
    gl::begin(GL_LINES);

    // x lines
    for (float i = minY; i <= maxY; i = i + cellSpacing)
        for (float j = minZ; j <= maxZ; j = j + cellSpacing)
        {
            gl::vertex(vec3(minX, i, j)); gl::vertex(vec3(maxX, i, j));
        }

    // y lines
    for (float i = minX; i <= maxX; i = i + cellSpacing)
        for (float j = minZ; j <= maxZ; j = j + cellSpacing)
        {
            gl::vertex(vec3(i, minY, j)); gl::vertex(vec3(i, maxY, j));
        }

    //z lines
    for (float i = minX; i <= maxX; i = i + cellSpacing)
        for (float j = minY; j <= maxY; j = j + cellSpacing)
        {
            gl::vertex(vec3(i, j, minZ)); gl::vertex(vec3(i, j, maxZ));
        }

    gl::end();
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
    // On current macOS the GL view gets a Retina-resolution drawable regardless, so opt in
    // explicitly; otherwise Cinder reports a content scale of 1 and ImGui renders at the
    // wrong scale.
    settings->setHighDensityDisplayEnabled(true);
}

CINDER_APP(FlockingApp, RendererGl(), prepareSettings)
