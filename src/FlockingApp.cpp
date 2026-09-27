#include "cinder/app/App.h"
#include "cinder/app/RendererGl.h"
#include "cinder/CinderImGui.h"
#include "cinder/Camera.h"
#include "cinder/gl/gl.h"
#include "sim/Flock.h"
#include "sim/Measure.h"
#include "sim/Predator.h"
#include <deque>
#include <memory>

using namespace ci;
using namespace ci::app;

// Real-time flock simulation with a predator. The physics lives in src/sim (see
// docs/SIMULATION_DESIGN.md); this file only drives it, draws it, and shows the controls and the
// validation panel.
class FlockingApp : public App {
 public:
	void setup();
	void update();
	void draw();
	void keyDown( KeyEvent event );
	void mouseDown( MouseEvent event );
	void mouseDrag( MouseEvent event );
	void mouseWheel( MouseEvent event );
	void resize();

 private:
    enum class ColorMode { TurnedSinceAttack, TurnRate, SpeedDeviation, Detection, Plain };

    void resetFlock();
    void applyPresets(bool resetState);
    void refreshMeasurements();
    void drawParamsWindow();
    void drawPredatorWindow();
    void drawValidationWindow();
    void drawGrid(const vec3 &center, float cellSpacing, float gridRadius);
    void drawCorrelationGraphs();
    void drawArrowHead(const vec3 &end, const vec3 &dir, float headLength, float headRadius);
    ColorA birdColor(int i) const;
    vec3   toVec(const sim::Vec3 &v) const { return vec3((float) v.x, (float) v.y, (float) v.z); }

    // simulation
    std::unique_ptr<sim::Flock> mFlock;
    sim::Predator               mPredator;
    sim::TurningPreset          mTurningPreset = sim::TurningPreset::Inertial;
    sim::SpeedPreset            mSpeedPreset   = sim::SpeedPreset::Marginal;
    int                         mN             = 300;
    uint64_t                    mSeed          = 1;
    bool                        mPaused        = false;
    float                       mTimeScale     = 0.5f;   // simulated seconds per real second

    // measurements
    sim::Stats                      mStats;
    double                          mLastMeasure = -1;
    std::unique_ptr<sim::FluctuationProbe> mProbe;
    sim::FluctuationProbe::Result   mProbeResult;
    double                          mProbeSince = 0, mProbeLastResult = 0;
    std::deque<sim::Vec3>           mPredatorTrail;
    // headings over the last ~0.25 s of simulated time, for a turn rate that ignores jitter
    std::deque<std::pair<double, std::vector<sim::Vec3>>> mHeadingHistory;
    std::vector<float>              mTurnRate;
    // time-averaged correlation functions for the graphs
    std::vector<double>             mCdirAvg, mCspAvg;

    // display
    ColorMode mColorMode   = ColorMode::TurnedSinceAttack;
    // flock heading when the latest attack started (reference for the "turned" color mode)
    sim::Vec3 mAttackHeading;
    int       mAttacksSeen = 0;
    bool      mDetectionSeen = false;
    bool      mShowGrid    = true;
    int       mHighlight   = -1;   // show this bird's neighbor links (-1: none)
    float     mBirdRadius  = 0.15f;
    float     mTailLength  = 0.6f;

	// camera
	CameraPersp mCam;
	quat        mSceneRotation;
	float       mCameraDistance = 35.0f;
	ivec2       mLastMousePos;
	vec3        mCenter;

	gl::BatchRef mSphereBatch;
	gl::BatchRef mConeBatch;
};


void FlockingApp::setup()
{
    mCam.setPerspective(60.0f, getWindowAspectRatio(), 0.5f, 3000.0f);

    ImGui::Initialize();
    // ImGui works in framebuffer pixels, so scale it up on high-density (Retina) displays
    float uiScale = getWindowContentScale();
    ImGui::GetStyle().ScaleAllSizes(uiScale);
    ImGui::GetStyle().FontScaleDpi = uiScale;

    mSphereBatch = gl::Batch::create(geom::Sphere().subdivisions(8), gl::getStockShader(gl::ShaderDef().color()));
    // unit cone along +y with its base at the origin, for arrowheads and the predator
    mConeBatch = gl::Batch::create(geom::Cone().base(1.0f).apex(0.0f).height(1.0f)
                                       .origin(vec3(0)).direction(vec3(0, 1, 0)),
                                   gl::getStockShader(gl::ShaderDef().color()));
    resetFlock();
}

void FlockingApp::resetFlock()
{
    mFlock = std::make_unique<sim::Flock>(sim::makePreset(mTurningPreset, mSpeedPreset, mN), mSeed++);
    mProbe = std::make_unique<sim::FluctuationProbe>(0.02, std::vector<double>{ 1.5 }, 2.0, 20.0);
    mProbeResult = {};
    mPredator    = sim::Predator(mSeed + 1000);
    mPredatorTrail.clear();
    mHeadingHistory.clear();
    mTurnRate.assign(mN, 0.0f);
    mCdirAvg.clear();
    mCspAvg.clear();
    mLastMeasure = -1;
    refreshMeasurements();
}

// Switch presets on the running flock (keeps positions and headings)
void FlockingApp::applyPresets(bool resetState)
{
    if (resetState) { resetFlock(); return; }
    mFlock->setParams(sim::makePreset(mTurningPreset, mSpeedPreset, mN));
    mProbe->clear();
    mProbeResult = {};
}

void FlockingApp::refreshMeasurements()
{
    mStats       = sim::measure(*mFlock);
    mLastMeasure = mFlock->time();

    // exponential moving average of the correlation functions (they are noisy snapshot to snapshot)
    auto blend = [](std::vector<double> &avg, const std::vector<double> &c) {
        if (avg.size() != c.size()) avg.resize(c.size(), 0.0);
        for (size_t b = 0; b < c.size(); b++) avg[b] = avg[b] == 0.0 ? c[b] : 0.85 * avg[b] + 0.15 * c[b];
    };
    blend(mCdirAvg, mStats.Cdir);
    blend(mCspAvg, mStats.Csp);
}

void FlockingApp::keyDown( KeyEvent event )
{
    if (event.getCode() == KeyEvent::KEY_SPACE) {
        mPaused = !mPaused;
        return;
    }
    switch (event.getChar()) {
        case 'a': mPredator.triggerAttack(*mFlock); break;
        case 'n': resetFlock(); break;
        case 'r': mSceneRotation = quat(); break;
        case 'g': mShowGrid = !mShowGrid; break;
        case 'c': mColorMode = ColorMode(((int) mColorMode + 1) % 5); break;
        case 'f': mPredator.params.falseAlarmRate = mPredator.params.falseAlarmRate > 0 ? 0.0 : 0.3; break;
        case '[': mTimeScale = glm::clamp(mTimeScale * 0.8f, 0.05f, 1.0f); break;
        case ']': mTimeScale = glm::clamp(mTimeScale * 1.25f, 0.05f, 1.0f); break;
        case '1': mTurningPreset = sim::TurningPreset::Overdamped; applyPresets(false); break;
        case '2': mTurningPreset = sim::TurningPreset::Inertial; applyPresets(false); break;
        case '3': mTurningPreset = sim::TurningPreset::NonlinearInertial; applyPresets(false); break;
        case '4': mSpeedPreset = sim::SpeedPreset::Stiff; applyPresets(false); break;
        case '5': mSpeedPreset = sim::SpeedPreset::NearCritical; applyPresets(false); break;
        case '6': mSpeedPreset = sim::SpeedPreset::Marginal; applyPresets(false); break;
    }
}

// Left-drag orbits the camera around the flock center; the wheel zooms. Events over the ImGui
// windows are consumed by ImGui.
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

void FlockingApp::update()
{
    if (!mPaused) {
        // advance in chunks of at most 1/60 s of simulated time, so the predator and the
        // measurements update smoothly even in real time
        double remaining = mTimeScale / 60.0;
        while (remaining > 1e-9) {
            double dt = std::min(remaining, 1.0 / 60.0);
            mPredator.update(*mFlock, dt);
            mFlock->advance(dt);
            remaining -= dt;

            mProbeSince += dt;
            if (mProbeSince >= 0.02) {
                mProbe->sample(*mFlock);
                mProbeSince = 0;
            }
        }
        // turn rate = heading change over the last ~0.25 s (fast wobbles cancel out)
        mHeadingHistory.push_back({ mFlock->time(), mFlock->u });
        while (mHeadingHistory.size() > 2 && mFlock->time() - mHeadingHistory[1].first > 0.25)
            mHeadingHistory.pop_front();
        const auto &old = mHeadingHistory.front();
        double      span = mFlock->time() - old.first;
        mTurnRate.assign(mFlock->size(), 0.0f);
        if (span > 0 && (int) old.second.size() == mFlock->size())
            for (int i = 0; i < mFlock->size(); i++)
                mTurnRate[i] = (float) (std::acos(glm::clamp(sim::dot(old.second[i], mFlock->u[i]), -1.0, 1.0)) / span);

        // reference heading for the "turned" colors: reset when an attack starts and again when the
        // first bird detects the predator, so the colors show the response, not earlier drift
        const auto &rep = mPredator.report();
        if (rep.attacks != mAttacksSeen || (!mDetectionSeen && rep.detectedEver > 0)) {
            mAttacksSeen   = rep.attacks;
            mDetectionSeen = rep.detectedEver > 0;
            mAttackHeading = sim::measure(*mFlock, false).meanHeading;
        }

        if (mPredator.active()) {
            mPredatorTrail.push_back(mPredator.position());
            if (mPredatorTrail.size() > 120) mPredatorTrail.pop_front();
        }
        else
            mPredatorTrail.clear();

        if (mFlock->time() - mLastMeasure > 0.25) refreshMeasurements();
        if (mFlock->time() - mProbeLastResult > 1.0) {
            mProbeResult     = mProbe->result();
            mProbeLastResult = mFlock->time();
        }
    }

    drawParamsWindow();
    drawPredatorWindow();
    drawValidationWindow();

    // camera follows the flock center
    mCenter          = toVec(mStats.center);
    vec3 camOffset   = mSceneRotation * vec3(0, 0, mCameraDistance);
    mCam.lookAt(mCenter - camOffset, mCenter, mSceneRotation * vec3(0, 1, 0));
}

ColorA FlockingApp::birdColor(int i) const
{
    switch (mColorMode) {
        case ColorMode::TurnedSinceAttack: {
            // angle from the flock's heading when the latest attack started: 0 white -> 60 deg red
            if (mAttacksSeen == 0) return ColorA(1, 1, 1, 1);
            float angle = (float) std::acos(glm::clamp(sim::dot(mFlock->u[i], mAttackHeading), -1.0, 1.0));
            float t     = glm::clamp(angle / 1.047f, 0.0f, 1.0f);
            return ColorA(1.0f, 1.0f - t, 1.0f - t, 1.0f);
        }
        case ColorMode::TurnRate: {
            // 0 rad/s white -> 2 rad/s red
            float rate = i < (int) mTurnRate.size() ? mTurnRate[i] : 0.0f;
            float t    = glm::clamp((rate - 0.8f) / 2.0f, 0.0f, 1.0f);   // dead zone for wobbles
            return ColorA(1.0f, 1.0f - t, 1.0f - t, 1.0f);
        }
        case ColorMode::SpeedDeviation: {
            // slower than the flock -> blue, faster -> red, +-20%
            float d = (float) ((mFlock->speed(i) - mStats.meanSpeed) / mStats.meanSpeed / 0.2);
            d       = glm::clamp(d, -1.0f, 1.0f);
            return d > 0 ? ColorA(1.0f, 1.0f - d, 1.0f - d, 1.0f) : ColorA(1.0f + d, 1.0f + d, 1.0f, 1.0f);
        }
        case ColorMode::Detection:
            if (mPredator.detecting(i)) return ColorA(1.0f, 0.15f, 0.15f, 1.0f);
            if (mPredator.falseAlarm(i)) return ColorA(1.0f, 0.6f, 0.0f, 1.0f);
            return ColorA(0.75f, 0.75f, 0.75f, 1.0f);
        case ColorMode::Plain: break;
    }
    return ColorA(1, 1, 1, 1);
}

void FlockingApp::draw()
{
    gl::clear(Color(0, 0, 0), true);
    gl::setMatrices(mCam);
    gl::enableDepthRead();
    gl::enableDepthWrite();

    if (mShowGrid) drawGrid(mCenter, 5.0f, 30.0f);

    const sim::Flock &f = *mFlock;
    // birds: spheres plus a short tail behind each one
    for (int i = 0; i < f.size(); i++) {
        gl::color(birdColor(i));
        gl::ScopedModelMatrix scpModel;
        gl::translate(toVec(f.x[i]));
        gl::scale(vec3(mBirdRadius));
        mSphereBatch->draw();
    }
    gl::begin(GL_LINES);
    for (int i = 0; i < f.size(); i++) {
        vec3 p = toVec(f.x[i]);
        gl::color(ColorA(0.8f, 0.8f, 0.8f, 1.0f));
        gl::vertex(p);
        gl::color(ColorA(0.6f, 0.0f, 0.0f, 1.0f));
        gl::vertex(p - toVec(f.u[i]) * mTailLength);
    }
    gl::end();

    // neighbor links of a highlighted bird (who it listens to)
    if (mHighlight >= 0 && mHighlight < f.size() && !f.directed.empty()) {
        vec3 p = toVec(f.x[mHighlight]);
        gl::color(ColorA(0.0f, 1.0f, 1.0f, 1.0f));
        gl::begin(GL_LINES);
        for (int j : f.directed[mHighlight]) {
            gl::vertex(p);
            gl::vertex(toVec(f.x[j]));
        }
        gl::end();
        for (int j : f.directed[mHighlight]) {
            vec3 d = toVec(f.x[j]) - p;
            drawArrowHead(p + d * 0.5f, d, 0.3f, 0.1f);
        }
    }

    // predator: a red cone pointing along its velocity, with a trail
    if (mPredator.active()) {
        gl::color(ColorA(1.0f, 0.1f, 0.1f, 1.0f));
        vec3 v = toVec(mPredator.velocity());
        drawArrowHead(toVec(mPredator.position()) + normalize(v) * 0.8f, v, 1.6f, 0.4f);
        gl::begin(GL_LINE_STRIP);
        for (const auto &q : mPredatorTrail) gl::vertex(toVec(q));
        gl::end();
    }

    // flock's mean velocity
    gl::color(ColorA(0.0f, 1.0f, 0.0f, 1.0f));
    vec3 mv = toVec(mStats.meanHeading) * (float) (mStats.meanSpeed * 0.3);
    gl::drawLine(mCenter, mCenter + mv);
    drawArrowHead(mCenter + mv, mv, 0.8f, 0.25f);

    drawCorrelationGraphs();
}

// Same arrowhead that gl::drawVector() draws for an arrow ending at 'end' along 'dir', but from
// a cached batch, which avoids drawVector's per-call vertex upload (slow on macOS).
void FlockingApp::drawArrowHead(const vec3 &end, const vec3 &dir, float headLength, float headRadius)
{
    gl::ScopedModelMatrix scpModel;
    gl::translate(end - normalize(dir) * headLength);
    gl::rotate(glm::rotation(vec3(0, 1, 0), normalize(dir)));
    gl::scale(vec3(headRadius, headLength, headRadius));
    mConeBatch->draw();
}

void FlockingApp::drawGrid(const vec3 &center, float cellSpacing, float gridRadius)
{
    gl::color(ColorA(0.13f, 0.13f, 0.13f, 1.0f));
    auto lo = [&](float c) { return std::floor((c - gridRadius) / cellSpacing) * cellSpacing; };
    auto hi = [&](float c) { return std::ceil((c + gridRadius) / cellSpacing) * cellSpacing; };
    float minX = lo(center.x), maxX = hi(center.x);
    float minY = lo(center.y), maxY = hi(center.y);
    float minZ = lo(center.z), maxZ = hi(center.z);

    // all lines go in a single batch; one gl::drawLine() per segment is very slow on macOS
    gl::begin(GL_LINES);
    for (float i = minY; i <= maxY; i += cellSpacing)
        for (float j = minZ; j <= maxZ; j += cellSpacing) {
            gl::vertex(vec3(minX, i, j)); gl::vertex(vec3(maxX, i, j));
        }
    for (float i = minX; i <= maxX; i += cellSpacing)
        for (float j = minZ; j <= maxZ; j += cellSpacing) {
            gl::vertex(vec3(i, minY, j)); gl::vertex(vec3(i, maxY, j));
        }
    for (float i = minX; i <= maxX; i += cellSpacing)
        for (float j = minY; j <= maxY; j += cellSpacing) {
            gl::vertex(vec3(i, j, minZ)); gl::vertex(vec3(i, j, maxZ));
        }
    gl::end();
}

// Correlation functions of heading (top) and speed (bottom) fluctuations vs distance, normalized
// to their first bin, with the zero crossing (correlation length) marked.
void FlockingApp::drawCorrelationGraphs()
{
    gl::disableDepthRead();
    gl::disableDepthWrite();
    gl::setMatricesWindow(getWindowSize());

    auto crossing = [&](const std::vector<double> &C) {
        for (size_t b = 1; b < C.size(); b++)
            if (C[b - 1] > 0 && C[b] <= 0)
                return (b - 0.5 + C[b - 1] / (C[b - 1] - C[b])) * mStats.binWidth;
        return 0.0;
    };
    auto graph = [&](const std::vector<double> &C, double xi, float top, ColorA color) {
        if (C.empty() || C[0] == 0) return;
        const float width = 260.0f, halfHeight = 50.0f;
        const float binW  = width / (float) C.size();
        gl::pushModelView();
        gl::translate(vec3(15.0f, top + halfHeight, 0.0f));
        gl::color(ColorA(0.7f, 0.7f, 0.7f, 1.0f));
        gl::begin(GL_LINES);
        gl::vertex(vec2(0, 0)); gl::vertex(vec2(width, 0));
        gl::vertex(vec2(0, -halfHeight)); gl::vertex(vec2(0, halfHeight));
        gl::end();
        gl::color(color);
        gl::begin(GL_TRIANGLES);
        for (size_t b = 0; b < C.size(); b++) {
            float h = (float) glm::clamp(C[b] / C[0], -1.0, 1.0) * -halfHeight;
            Rectf bar(b * binW + 0.1f * binW, 0.0f, (b + 0.9f) * binW, h);
            gl::vertex(bar.getUpperLeft()); gl::vertex(bar.getUpperRight()); gl::vertex(bar.getLowerRight());
            gl::vertex(bar.getUpperLeft()); gl::vertex(bar.getLowerRight()); gl::vertex(bar.getLowerLeft());
        }
        gl::end();
        if (xi > 0) {
            float x = (float) (xi / mStats.binWidth) * binW;
            gl::color(ColorA(1.0f, 1.0f, 0.0f, 1.0f));
            gl::drawLine(vec2(x, -halfHeight), vec2(x, halfHeight));
        }
        gl::popModelView();
    };
    float h = (float) getWindowHeight();
    graph(mCdirAvg, crossing(mCdirAvg), h - 250.0f, ColorA(0.3f, 0.8f, 0.3f, 1.0f));
    graph(mCspAvg, crossing(mCspAvg), h - 130.0f, ColorA(0.3f, 0.3f, 1.0f, 1.0f));
}

void FlockingApp::drawParamsWindow()
{
    float uiScale = getWindowContentScale();
    ImGui::SetNextWindowPos(ImVec2(10 * uiScale, 10 * uiScale), ImGuiCond_FirstUseEver);
    ImGui::Begin("Flock", nullptr, ImGuiWindowFlags_AlwaysAutoResize);
    ImGui::PushTextWrapPos(ImGui::GetFontSize() * 26.0f);

    // presets
    const sim::TurningPreset turnings[] = { sim::TurningPreset::Overdamped, sim::TurningPreset::Inertial,
                                            sim::TurningPreset::NonlinearInertial };
    const char *turnHelp[] = {
        "The paper's Appendix G dynamics: no turning inertia. Same snapshot statistics as B. "
        "Published model results and turning data say this can't carry a turn across a flock. [model]",
        "Inertial Spin Model (Cavagna et al. 2015): turns travel as waves, matching starling data "
        "(Attanasi et al. 2014). Known mismatch: its small spontaneous wobbles oscillate, while real "
        "ones don't (Cavagna et al., PRL 2026). [model]",
        "ISM + quartic alignment (Cavagna et al., PRL 2026). Newest and least tested. In our tests at "
        "realistic polarization, turns sometimes don't spread at all. [new, experimental]" };
    ImGui::Text("Turning dynamics");
    for (int k = 0; k < 3; k++) {
        if (ImGui::RadioButton(sim::toString(turnings[k]), mTurningPreset == turnings[k])) {
            mTurningPreset = turnings[k];
            applyPresets(false);
        }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", turnHelp[k]);
    }
    const sim::SpeedPreset speeds[] = { sim::SpeedPreset::Stiff, sim::SpeedPreset::NearCritical,
                                        sim::SpeedPreset::Marginal };
    const char *speedHelp[] = {
        "Stiff individual speed control, g/(J n_c) = 1: a bird's speed change stays local. [ours]",
        "Near-critical linear speed control, g/(J n_c) = 1e-3 (the paper's value). Speed changes spread "
        "flock-wide; group speed drifts slowly (Cavagna et al. 2022). [model]",
        "Marginal (flat-bottomed) speed control (Cavagna et al. 2022): flock-wide speed correlations "
        "with a stable group speed. Best supported by data. [data + model]" };
    ImGui::Separator();
    ImGui::Text("Speed control");
    for (int k = 0; k < 3; k++) {
        if (ImGui::RadioButton(sim::toString(speeds[k]), mSpeedPreset == speeds[k])) {
            mSpeedPreset = speeds[k];
            applyPresets(false);
        }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", speedHelp[k]);
    }
    ImGui::TextDisabled("Hover a preset for its source and status.");

    ImGui::Separator();
    ImGui::SliderFloat("Time scale", &mTimeScale, 0.05f, 1.0f, "%.2fx", ImGuiSliderFlags_Logarithmic);
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Simulated seconds per real second. A real turn crosses a flock in under a second.");
    ImGui::Checkbox("Paused", &mPaused);
    ImGui::SameLine();
    if (ImGui::Button("New flock")) resetFlock();
    ImGui::SliderInt("Birds", &mN, 100, 1000);
    if (ImGui::IsItemDeactivatedAfterEdit()) resetFlock();

    const char *colorModes[] = { "Turned since attack", "Turn rate", "Speed vs. flock",
                                 "Predator detection", "Plain" };
    int         cm           = (int) mColorMode;
    if (ImGui::Combo("Color", &cm, colorModes, 5)) mColorMode = ColorMode(cm);
    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Turned since attack: angle from the flock's heading when the last attack "
                          "began (white 0, red 60+ deg). Speed vs. flock: blue slower, red faster (+-20%%).");
    ImGui::Checkbox("Grid", &mShowGrid);
    ImGui::SameLine();
    if (ImGui::Button("Reset view")) mSceneRotation = quat();
    ImGui::SliderInt("Show neighbors of bird", &mHighlight, -1, mN - 1);

    ImGui::Separator();
    ImGui::Text("t = %.1f s   polarization %.3f", mFlock->time(), mStats.polarization);
    ImGui::Text("speed %.1f m/s   size L = %.0f m", mStats.meanSpeed, mStats.L);
    ImGui::TextDisabled("Graphs, lower left: heading (green) and speed (blue) correlation vs. "
                        "distance; yellow line = correlation length.");

    if (ImGui::CollapsingHeader("Keyboard shortcuts")) {
        ImGui::TextUnformatted(
            "space  pause / run\n"
            "a      predator attack\n"
            "f      toggle false alarms\n"
            "1 2 3  turning preset A / B / C\n"
            "4 5 6  speed preset S1 / S2 / S3\n"
            "[ ]    slower / faster time\n"
            "c      cycle color mode\n"
            "n      new flock\n"
            "g      toggle grid\n"
            "r      reset view\n"
            "drag   rotate view, wheel: zoom");
    }
    ImGui::PopTextWrapPos();
    ImGui::End();
}

void FlockingApp::drawPredatorWindow()
{
    float uiScale = getWindowContentScale();
    ImGui::SetNextWindowPos(ImVec2(getWindowWidth() * uiScale - 360 * uiScale, 10 * uiScale), ImGuiCond_FirstUseEver);
    ImGui::Begin("Predator", nullptr, ImGuiWindowFlags_AlwaysAutoResize);
    ImGui::PushTextWrapPos(ImGui::GetFontSize() * 24.0f);
    auto &pp = mPredator.params;

    if (ImGui::Button("Attack now (a)")) mPredator.triggerAttack(*mFlock);
    ImGui::SameLine();
    ImGui::Checkbox("Auto", &pp.autoAttack);
    float interval = (float) pp.autoInterval;
    if (ImGui::SliderFloat("Interval (s)", &interval, 5.0f, 60.0f, "%.0f")) pp.autoInterval = interval;

    float radius = (float) pp.detectRadius, hf = (float) pp.headingField, sf = (float) pp.speedField;
    if (ImGui::SliderFloat("Detection radius (m)", &radius, 1.0f, 15.0f, "%.1f")) pp.detectRadius = radius;
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Only birds this close sense the predator directly; the rest learn from their neighbors.");
    if (ImGui::SliderFloat("Turn-away push", &hf, 0.0f, 3.0f, "%.2f x J n_c")) pp.headingField = hf;
    if (ImGui::SliderFloat("Speed-up push", &sf, 0.0f, 1.0f, "%.2f x J n_c")) pp.speedField = sf;
    float rate = (float) pp.falseAlarmRate;
    if (ImGui::SliderFloat("False alarms (/s)", &rate, 0.0f, 2.0f, "%.2f")) pp.falseAlarmRate = rate;
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("A random bird startles as if it saw a predator. In fish schools, false alarms spread like real ones (Rosenthal 2015); a more responsive group pays for it (Poel 2022).");

    ImGui::Separator();
    const auto &r = mPredator.report();
    if (r.attacks == 0)
        ImGui::TextDisabled("No attack yet.");
    else {
        ImGui::Text("Attack %d%s", r.attacks, r.finished ? "" : " (in progress)");
        ImGui::Text("Sensed predator directly: %d birds (%.0f%%)", r.detectedEver, 100.0 * r.detectedEver / mFlock->size());
        ImGui::Text("Turned > 30 deg within 2 s: %.0f%%", 100.0 * r.fractionTurned);
        ImGui::Text("Turn front speed: %.0f m/s", r.frontSpeed);
        ImGui::Text("Lowest polarization: %.2f   pieces: %d", r.minPolarization, r.maxFragments);
        ImGui::Text("Captures: %d", r.captures);
        ImGui::TextDisabled("Model predictions, not measurements of real starlings. Whether critical "
                            "flocks evade better is untested.");
    }
    ImGui::PopTextWrapPos();
    ImGui::End();
}

// Each row compares a live measurement with a published range (docs/SIMULATION_DESIGN.md sec. 5)
void FlockingApp::drawValidationWindow()
{
    float uiScale = getWindowContentScale();
    ImGui::SetNextWindowPos(ImVec2(getWindowWidth() * uiScale - 460 * uiScale, 330 * uiScale), ImGuiCond_FirstUseEver);
    ImGui::Begin("Validation vs. real starling flocks", nullptr, ImGuiWindowFlags_AlwaysAutoResize);

    enum Verdict { Pass, Fail, NA };
    auto row = [](const char *name, const std::string &value, const char *target, Verdict v, const char *source) {
        ImGui::TableNextRow();
        ImGui::TableNextColumn(); ImGui::TextUnformatted(name);
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", source);
        ImGui::TableNextColumn(); ImGui::TextUnformatted(value.c_str());
        ImGui::TableNextColumn(); ImGui::TextUnformatted(target);
        ImGui::TableNextColumn();
        if (v == Pass) ImGui::TextColored(ImVec4(0.3f, 1.0f, 0.3f, 1.0f), "pass");
        else if (v == Fail) ImGui::TextColored(ImVec4(1.0f, 0.35f, 0.35f, 1.0f), "fail");
        else ImGui::TextDisabled("n/a");
    };
    auto in = [](double v, double lo, double hi) { return v >= lo && v <= hi ? Pass : Fail; };
    auto fmt = [](const char *f, double v) { char b[64]; snprintf(b, sizeof b, f, v); return std::string(b); };

    if (ImGui::BeginTable("validation", 4, ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingFixedFit)) {
        ImGui::TableSetupColumn("Observable");
        ImGui::TableSetupColumn("Sim");
        ImGui::TableSetupColumn("Real flocks");
        ImGui::TableSetupColumn("");
        ImGui::TableHeadersRow();
        const sim::Stats &s = mStats;
        row("Polarization", fmt("%.3f", s.polarization), "0.84-0.995", in(s.polarization, 0.84, 0.995),
            "Bialek et al. 2014 Table I; Cavagna, Giardina, Grigera 2018 review");
        row("Interacting neighbors", fmt("%.1f", s.avgNeighbors), "6-8", in(s.avgNeighbors, 6.0, 8.0),
            "Ballerini et al. 2008 (topological interaction)");
        row("Group speed (m/s)", fmt("%.1f", s.meanSpeed), "9.5-14.5", in(s.meanSpeed, 9.5, 14.5),
            "Review 2018: 12 +- 2.5 m/s, independent of flock size");
        row("Speed spread (SD/mean)", fmt("%.2f", s.speedSpread), "0.07-0.2", in(s.speedSpread, 0.07, 0.2),
            "Bialek et al. 2014 Fig. 1a (one snapshot); review figures. Debated.");
        row("Nearest neighbor (m)", fmt("%.2f", s.nearestDist), "0.68-1.51", in(s.nearestDist, 0.68, 1.51),
            "Ballerini et al. 2008 benchmark");
        row("Q_int (paper eq. 1)", fmt("%.3f", s.Qint), "0.006-0.13", in(s.Qint, 0.006, 0.13),
            "Bialek et al. 2014 Table I");
        row("Heading corr. length / L", fmt("%.2f", s.L > 0 ? s.xiDir / s.L : 0), "0.1-0.4",
            s.L > 0 ? in(s.xiDir / s.L, 0.1, 0.4) : NA,
            "Cavagna et al. 2010, 2022. One flock size can't prove 'scale-free': that needs "
            "correlation length proportional to L across sizes (run flocksim_cli at several N).");
        row("Speed corr. length / L", fmt("%.2f", s.L > 0 ? s.xiSp / s.L : 0), "0.1-0.4",
            s.L > 0 ? in(s.xiSp / s.L, 0.1, 0.4) : NA,
            "Cavagna et al. 2010, 2022. Same caveat as above.");
        row("Pieces (fragments)", fmt("%.0f", s.fragments), "1", s.fragments == 1 ? Pass : Fail,
            "Cohesive flocks");

        const auto &r = mPredator.report();
        bool haveTurn = r.attacks > 0 && r.frontSpeed > 0;
        row("Turn front speed (m/s)", haveTurn ? fmt("%.0f", r.frontSpeed) : "-", "20-40",
            haveTurn ? in(r.frontSpeed, 20.0, 40.0) : NA,
            "Attanasi et al., Nature Physics 2014 (published values; the preprint lists 9.4-21.3). "
            "Measured from the last attack; noisy.");
        bool haveProbe = mProbeResult.valid && !mProbeResult.minNormalized.empty();
        double minC    = haveProbe ? mProbeResult.minNormalized[0] : 0;
        double half    = haveProbe ? mProbeResult.halfLife[0] : 0;
        row("Spontaneous wobbles", haveProbe ? (minC < -0.1 ? "oscillate" : "overdamped") : "measuring",
            "overdamped", haveProbe ? (minC < -0.1 ? Fail : Pass) : NA,
            "Cavagna et al., PRL 2026: no spin-wave peaks in real flocks (k ~ 1.5 /m). Needs ~20 s "
            "without attacks.");
        row("Wobble half-life (s)", haveProbe && half > 0 ? fmt("%.2f", half) : "-", "~0.2-0.4",
            haveProbe && half > 0 ? in(half, 0.15, 0.5) : NA,
            "Read off Cavagna et al. PRL 2026 Fig. S3c (omega_L ~ 2-3 /s at k ~ 1.5-2 /m); "
            "Mora et al. 2016: relaxation 0.1-0.5 s. Approximate.");
        ImGui::EndTable();
    }
    ImGui::TextDisabled("Hover a row for its source. Known mismatch not checked: real flocks are flat\n"
                        "(axes about 1 : 2.8 : 5.6); model flocks are round.");
    ImGui::End();
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
