#include "Measure.h"
#include <algorithm>
#include <cmath>
#include <numeric>

namespace sim {

namespace {
const int kMinFragment = 5;   // groups smaller than this don't count as fragments

// first zero crossing (from positive to negative) of a binned function, interpolated
double zeroCrossing(const std::vector<double> &C, const std::vector<int> &count, double binWidth)
{
    int prev = -1;
    for (int b = 0; b < (int) C.size(); b++) {
        if (count[b] == 0) continue;
        if (prev >= 0 && C[prev] > 0 && C[b] <= 0) {
            double r0 = (prev + 0.5) * binWidth, r1 = (b + 0.5) * binWidth;
            return r0 + (r1 - r0) * C[prev] / (C[prev] - C[b]);
        }
        prev = b;
    }
    return 0;
}

int findRoot(std::vector<int> &parent, int a)
{
    while (parent[a] != a) a = parent[a] = parent[parent[a]];
    return a;
}
} // namespace

Stats measure(const Flock &flock, bool correlations, double binWidth)
{
    Stats      st;
    const int  N = flock.size();
    if (N == 0) return st;

    Vec3   sumU, sumX;
    double sumW = 0, sumW2 = 0, sumNb = 0;
    for (int i = 0; i < N; i++) {
        sumU += flock.u[i];
        sumX += flock.x[i];
        sumW += flock.w[i];
        sumW2 += flock.w[i] * flock.w[i];
        sumNb += flock.directed.empty() ? 0 : flock.directed[i].size();
    }
    Vec3   meanU = sumU / N;
    double meanW = sumW / N;
    st.polarization = length(meanU);
    st.meanHeading  = normalize(meanU);
    st.center       = sumX / N;
    st.meanSpeed    = meanW * flock.params().v0;
    st.speedSpread  = std::sqrt(std::max(0.0, sumW2 / N - meanW * meanW)) / meanW;
    st.avgNeighbors = sumNb / N;

    // Q_int = 1/(2 v0^2 N) sum_i 1/n_i sum_{j in N_i} |v_i - v_j|^2   (paper eq. 1)
    if (!flock.directed.empty()) {
        double q = 0;
        for (int i = 0; i < N; i++) {
            if (flock.directed[i].empty()) continue;
            Vec3   vi = flock.u[i] * flock.w[i];
            double s  = 0;
            for (int j : flock.directed[i]) s += length2(vi - flock.u[j] * flock.w[j]);
            q += s / flock.directed[i].size();
        }
        st.Qint = q / (2.0 * N * meanW * meanW);
    }

    // flock size, fragments and (optionally) correlation functions share the pair loop
    double              L2 = 0;
    std::vector<double> nn2(N, 1e300);
    for (int i = 0; i < N; i++)
        for (int j = i + 1; j < N; j++) {
            double d2 = length2(flock.x[i] - flock.x[j]);
            L2        = std::max(L2, d2);
            nn2[i]    = std::min(nn2[i], d2);
            nn2[j]    = std::min(nn2[j], d2);
        }
    st.L = std::sqrt(L2);
    double sumNN = 0;
    for (int i = 0; i < N; i++) sumNN += std::sqrt(nn2[i]);
    st.nearestDist = N > 1 ? sumNN / N : 0;

    // fragments: connected components of the neighbor graph, keeping only links shorter than
    // 3 force-free distances (topological links can otherwise bridge a real gap)
    std::vector<int> parent(N);
    std::iota(parent.begin(), parent.end(), 0);
    double maxLink2 = 9.0 * flock.params().rEq * flock.params().rEq;
    if (!flock.directed.empty())
        for (int i = 0; i < N; i++)
            for (int j : flock.directed[i])
                if (length2(flock.x[i] - flock.x[j]) < maxLink2)
                    parent[findRoot(parent, i)] = findRoot(parent, j);
    std::vector<int> compSize(N, 0);
    for (int i = 0; i < N; i++) compSize[findRoot(parent, i)]++;
    st.fragments = (int) std::count_if(compSize.begin(), compSize.end(),
                                       [](int s) { return s >= kMinFragment; });

    if (correlations) {
        st.binWidth = binWidth;
        int nBins   = std::max(1, (int) std::ceil(st.L / binWidth) + 1);
        st.Cdir.assign(nBins, 0.0);
        st.Csp.assign(nBins, 0.0);
        std::vector<int> count(nBins, 0);
        std::vector<Vec3>   du(N);
        std::vector<double> dw(N);
        for (int i = 0; i < N; i++) {
            du[i] = flock.u[i] - meanU;
            dw[i] = flock.w[i] - meanW;
        }
        for (int i = 0; i < N; i++)
            for (int j = i + 1; j < N; j++) {
                int b = (int) (length(flock.x[i] - flock.x[j]) / binWidth);
                if (b >= nBins) continue;
                st.Cdir[b] += dot(du[i], du[j]);
                st.Csp[b] += dw[i] * dw[j];
                count[b]++;
            }
        for (int b = 0; b < nBins; b++)
            if (count[b]) {
                st.Cdir[b] /= count[b];
                st.Csp[b] /= count[b];
            }
        st.xiDir = zeroCrossing(st.Cdir, count, binWidth);
        st.xiSp  = zeroCrossing(st.Csp, count, binWidth);
    }
    return st;
}

// ---------------------------------------------------------------------------------------------

void TurnTracker::start(const Flock &flock, const Vec3 &origin, double thresholdDeg,
                        const Vec3 &turnDirection)
{
    const int N = flock.size();
    Vec3      sumU;
    for (int i = 0; i < N; i++) sumU += flock.u[i];
    mRefHeading   = normalize(sumU);
    mOrigin       = origin;
    mCosThreshold = std::cos(thresholdDeg * 3.14159265358979 / 180.0);
    mSinThreshold = std::sin(thresholdDeg * 3.14159265358979 / 180.0);
    mTurnDir      = length2(turnDirection) > 0 ? normalize(turnDirection) : Vec3();
    mT0 = mLastTime = flock.time();
    mDist.resize(N);
    mArrival.assign(N, -1.0);
    mPeakRate.assign(N, 0.0);
    for (int i = 0; i < N; i++) mDist[i] = length(flock.x[i] - origin);
    mActive = true;
}

void TurnTracker::update(const Flock &flock)
{
    if (!mActive) return;
    mLastTime = flock.time();
    for (int i = 0; i < flock.size(); i++) {
        mPeakRate[i] = std::max(mPeakRate[i], flock.turnRate(i));
        if (mArrival[i] >= 0) continue;
        bool turned = length2(mTurnDir) > 0 ? dot(flock.u[i], mTurnDir) > mSinThreshold
                                            : dot(flock.u[i], mRefHeading) < mCosThreshold;
        if (turned) mArrival[i] = flock.time() - mT0;
    }
}

TurnTracker::Result TurnTracker::result() const
{
    Result r;
    r.elapsed = mLastTime - mT0;
    std::vector<int> idx;
    for (int i = 0; i < (int) mArrival.size(); i++)
        if (mArrival[i] >= 0) idx.push_back(i);
    r.turned         = (int) idx.size();
    r.fractionTurned = mArrival.empty() ? 0 : double(idx.size()) / mArrival.size();
    if (idx.size() < 5) return r;

    // least-squares fit of distance vs arrival time
    double st = 0, sd = 0, stt = 0, std_ = 0, sdd = 0;
    for (int i : idx) {
        st += mArrival[i];
        sd += mDist[i];
        stt += mArrival[i] * mArrival[i];
        std_ += mArrival[i] * mDist[i];
        sdd += mDist[i] * mDist[i];
    }
    double n = idx.size(), vt = stt - st * st / n, vd = sdd - sd * sd / n, cov = std_ - st * sd / n;
    if (vt > 0) {
        r.frontSpeed = cov / vt;
        r.fitR2      = vd > 0 ? cov * cov / (vt * vd) : 0;
    }

    // attenuation: peak turn rate of the farthest third vs the nearest third (among all birds)
    std::vector<int> all(mDist.size());
    std::iota(all.begin(), all.end(), 0);
    std::sort(all.begin(), all.end(), [&](int a, int b) { return mDist[a] < mDist[b]; });
    size_t third = all.size() / 3;
    double nearRate = 0, farRate = 0;
    for (size_t k = 0; k < third; k++) {
        nearRate += mPeakRate[all[k]];
        farRate += mPeakRate[all[all.size() - 1 - k]];
    }
    r.attenuation = nearRate > 0 ? farRate / nearRate : 0;

    double maxD = 0;
    for (int i : idx) maxD = std::max(maxD, mDist[i]);
    int nb = (int) (maxD / r.profileBin) + 1;
    r.profileTime.assign(nb, 0.0);
    r.profileCount.assign(nb, 0);
    for (int i : idx) {
        int b = (int) (mDist[i] / r.profileBin);
        r.profileTime[b] += mArrival[i];
        r.profileCount[b]++;
    }
    for (int b = 0; b < nb; b++)
        if (r.profileCount[b]) r.profileTime[b] /= r.profileCount[b];

    // count-weighted fit over bins holding at least 5 birds
    double W = 0, bt = 0, bd = 0, btt = 0, btd = 0;
    for (int b = 0; b < nb; b++) {
        if (r.profileCount[b] < 5) continue;
        double wgt = r.profileCount[b], t = r.profileTime[b], d = (b + 0.5) * r.profileBin;
        W += wgt; bt += wgt * t; bd += wgt * d; btt += wgt * t * t; btd += wgt * t * d;
    }
    if (W > 0) {
        double var = btt / W - (bt / W) * (bt / W);
        if (var > 0) r.frontSpeedBinned = (btd / W - (bt / W) * (bd / W)) / var;
    }
    return r;
}

// ---------------------------------------------------------------------------------------------

FluctuationProbe::FluctuationProbe(double sampleInterval, std::vector<double> kValues,
                                   double maxLag, double window)
    : mInterval(sampleInterval), mK(std::move(kValues))
{
    const double s = std::sqrt(0.5);
    mDirs          = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { s, s, 0 }, { s, 0, s }, { 0, s, s } };
    mMaxLagSamples = (size_t) std::lround(maxLag / sampleInterval);
    mWindowSamples = (size_t) std::lround(window / sampleInterval);
    mSeries.resize(mK.size() * mDirs.size());
}

void FluctuationProbe::clear()
{
    for (auto &s : mSeries) s.clear();
}

void FluctuationProbe::sample(const Flock &flock)
{
    const int N = flock.size();
    Vec3      sumU, sumX;
    for (int i = 0; i < N; i++) {
        sumU += flock.u[i];
        sumX += flock.x[i];
    }
    Vec3 meanU = sumU / N, c = sumX / N;
    for (size_t a = 0; a < mK.size(); a++)
        for (size_t d = 0; d < mDirs.size(); d++) {
            Vec3                                 kv = mDirs[d] * mK[a];
            std::array<std::complex<double>, 3> A{};
            for (int i = 0; i < N; i++) {
                std::complex<double> ph = std::polar(1.0, dot(kv, flock.x[i] - c));
                Vec3                 du = flock.u[i] - meanU;
                A[0] += ph * du.x;
                A[1] += ph * du.y;
                A[2] += ph * du.z;
            }
            auto &series = mSeries[a * mDirs.size() + d];
            series.push_back(A);
            if (series.size() > mWindowSamples) series.pop_front();
        }
}

FluctuationProbe::Result FluctuationProbe::result() const
{
    Result r;
    r.k = mK;
    for (size_t a = 0; a < mK.size(); a++) {
        std::vector<double> C(mMaxLagSamples + 1, 0.0);
        std::vector<int>    n(mMaxLagSamples + 1, 0);
        for (size_t d = 0; d < mDirs.size(); d++) {
            const auto &s = mSeries[a * mDirs.size() + d];
            if (s.size() <= mMaxLagSamples * 2) continue;
            for (size_t lag = 0; lag <= mMaxLagSamples; lag++)
                for (size_t t = 0; t + lag < s.size(); t++) {
                    double v = 0;
                    for (int c = 0; c < 3; c++) v += std::real(s[t][c] * std::conj(s[t + lag][c]));
                    C[lag] += v;
                    n[lag]++;
                }
        }
        if (n[0] == 0 || C[0] <= 0) {
            r.minNormalized.push_back(0);
            r.halfLife.push_back(0);
            continue;
        }
        double c0 = C[0] / n[0], minv = 1.0, half = -1;
        for (size_t lag = 1; lag <= mMaxLagSamples; lag++) {
            double v = (C[lag] / n[lag]) / c0;
            minv     = std::min(minv, v);
            if (half < 0 && v < 0.5) half = lag * mInterval;
        }
        r.minNormalized.push_back(minv);
        r.halfLife.push_back(half);
        r.valid = true;
    }
    return r;
}

} // namespace sim
