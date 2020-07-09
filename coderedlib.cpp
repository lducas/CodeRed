#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <set>
#include <bitset>
#include <assert.h>
#include <time.h>
#include <unistd.h>

using namespace std;

#define max_n maxn   // Hardcode maximal code length as a compilation parameter maxn

typedef bitset<max_n> binvec;   
typedef vector<binvec> binmat;

int skip = 0;

size_t n, k;        // The length n and dimension k of the code.

binmat B;           // The basis of the code
binmat E;           // The epipodal matrix
binmat P;           // The cumulative projector matrix P[i] = &_{j<i} ~ B[j] (has length k+1)

vector<size_t> l; // The profile of the basis (list of epipodal length)

inline int64_t popcnt(binvec& t)
{
    int64_t ham1 = 0;
    int64_t ham2 = 0;
    uint64_t * t_ = (uint64_t *) &t;

    for (int j = 0; j < maxn/64; j+=4)
    {
        ham1 +=__builtin_popcountll(t_[j+0]);
        ham2 +=__builtin_popcountll(t_[j+1]);
        ham1 +=__builtin_popcountll(t_[j+2]);
        ham2 +=__builtin_popcountll(t_[j+3]);
    }
return ham1 + ham2;
}

inline int64_t AND_popcnt(binvec& t, binvec&  e)
{
    int64_t ham1 = 0;
    int64_t ham2 = 0;
    uint64_t * t_ = (uint64_t *) &t;
    uint64_t * e_ = (uint64_t *) &e;

    for (int j = 0; j < maxn/64; j+=4)
    {
        ham1 +=__builtin_popcountll(t_[j+0] & e_[j+0]);
        ham2 +=__builtin_popcountll(t_[j+1] & e_[j+1]);
        ham1 +=__builtin_popcountll(t_[j+2] & e_[j+2]);
        ham2 +=__builtin_popcountll(t_[j+3] & e_[j+3]);
    }
return ham1 + ham2;
}

// Update the epipodal vectors from beg to end, assuming it is up to date up to beg already.
void UpdateEP(size_t beg, size_t end)
{
    assert(beg <= end);
    assert(end <= k);

    P[0].set();
    for (int i = beg; i < end; ++i)
    {        
        E[i] = B[i] & P[i];
        l[i] = E[i].count();
        P[i+1] = P[i] & ~B[i];
    }
}

void UpdateEP()
{
    UpdateEP(0, k);
}



// Update the epipodal vectors from beg to end,
void SizeRedBasis(size_t beg, size_t end)
{
    assert(beg <= end);
    assert(end <= k);

    for (int j = end-1; j >= (int) beg; --j)
    {
        for (int i = j-1; i >= (int) beg; --i)
        {
            if ((B[j] & P[i]).count() > ((B[j]^B[i]) & P[i]).count()) B[j] ^= B[i];
        }
    }
}

// Apply a random transformation on the basis
void Randomize(bool light=true)
{
    size_t steps = light ? 3*k : k*k;

    for (size_t t = 0; t < steps; ++t)
    {
        size_t i = rand() % k;
        size_t j = rand() % k;
        if (i==j) continue;
        B[i] ^= B[j];
    }
    UpdateEP(0, k);
}

// put the basis in systematic form, according to a random information set.
void Systematize()
{
    for (int i = 0; i < k; ++i)
    {
        size_t pivot = rand() % n;
        while (!B[i][pivot]) pivot = rand() % n;
        for (int j = 0; j < k; ++j)
        {
            if (i==j) continue;
            if (B[j][pivot]) B[j] ^= B[i];
        }
    }
    UpdateEP(0, k);
}

void EpiSort()
{
    binvec p;
    p.set();

    for (int i = 0; i < k; ++i)
    {
        size_t best_w=n, best_j=i;
        for (int j = i; j < k; ++j)
        {
            int w= (p&B[j]).count();

            if (w < best_w)
            {
                best_w = w;
                best_j = j;
            }
        }
        if (i != best_j) swap(B[i], B[best_j]);
        p &= ~B[i];
    }

    UpdateEP(0, k);
}




// Put the basis B into semi-systematic form, only permuting its Epipodal matrix E.
// E and P are not maintained during the computation but simply recomputed at the end.
void SemiSystematize()
{
    SizeRedBasis(0, k);
    for (int i = 0; i < k; ++i)
    {
        for (int j = 0; j < k-1; ++j)
        {
            if ((l[j]==1) & (l[j+1]>1))
            {
                swap(l[j], l[j+1]);
                swap(B[j], B[j+1]);
            }
        }
    }
    UpdateEP();
}

void export_mat(/*output*/ char* M_, /*input*/ binmat& M)
{
    size_t i = 0;
    for (auto& v : M)
    {
        for (size_t k = 0; k < n; ++k)
        {
            M_[i] = v[k];
            i++;
        }
    }
}


void LLL(size_t beg, size_t end)
{
    assert(end <= k);
    size_t i = beg;
    binvec p;

    // Loop invariant: the basis is LLL-reduced from beg to i.
    while(i+1 < end)
    {
        // define the projection
        p = P[i];

        // Local size-reduction
        if (((B[i+1]^B[i]) & p).count() < ((B[i+1]) & p).count()) B[i+1] ^= B[i];

        //Lovasz condition
        if ((B[i+1] & p).count() < (B[i] & p).count())
        {
            swap(B[i+1], B[i]);

            // Update auxiliary data
            E[i] = B[i] & P[i];
            l[i] = E[i].count();
            P[i+1] = P[i] & ~B[i];

            E[i+1] = B[i+1] & P[i+1];
            l[i+1] = E[i+1].count();
            P[i+2] = P[i+1] & ~B[i+1];

            if (i > beg) 
            {
                --i;
                continue;
            }
        }
        ++i;
    }
}



void KillTwos()
{
    for (int i = 0; i < k; ++i)
    {
        if (l[i] != 2) continue;
        for (int j = i+1; j < k; ++j)
            {
                if (l[j] != 2) continue;
                if ((B[j] & P[i]).count() != 3) continue;
                swap(B[i], B[j]);
                UpdateEP(i, j+1);
                LLL(i+1, k);
                break;
            }
    }
    SizeRedBasis(0, k);
}


// set-up auxiliary data for enumeration 
// First and last element is only a helper for streamlining iteration.
vector<int> start(size_t beg, size_t end, size_t w)
{
    vector<int> res;
    res.push_back(beg-1);
    for (size_t i = beg; i < beg+w; ++i) res.push_back(-1);
    res.push_back(end);
    return res;
}

// An helper function to enumerate targets of weight w as follows.
// example enumerating 3 choose 5
// 0 1 2 3 4 5
// 01 02 12 03 13 23 04 14 24 34
// 012 013 023 123 014 024 124 034 134 234
inline bool next(binvec& t, vector<int>& e)
{
    for (size_t i = 1; i < e.size()-1; ++i)
    {
        if (e[i] >= 0)
        {
            // clear codeword from target
            t ^= B[e[i]];
            ++e[i];
            if (i > 1) e[i] += skip;    // Only search a fraction of the space for other indices
        }
        else
        {
            e[i] = e[i-1]+1;
        }

        if ((e[i] < e[i+1]) | ((e[i+1] < 0) & (e[i] < e[e.size()-1])) )
        {
            // add the next codeword
            t ^= B[e[i]];
            return true;            
        }
        else
        {
            // reset the coordinate 
            e[i] = e[i-1]+1;
            // add the codeword
            t ^= B[e[i]];
            // move on to the next coordinate
        }
    }
    return false;
}

void TestEnum(int p, int w)
{
    vector<int> enumerator = start(p, k, w);
    binvec t;
    t.reset();

    while(next(t, enumerator))
    {
        cerr << t << endl;
    }

}

bool LB(binvec& tt, size_t w2, int goal_w, uint64_t* stats)
{
    binvec t = tt;
    vector<int> enumerator = start(0, k, w2);

    // If no goal set, just return the best visited solution
    int best_w = goal_w > 0 ? goal_w + 1 : tt.count();
    if (best_w==0) best_w=n;

    while(next(t, enumerator))
    {
        size_t w = popcnt(t);
        if (stats) stats[w]++;
        if (w >= best_w) continue;
        if (w == 0) continue;
        tt = t;
        if (goal_w > 0) return true;
        best_w = w;
    }
    return (goal_w==0);
}


// Size-reduce the target word t with respect to a (segment) of B
inline void SizeRed(binvec& t, size_t beg, size_t end)
{
    for (int i = end-1; i >= (int) beg; --i)
    {
        // This is the most critical peace: helping the compiler
        // For some reason using (t & P[i]).count() gets slow for n > 1024.
        int64_t ham = (t, E[i]).count();
        if (2*ham > l[i]) t ^= B[i];
    }
}



// Making the critical data contiguous
vector<binvec> stream_SR; 
inline void StreamSizeRed(array<binvec,4>& ts, size_t k1)
{
    for (int i = 0; i < k1; ++i)
    {
        // This is the most critical loop: helping the compiler.
        // For some reason using (t & P[i]).count() gets slow for n > 1024.
        for (int j = 0; j < 4; ++j)
        {
            int64_t ham = AND_popcnt(ts[j], stream_SR[2*i]);
            if (2*ham > l[k1 -1 - i]) ts[j] ^= stream_SR[2*i+1];
        }
    }
}

// Lee-Brickell-Babai
bool LBB(binvec& tt, size_t k1, size_t w2, int goal_w, uint64_t* stats)
{
    binvec t = tt;
    array<binvec,4> ts;
    vector<int> enumerator = start(k1, k, w2);
    stream_SR.resize(2 * k1);
    for (int i = 0; i < k1; ++i)
    {
        stream_SR[2*i] = E[k1 - 1 - i];
        stream_SR[2*i+1] = B[k1 - 1 - i];
    }


    // If no goal set, just return the best visited solution
    int best_w = goal_w > 0 ? goal_w + 1 : tt.count();
    if (best_w==0) best_w=n;

    bool notover = true;
    while(notover)
    {

        notover &= next(t, enumerator);
        ts[0] = t;
        notover &= next(t, enumerator);
        ts[1] = t;
        notover &= next(t, enumerator);
        ts[2] = t;
        notover &= next(t, enumerator);
        ts[3] = t;

        StreamSizeRed(ts, k1);

        for (int i = 0; i < 4; ++i)
        {
            size_t w = popcnt(ts[i]);
            if (stats) stats[w]++;        
            if (w >= best_w) continue;
            if (w == 0) continue;
            tt = ts[i];
            if (goal_w > 0) return true;
            best_w = w;
        }
    }
    return (goal_w==0);
}




extern "C" 
{
    void _setup(/*input*/ size_t k_, size_t n_, char* B_, long seed=0)
    {
        k = k_;
        n = n_;
        assert(n <= max_n);
        if (seed==0) seed = time(NULL)+99997*getpid()+123*clock();
        srand(seed);

        B.clear();
        P.clear();
        E.clear();
        l.clear();

        binvec v, zero;
        zero.reset();
        
        l.resize(k);
        P.push_back(zero);

        size_t i = 0;
        for (size_t j = 0; j < k; ++j)
        {
            v.reset();
            for (size_t k = 0; k < n; ++k)
            {
                v[k] = B_[i];
                i++;
            }
            B.push_back(v);
            E.push_back(zero);
            P.push_back(zero);
        }
        UpdateEP();    
    }

    void _export_all(/*output*/ char* B_, char* E_, char* P_, long* l_)
    {
        export_mat(B_, B);
        export_mat(E_, E);
        export_mat(P_, P);

        for (int i = 0; i < k; ++i)
        {
            l_[i] = l[i];
        }
    }

    void _LLL()
    {
        LLL(0, k);
    }

    bool _LBB(char* tt_, size_t k1, size_t w2, int goal_w, uint64_t* stats)
    {
        binvec tt;
        for (int i = 0; i < n; ++i) tt[i] = tt_[i];
        bool res = LBB(tt, k1, w2, goal_w, stats);
        for (int i = 0; i < n; ++i) tt_[i] = tt[i];
        return res;
    }

    bool _LB(char* tt_, size_t w2, int goal_w, uint64_t* stats)
    {
        binvec tt;
        for (int i = 0; i < n; ++i) tt[i] = tt_[i];
        bool res = LB(tt, w2, goal_w, stats);
        for (int i = 0; i < n; ++i) tt_[i] = tt[i];
        return res;
    }

    void _TestEnum(int p, int w)
    {
        TestEnum(p, w);
    }

    void _SizeRedBasis()
    {
        SizeRedBasis(0, k);
    }
    void _SizeRed(char* tt_)
    {
        binvec tt;
        for (int i = 0; i < n; ++i) tt[i] = tt_[i];
        SizeRed(tt, 0, k);
        for (int i = 0; i < n; ++i) tt_[i] = tt[i];
    }

    void _Systematize()
    {
        Systematize();
    }

    void _EpiSort()
    {
        EpiSort();
    }


    void _SemiSystematize()
    {
        SemiSystematize();
    }

    void _KillTwos()
    {
        KillTwos();
    }

    void _Randomize(int light=1)
    {
        Randomize(light);
    }
    void _set_skip(int skip_)
    {
        skip = skip_;
    }
}