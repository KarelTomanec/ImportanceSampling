#pragma once
#include <iostream>
#include <algorithm>
#include <random>

// The cost of sampling - should be measured and set
double costBRDF = 1.0, costLight = 1.0, referenceEfficiency = 1.0;

double INV_PI = 0.31830988618379067154;
double INV_2_PI = 0.15915494309189533577;

// Image
const int screenWidth = 600;
const int screenHeight = 600;
vec3 image[screenWidth * screenHeight]; // computed image
vec3 reference[screenWidth * screenHeight]; // reference image
// alpha of light source sampling, 0 .. BRDF only, 1.0 .. light only
double weight[screenWidth * screenHeight] = { 0 };

// Pseudocolor mapping of a scalar value
// http://www.kennethmoreland.com/color-maps/, CoolWarmFloat33.csv 
float pscols[4 * 33] = { // 33 colors RGB
  0,0.2298057,0.298717966,0.753683153, 0.03125,0.26623388,0.353094838,0.801466763,
  0.0625,0.30386891,0.406535296,0.84495867, 0.09375,0.342804478,0.458757618,0.883725899,
  0.125,0.38301334,0.50941904,0.917387822, 0.15625,0.424369608,0.558148092,0.945619588,
  0.1875,0.46666708,0.604562568,0.968154911, 0.21875,0.509635204,0.648280772,0.98478814,
  0.25,0.552953156,0.688929332,0.995375608, 0.28125,0.596262162,0.726149107,0.999836203,
  0.3125,0.639176211,0.759599947,0.998151185, 0.34375,0.681291281,0.788964712,0.990363227,
  0.375,0.722193294,0.813952739,0.976574709, 0.40625,0.761464949,0.834302879,0.956945269,
  0.4375,0.798691636,0.849786142,0.931688648, 0.46875,0.833466556,0.860207984,0.901068838,
  0.5,0.865395197,0.86541021,0.865395561, 0.53125,0.897787179,0.848937047,0.820880546,
  0.5625,0.924127593,0.827384882,0.774508472, 0.59375,0.944468518,0.800927443,0.726736146,
  0.625,0.958852946,0.769767752,0.678007945, 0.65625,0.96732803,0.734132809,0.628751763,
  0.6875,0.969954137,0.694266682,0.579375448, 0.71875,0.966811177,0.650421156,0.530263762,
  0.75,0.958003065,0.602842431,0.481775914, 0.78125,0.943660866,0.551750968,0.434243684,
  0.8125,0.923944917,0.49730856,0.387970225, 0.84375,0.89904617,0.439559467,0.343229596,
  0.875,0.869186849,0.378313092,0.300267182, 0.90625,0.834620542,0.312874446,0.259301199,
  0.9375,0.795631745,0.24128379,0.220525627, 0.96875,0.752534934,0.157246067,0.184115123,
  1.0,0.705673158,0.01555616,0.150232812
};

// Structure to compute variance
struct SVAR {
	unsigned int cnt; // the number of samples taken to compute statistics
	double mean; // mean value1
	double M2; // sum for variance1
	void Reset() { cnt = 0; mean = M2 = 0; }
	// Statistical support
	SVAR() { Reset(); }
	// add a single sample
	void Update(const double newSampleValue) {
		cnt++;
		double delta = newSampleValue - mean;
		mean += delta / (double)cnt;
		M2 += delta * (newSampleValue - mean);
	}
	// It returns unbiased sample variance (so not for finite population)
	double Evaluate() { return (double)M2 / ((double)cnt - 1); }
};

const int rainbowPSC = 0; // 0 .. use CoolWarm mapping, 1 .. use rainbow color mapping
const int showBargraph = 1; // 0/1 .. dont use/use bargraph on the right for color mapping

enum Show { DIFF, WEIGHT, WEIGHT_PSEUDOCOLOR };
const Show showFlag = WEIGHT_PSEUDOCOLOR;

// Compute random number in range [0,1], uniform distribution
//double drandom() { return (double)rand() / RAND_MAX; }
double drandom() {
	static std::uniform_real_distribution<double> distribution(0.0, 1.0);
	static std::mt19937 generator;
	return distribution(generator);

}



inline float clamp(float val, float low, float high) {
	if (val < low) return low;
	else if (val > high) return high;
	else return val;
	
}


inline float SphericalTheta(const vec3& v) {
	//return acos(clamp(v.z, -1, 1)); // pbrt
	return acos(clamp(v.y, -1, 1)); // funkcni
}

inline float SphericalPhi(const vec3& v) {
	//float p = atan2(v.y, v.x); // pbrt
	float p = atan2(v.z, v.x) - M_PI / 2.0f; // funkcni
	return (p < 0) ? (p + 2.0f * M_PI) : p;
}

struct Distribution1D {

	Distribution1D(const float* f, int n) {
		count = n;
		func = new float[n];
		memcpy(func, f, n * sizeof(float));
		cdf = new float[n + 1];

		cdf[0] = 0.;
		for (int i = 1; i < count + 1; ++i) {
			cdf[i] = cdf[i - 1] + func[i - 1] / n;
		}


		funcInt = cdf[count];
		if (funcInt == 0.0f) {
			for (int i = 1; i < n + 1; ++i) {
				cdf[i] /= float(n);
			}
		}
		else {

			for (int i = 1; i < n + 1; ++i) {
				cdf[i] /= funcInt;
			}
		}

	}
	~Distribution1D() {
		delete[] func;
		delete[] cdf;

	}
	float SampleContinuous(float u, float* pdf, int* off = NULL) const {

		float* ptr = std::lower_bound(cdf, cdf + count + 1, u);
		int offset = max(0, int(ptr - cdf - 1));
		if (off) *off = offset;

		//assert(offset < count);
		//assert(u >= cdf[offset] && u < cdf[offset + 1]);

		float du = (u - cdf[offset]) / (cdf[offset + 1] - cdf[offset]);
		assert(!isnan(du));

		if (pdf) *pdf = func[offset] / funcInt;

		return (offset + du) / count;

	}
private:
	friend struct Distribution2D;
	float* func;
	float *cdf;
	float funcInt;
	int count;

};

struct Distribution2D {

	Distribution2D(const float* func, int nu, int nv) {
		pConditionalV.reserve(nv);
		for (int v = 0; v < nv; ++v) {

			pConditionalV.push_back(new Distribution1D(&func[v * nu], nu));

		}

		std::vector<float> marginalFunc;
		marginalFunc.reserve(nv);
		for (int v = 0; v < nv; ++v) {
			marginalFunc.push_back(pConditionalV[v]->funcInt);
		}
		pMarginal = new Distribution1D(&marginalFunc[0], nv);
	}
	~Distribution2D() {
		delete pMarginal;
		for (uint32_t i = 0; i < pConditionalV.size(); ++i) {
			delete pConditionalV[i];
		}
	}
	void SampleContinuous(float u0, float u1, float uv[2],
		float* pdf) const {
		float pdfs[2];
		int v;
		uv[1] = pMarginal->SampleContinuous(u1, &pdfs[1], &v);
		uv[0] = pConditionalV[v]->SampleContinuous(u0, &pdfs[0]);
		*pdf = pdfs[0] * pdfs[1];

	}
	float Pdf(float u, float v) const {
		int iu = clamp((int)(u * pConditionalV[0]->count), 0,
			pConditionalV[0]->count - 1);
		int iv = clamp((int)(v * pMarginal->count), 0,
			pMarginal->count - 1);
		if (pConditionalV[iv]->funcInt * pMarginal->funcInt == 0.f) return 0.f;
		return (pConditionalV[iv]->func[iu] * pMarginal->func[iv]) / (pConditionalV[iv]->funcInt * pMarginal->funcInt);

	}
private:
	// Distribution2D Private Data
	std::vector<Distribution1D*> pConditionalV;
	Distribution1D* pMarginal;

};



/*
struct Distribution1D {
public:

	//std::vector<float> func;
	//std::vector<float> cdf;
	float* func;
	float* cdf;
	float funcInt;
	float invFuncInt;
	float invCount;
	int count;

public:

	Distribution1D() = default;

	Distribution1D(float* f, int n) {

		func = new float[n];
		cdf = new float[n + 1];
		count = n;

		memcpy(func, f, n * sizeof(float));

		cdf[0] = 0.0f;

		for (int i = 1; i < n + 1; i++)
		{
			cdf[i] = cdf[i - 1] + func[i - 1] / n;
		}

		funcInt = cdf[n];
		invFuncInt = 1.0 / funcInt;
		invCount = 1.0 / count;

		for (int i = 1; i < n + 1; i++)
		{
			cdf[i] *= invFuncInt;
		}
	}

	~Distribution1D() {
		delete[] func;
		delete[] cdf;
	}

	//float sample(float u, float* pdf) {
	//    // find surrouding cdf segments
	//    auto ptr = std::lower_bound(cdf.begin(), cdf.end(), u);
	//    int offset = int(ptr - cdf.begin() - 1);
	//    // return offset along current cdf segment
	//    //std::cout << offset << std::endl;
	//    if (offset < 0 || offset > cdf.size() - 2) return 0;
	//    u = (u - cdf[offset]) / (cdf[offset + 1] - cdf[offset]);
	//    *pdf = func[offset] * invFuncInt;
	//    return offset + u;
	//}

	float sample(float u, float* pdf) {
		// find surrouding cdf segments
		float* ptr = std::lower_bound(cdf, cdf + count + 1, u);
		int offset = int(ptr - cdf - 1);
		// return offset along current cdf segment
		//std::cout << offset << std::endl;
		//if (offset < 0 || offset > count - 2) return 0;
		u = (u - cdf[offset]) / (cdf[offset + 1] - cdf[offset]);
		*pdf = func[offset] * invFuncInt;
		return offset + u;
	}


};*/
