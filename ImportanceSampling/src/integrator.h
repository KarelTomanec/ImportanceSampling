#pragma once
#include "scene.h"


class Integrator {
public:

	// Computes a set of measurements of the scene lighting
	virtual void Render(const Scene &scene) = 0;
};

class SamplerIntegrator {
public:

	SamplerIntegrator(std::shared_ptr<const Camera> camera, std::shared_ptr<Sampler> sampler)
		: camera(camera), sampler(sampler) {}

	void Render(const Scene& scene);

	virtual Spectrum Li(const RayDifferential& ray, const Scene& scene,
		Sampler& sampler, MemoryArena& arena, int depth = 0) const = 0;

};