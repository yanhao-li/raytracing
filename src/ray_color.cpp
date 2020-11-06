#include "ray_color.h"

using namespace Eigen;

Vector3d ray_color(const Scene &scene, const Ray &ray, const Object &obj, const Intersection &hit, int max_bounce) {
	// Material for hit object
	const Material &mat = obj.material;

	// Ambient light contribution
	Vector3d ambient_color = obj.material.ambient_color.array() * scene.ambient_light.array();

	// Punctual lights contribution (direct lighting)
	Vector3d lights_color(0, 0, 0);
	for (const Light &light : scene.lights) {
		Vector3d Li = (light.position - hit.position).normalized();
		Vector3d N = hit.normal;

		// TODO: (Assignment 2, shadow rays)

		// Diffuse contribution
		Vector3d diffuse = mat.diffuse_color * std::max(Li.dot(N), 0.0);

		// TODO: (Assignment 2, specular contribution)
		Vector3d specular(0, 0, 0);

		// Attenuate lights according to the squared distance to the lights
		Vector3d D = light.position - hit.position;
		lights_color += (diffuse + specular).cwiseProduct(light.intensity) /  D.squaredNorm();
	}

	// TODO: (Assignment 2, reflected ray)
	Vector3d reflection_color(0, 0, 0);

	// TODO: (Assignment 2, refracted ray)
	Vector3d refraction_color(0, 0, 0);

	// Rendering equation
	Vector3d C = ambient_color + lights_color + reflection_color + refraction_color;

	return C;
}

bool is_light_visible(const Scene &scene, const Ray &ray, const Light &light) {
	// TODO: (Assignment 2, shadow ray)
	return true;
}