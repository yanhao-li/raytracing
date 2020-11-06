#include "types.h"
#include "ray_color.h"

// Function declaration here (could be put in a header file)
Object * find_nearest_object(const Scene &scene, const Ray &ray, Intersection &closest_hit);
Vector3d shoot_ray(const Scene &scene, const Ray &ray, int max_bounce);