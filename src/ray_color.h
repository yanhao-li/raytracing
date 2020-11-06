#include "types.h"

Eigen::Vector3d ray_color(const Scene &scene, const Ray &ray, const Object &object, const Intersection &hit, int max_bounce);
bool is_light_visible(const Scene &scene, const Ray &ray, const Light &light);