#include "raytracing.h"

Object* find_nearest_object(const Scene& scene, const Ray& ray,
                            Intersection& closest_hit) {
  int closest_index = -1;
  // TODO: - Done.
  //
  // Find the object in the scene that intersects the ray first
  // The function must return 'nullptr' if no object is hit, otherwise it must
  // return a pointer to the hit object, and set the parameters of the argument
  // 'hit' to their expected values.
	double closest_dist = INT_MAX;

  for (int i = 0; i < scene.objects.size(); i++) {
    Intersection hit;

    if (scene.objects[i]->intersect(ray, hit)) {
   		// calculate distance between hit and ray origin
			double dist = (hit.position - ray.origin).norm();
			if (dist < closest_dist) {
				closest_dist = dist;
				closest_hit = hit;
				closest_index = i;
			}
		};
  }

  if (closest_index < 0) {
    // Return a NULL pointer
    return nullptr;
  } else {
    // Return a pointer to the hit object. Don't forget to set 'closest_hit'
    // accordingly!
    return scene.objects[closest_index].get();
  }
}

Vector3d shoot_ray(const Scene &scene, const Ray &ray, int max_bounce) {
	Intersection hit;
	if (Object * obj = find_nearest_object(scene, ray, hit)) {
		// 'obj' is not null and points to the object of the scene hit by the ray
		return ray_color(scene, ray, *obj, hit, max_bounce);
	} else {
		// 'obj' is null, we must return the background color
		return scene.background_color;
	}
}
