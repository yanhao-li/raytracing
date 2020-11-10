#include "types.h"

bool intersect_box(const Ray &ray, const AlignedBox3d &box) {
	// TODO: (Assignment 3)
	//
	// Compute whether the ray intersects the given box.
	// There is no need to set the resulting normal and ray parameter, since
	// we are not testing with the real surface here anyway.
	double tmin = 0;
	double tmax = std::numeric_limits<double>::max();
	Vector3d bottom = box.min();
	Vector3d top = box.max();
	for (int i = 0; i < 3; i++) {
		if (std::abs(ray.direction(i)) < std::numeric_limits<double>::min() && (ray.origin[i] < bottom[i] || ray.origin[i] > top[i])) {
			return false;
		} else {
			double t1 = (bottom[i] - ray.origin[i]) / ray.direction[i];
			double t2 = (top[i] - ray.origin[i]) / ray.direction[i];
			tmin = std::max(std::min(t1, t2), tmin);
			tmax = std::min(std::max(t1, t2), tmax);
		}
	}
	if (tmin <= tmax) return true;
	return false;
}

bool intersect_triangle(const Ray &ray, const Vector3d &a, const Vector3d &b, const Vector3d &c, Intersection &hit) {
	// TODO: (Assignment 3) - Done
	//
	// Compute whether the ray intersects the given triangle.
	// If you have done the parallelogram case, this should be very similar to it.
	Matrix3d M;
	M.col(0) = a - b;
	M.col(1) = a - c;
	M.col(2) = ray.direction;

	Vector3d y;
	y = a - ray.origin;

	Vector3d x = M.colPivHouseholderQr().solve(y);
	double u = x(0), v = x(1), t = x(2);

	if (u > 0 && v > 0 && u + v < 1 && t > 0) {
		hit.ray_param = t;
		hit.position = ray.origin + ray.direction * hit.ray_param;
		hit.normal = (b - a).cross(c - a).normalized();
		return true;
	}

	return false;
}

bool Sphere::intersect(const Ray& ray, Intersection& hit) {
  // TODO: - Done
  //
  // Compute the intersection between the ray and the sphere
  // If the ray hits the sphere, set the result of the intersection in the
  // struct 'hit'

	// vector (e - c)
  Vector3d vec = ray.origin - this->position;  

	double a = ray.direction.dot(ray.direction);    
  double b = 2 * ray.direction.dot(vec);
	double c = (vec).dot(vec) - pow(this->radius, 2);  

	// discriminant, B^2 - 4AC
  double delta = pow(b, 2) - 4 * a * c;

  if (delta >= 0) {
		double t1 = (-b + sqrt(delta)) / (2 * a);
		double t2 = (-b - sqrt(delta)) / (2 * a);

    hit.ray_param = std::min(t1, t2);
    if (hit.ray_param < 0) {
        hit.ray_param = std::max(t1, t2);
    }
    if (hit.ray_param < 0) {
        return false;
    }

		hit.position = ray.origin + hit.ray_param * ray.direction;
		hit.normal = (hit.position - this->position).normalized();

    return true;
  }
  return false;
}

bool Parallelogram::intersect(const Ray &ray, Intersection &hit) {
	// TODO: (Assignment 2)
	return false;
}

bool Mesh::intersect(const Ray &ray, Intersection &closest_hit) {
	// TODO: (Assignment 3) - Done

	bool ret = false;
	closest_hit.ray_param = std::numeric_limits<double>::max();
	// Method (1): Traverse every triangle and return the closest hit.
	// for (int i = 0; i < facets.rows(); i++) {
	// 	Vector3d a =  vertices.row(facets(i, 0));
	// 	Vector3d b = vertices.row(facets(i, 1));
	// 	Vector3d c = vertices.row(facets(i, 2));

	// 	Intersection hit;
	// 	if (intersect_triangle(ray, a, b, c, hit)) {
	// 		if (hit.ray_param < closest_hit.ray_param) {
	// 			closest_hit = hit;
	// 			ret = true;
	// 		}
	// 	};
	// }

	// Method (2): Traverse the BVH tree and test the intersection with a
	// triangles at the leaf nodes that intersects the input ray.

	std::stack<int> s;
	s.push(bvh.root);
	while(!s.empty()) {
		AABBTree::Node &node = bvh.nodes[s.top()];
		s.pop();
		if (intersect_box(ray, node.bbox)) {
			if (node.left == -1 && node.right == -1) {
				Vector3d a = vertices.row(facets(node.triangle, 0)).transpose();
				Vector3d b = vertices.row(facets(node.triangle, 1)).transpose();
				Vector3d c = vertices.row(facets(node.triangle, 2)).transpose();

				Intersection hit;
				if (intersect_triangle(ray, a, b, c, hit)) {
					if (hit.ray_param < closest_hit.ray_param) {
						closest_hit = hit;
						ret = true;
					}
				}
			} else {
				s.push(node.left);
				s.push(node.right);
			}
		}
	}

	return ret;
}