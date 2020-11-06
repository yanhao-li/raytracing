// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"
#include "utils.h"
#include "ray_color.h"
#include "types.h"
#include "load.h"

// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;

Mesh::Mesh(const std::string &filename) {
	// Load a mesh from a file (assuming this is a .off file), and create a bvh
	load_off(filename, vertices, facets);
	bvh = AABBTree(vertices, facets);
}

////////////////////////////////////////////////////////////////////////////////
// BVH Implementation
////////////////////////////////////////////////////////////////////////////////

// Bounding box of a triangle
AlignedBox3d bbox_triangle(const Vector3d &a, const Vector3d &b, const Vector3d &c) {
	AlignedBox3d box;
	box.extend(a);
	box.extend(b);
	box.extend(c);
	return box;
}

AABBTree::AABBTree(const MatrixXd &V, const MatrixXi &F) {
	// Compute the centroids of all the triangles in the input mesh
	MatrixXd centroids(F.rows(), V.cols());
	centroids.setZero();
	for (int i = 0; i < F.rows(); ++i) {
		for (int k = 0; k < F.cols(); ++k) {
			centroids.row(i) += V.row(F(i, k));
		}
		centroids.row(i) /= F.cols();
	}

	// TODO: (Assignment 3)

	// Method (1): Top-down approach.
	// Split each set of primitives into 2 sets of roughly equal size,
	// based on sorting the centroids along one direction or another.

	// Method (2): Bottom-up approach.
	// Merge nodes 2 by 2, starting from the leaves of the forest, until only 1 tree is left.
}

////////////////////////////////////////////////////////////////////////////////

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

// -----------------------------------------------------------------------------

bool intersect_triangle(const Ray &ray, const Vector3d &a, const Vector3d &b, const Vector3d &c, Intersection &hit) {
	// TODO: (Assignment 3)
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

bool intersect_box(const Ray &ray, const AlignedBox3d &box) {
	// TODO: (Assignment 3)
	//
	// Compute whether the ray intersects the given box.
	// There is no need to set the resulting normal and ray parameter, since
	// we are not testing with the real surface here anyway.
	return false;
}

bool Mesh::intersect(const Ray &ray, Intersection &closest_hit) {
	// TODO: (Assignment 3)

	// Method (1): Traverse every triangle and return the closest hit.
	for (int i = 0; i < facets.rows(); i++) {
		std::cout << facets.row(i) << std::endl;
	}

	// Method (2): Traverse the BVH tree and test the intersection with a
	// triangles at the leaf nodes that intersects the input ray.

	return false;
}

////////////////////////////////////////////////////////////////////////////////
// Define ray-tracing functions
////////////////////////////////////////////////////////////////////////////////

// Function declaration here (could be put in a header file)
Object * find_nearest_object(const Scene &scene, const Ray &ray, Intersection &closest_hit);
Vector3d shoot_ray(const Scene &scene, const Ray &ray, int max_bounce);

// -----------------------------------------------------------------------------

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

////////////////////////////////////////////////////////////////////////////////

void render_scene(const Scene &scene) {
	std::cout << "Simple ray tracer." << std::endl;

	// int w = 640;
	// int h = 480;
	int w = 8;
	int h = 6;
	MatrixXd R = MatrixXd::Zero(w, h);
	MatrixXd G = MatrixXd::Zero(w, h);
	MatrixXd B = MatrixXd::Zero(w, h);
	MatrixXd A = MatrixXd::Zero(w, h); // Store the alpha mask

	// The camera always points in the direction -z
	// The sensor grid is at a distance 'focal_length' from the camera center,
	// and covers an viewing angle given by 'field_of_view'.
	double aspect_ratio = double(w) / double(h);
	double scale_y = tan(scene.camera.field_of_view / 2) * scene.camera.focal_length; // TODO: Stretch the pixel grid by the proper amount here
	double scale_x = scale_y * aspect_ratio; //

	// The pixel grid through which we shoot rays is at a distance 'focal_length'
	// from the sensor, and is scaled from the canonical [-1,1] in order
	// to produce the target field of view.
	// grid_origin - camera_position
	Vector3d grid_origin(-1 * scale_x, 1 * scale_y, 0);
	Vector3d x_displacement(2.0 / w * scale_x, 0, 0);
	Vector3d y_displacement(0, -2.0 / h * scale_y, 0);

	for (unsigned i = 0; i < w; ++i) {
		std::cout << std::fixed << std::setprecision(2);
		std::cout << "Ray tracing: " << (100.0 * i) / w << "%\r" << std::flush;
		for (unsigned j = 0; j < h; ++j) {
			// TODO: (Assignment 2, depth of field)
			Vector3d shift = grid_origin + (i + 0.5) * x_displacement + (j + 0.5) * y_displacement;
			// Vector3d noise(scene.camera.lens_radius * (1.0f * rand() / RAND_MAX), scene.camera.lens_radius * (1.0f * rand() / RAND_MAX), scene.camera.position[2]);

			// Prepare the ray
			Ray ray;

			if (scene.camera.is_perspective) {
				// Perspective camera
        ray.origin = scene.camera.position;
        ray.direction = shift - scene.camera.position;
			} else {
				// Orthographic camera
				ray.origin = scene.camera.position + Vector3d(shift[0], shift[1], 0);
				ray.direction = Vector3d(0, 0, -1);
			}

			int max_bounce = 5;
			Vector3d C = shoot_ray(scene, ray, max_bounce);
			R(i, j) = C(0);
			G(i, j) = C(1);
			B(i, j) = C(2);
			A(i, j) = 1;
		}
	}

	std::cout << "Ray tracing: 100%  " << std::endl;

	// Save to png
	const std::string filename("raytrace.png");
	write_matrix_to_png(R, G, B, A, filename);
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {
	if (argc < 2) {
		std::cerr << "Usage: " << argv[0] << " scene.json" << std::endl;
		return 1;
	}
	Scene scene = load_scene(argv[1]);
	render_scene(scene);
	return 0;
}
