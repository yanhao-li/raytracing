// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"
#include "utils.h"
#include "types.h"
#include "load.h"
#include "raytracing.h"
#include <chrono>

// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;

void render_scene(const Scene &scene) {
	std::chrono::steady_clock::time_point timer_begin = std::chrono::steady_clock::now();
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
	std::chrono::steady_clock::time_point timer_end = std::chrono::steady_clock::now();
	std::cout << "Time cost = " << std::chrono::duration_cast<std::chrono::seconds>(timer_end - timer_begin).count() << " seconds" << std::endl;
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
