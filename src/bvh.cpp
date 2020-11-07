#include "types.h"

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
	// F => collection of facets (triangles) [[index1, index2, index3], [index4, index5, index6],...]
	// V => collection of vertices [[x1, y1, z1], [x2, y2, z2],...]
	// centroids[i][j] => 
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