#include "types.h"
#include <numeric>

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

	// Top-down approach.
	// Split each set of primitives into 2 sets of roughly equal size,
	// based on sorting the centroids along one direction or another.
	// Find the longest dimention
	std::vector<int> triangles(F.rows()); 
	std::iota(triangles.begin(), triangles.end(), 0);

	// function for build the tree for range [l, r) of triangles
	// return index of root
	std::function<int(int, int, int)> build = [&](int l, int r, int parent) {
		if (r - l == 0) {
			return -1;
		}
		else if (r - l == 1) { // leaf node
			Node node;
			Vector3d a = V.row(F(l, 0));
			Vector3d b = V.row(F(l, 1));
			Vector3d c = V.row(F(l, 2));
			node.bbox = bbox_triangle(a, b, c);
			node.parent = parent;
			node.left = -1;
			node.right = -1;
			node.triangle = l;
			nodes.push_back(node);
			return (int) (nodes.size() - 1);
		}
		else {
			AlignedBox3d cbox;
			// generate the aligned box for all centroids
			for (int i = l; i < r; i++) {
				Vector3d centroid = centroids.row(i).transpose();
				cbox.extend(centroid);
			}

			// use the diagonal vector to find the longest side
			Vector3d diagonal = cbox.diagonal();
			int longest_side = 0;
			for (int i = 1; i < 3; i++) {
				if (diagonal(i) > diagonal(longest_side)) {
					longest_side = i;
				}
			}

			std::sort(triangles.begin() + l, triangles.begin() + r, [&](int i1, int i2) {
				return centroids(i1, longest_side) < centroids(i2, longest_side);
			});

			int midpoint = (l + r) / 2;
			int cur = nodes.size();
			nodes.resize(cur + 1);
			int left = build(l, midpoint, cur);
			int right = build(midpoint, r, cur);
			Node &node = nodes[cur];
			node.left = left;
			node.right = right;
			node.parent = parent;
			node.triangle = -1;
			node.bbox = nodes[node.left].bbox.extend(nodes[node.right].bbox);

			return cur;
		}

	};

	root = build(0, triangles.size(), -1);
}