#pragma once
#include <iostream>
#include <vector>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <Eigen/Core>
#include <Eigen/Geometry>

typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Point_3 Point3;
typedef K::Triangle_3 Triangle;
typedef std::vector<Triangle>::iterator Iterator;
typedef CGAL::AABB_triangle_primitive<K, Iterator> Primitive;
typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;
class FindIndex
{
public:
	FindIndex(const Eigen::MatrixXd& s_V_, const Eigen::MatrixXi& s_F_);
	~FindIndex();

public:
	int query_point(const Eigen::Vector3d& p_);

private:

	std::vector<Triangle> triangles;
	Tree tree;
	const Eigen::MatrixXd& s_V;
	const Eigen::MatrixXi& s_F;
};

