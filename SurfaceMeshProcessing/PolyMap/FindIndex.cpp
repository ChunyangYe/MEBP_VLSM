#include "FindIndex.h"



FindIndex::FindIndex(const Eigen::MatrixXd& s_V_, const Eigen::MatrixXi& s_F_):s_V(s_V_),s_F(s_F_)
{
	triangles.reserve(s_F.rows());

	for (size_t i = 0; i < s_F.rows(); i++)
	{
		Point3 a(s_V(s_F(i, 0), 0), s_V(s_F(i, 0), 1), 0.0);
		Point3 b(s_V(s_F(i, 1), 0), s_V(s_F(i, 1), 1), 0.0);
		Point3 c(s_V(s_F(i, 2), 0), s_V(s_F(i, 2), 1), 0.0);
		triangles.emplace_back(a, b, c);
	}
	// constructs AABB tre
	tree.insert(triangles.begin(), triangles.end());
	//tree.accelerate_distance_queries();
}


FindIndex::~FindIndex()
{
}

int FindIndex::query_point(const Eigen::Vector3d & p_)
{
	// compute closest point and squared distance
	Point3 point_query(p_[0],p_[1],p_[2]);
	auto& pp = tree.any_intersected_primitive(point_query);
	if (pp.has_value())
	{
		int id_ = std::distance(triangles.begin(), pp.get());
		return id_;
	}
	else
	{
		return -1;
	}

}
