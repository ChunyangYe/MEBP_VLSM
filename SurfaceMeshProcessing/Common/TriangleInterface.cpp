#include "TriangleInterface.h"

/*bnd_pts是逆时针顺序存储的边界点的二维坐标 n X 2，area_threshold生成三角形的最大面积约束，
	pts生成的所有网格顶点坐标，注意前n个依然和bnd_pts保持一致，后面是新加入的点，
	FV存储的是每个面的三个顶点id，nf X 3*/
void triangulate(const Eigen::MatrixXd & bnd_pts, const double & area_threshold, Eigen::MatrixXd & pts, Eigen::MatrixXi & FV)
{
	int p_num = bnd_pts.rows();
	Eigen::MatrixXi E;
	E.resize(p_num, 2);

	for (int i = 0; i < p_num - 1; i++)
	{
		E.row(i) << i, i + 1;
	}
	E.row(p_num - 1) << p_num - 1, 0;

	triangulateio in;

	in.numberofpoints = p_num;
	in.pointlist = (double *)calloc(bnd_pts.size(), sizeof(double));

	for (size_t i = 0; i < p_num; i++)
	{
		in.pointlist[2 * i] = bnd_pts(i, 0);
		in.pointlist[2 * i+1] = bnd_pts(i, 1);
	}

	in.numberofpointattributes = 0;
	in.pointmarkerlist = (int *)calloc(p_num, sizeof(int));
	for (unsigned i = 0; i < p_num; ++i) in.pointmarkerlist[i] = 1;

	in.trianglelist = NULL;
	in.numberoftriangles = 0;
	in.numberofcorners = 0;
	in.numberoftriangleattributes = 0;
	in.triangleattributelist = NULL;

	in.numberofsegments = E.size() ? E.rows() : 0;
	in.segmentlist = (int*)calloc(E.size(), sizeof(int));

	for (size_t i = 0; i < E.rows(); i++)
	{
		in.segmentlist[2 * i] = E(i, 0);
		in.segmentlist[2 * i + 1] = E(i, 1);
	}

	in.segmentmarkerlist = (int*)calloc(E.rows(), sizeof(int));
	for (unsigned i = 0; i < E.rows(); ++i) in.segmentmarkerlist[i] = 1;

	in.numberofholes = 0;

	in.numberofregions = 0;

	triangulateio out;
	out.pointlist = NULL;
	out.trianglelist = NULL;
	out.segmentlist = NULL;
	out.segmentmarkerlist = NULL;
	out.pointmarkerlist = NULL;
	std::string full_flags;

	if (area_threshold == 0.0)
	{
		full_flags = "qYYQpz";
		printf("NO area threshold constraints in the triangulation;\n");
	}
	else
	{
		full_flags = "qYYa" + std::to_string(area_threshold) + "Qpz";
	}

	triangulate(const_cast<char*>(full_flags.c_str()), &in, &out, 0);

	pts.resize(out.numberofpoints, 2);
	FV.resize(out.numberoftriangles, 3);

	for (size_t i = 0; i < out.numberofpoints; i++)
	{
		pts.row(i)<<out.pointlist[2 * i], out.pointlist[2 * i + 1];
	}

	for (size_t i = 0; i < out.numberoftriangles; i++)
	{
		FV.row(i)<< out.trianglelist[3 * i], out.trianglelist[3 * i + 1], out.trianglelist[3 * i + 2];
	}
	// Cleanup in
	free(in.pointlist);
	free(in.pointmarkerlist);
	free(in.segmentlist);
	free(in.segmentmarkerlist);
	// Cleanup out
	free(out.pointlist);
	free(out.trianglelist);
	free(out.segmentlist);
	free(out.segmentmarkerlist);
	free(out.pointmarkerlist);
}


/*bnd_pts是逆时针顺序存储的边界点的二维坐标 n X 2,注意bnd_pts有多条边界线，按条地存储，
	E存储的是上述边界上的边的连接关系，比如E（0,1;1,2;2,0）表示有三条边分别是0-1,1-2,2-0；
	hole存储的是洞的坐标（一个洞用洞内部一点坐标代表），洞内部不生成网格。这个函数用于生成scaffold。
	pts生成的所有网格顶点坐标，注意前n个依然和bnd_pts保持一致，后面是新加入的点，
	FV存储的是每个面的三个顶点id，nf X 3*/

void triangulate(const Eigen::MatrixXd & bnd_pts, const Eigen::MatrixXi & E, const Eigen::MatrixXd & hole, Eigen::MatrixXd & pts, Eigen::MatrixXi & FV)
{
	int p_num = bnd_pts.rows();
	triangulateio in;

	in.numberofpoints = p_num;
	in.pointlist = (double *)calloc(bnd_pts.size(), sizeof(double));

	for (size_t i = 0; i < p_num; i++)
	{
		in.pointlist[2 * i] = bnd_pts(i, 0);
		in.pointlist[2 * i + 1] = bnd_pts(i, 1);
	}

	in.numberofpointattributes = 0;
	in.pointmarkerlist = (int *)calloc(p_num, sizeof(int));
	for (unsigned i = 0; i < p_num; ++i) in.pointmarkerlist[i] = 1;

	in.trianglelist = NULL;
	in.numberoftriangles = 0;
	in.numberofcorners = 0;
	in.numberoftriangleattributes = 0;
	in.triangleattributelist = NULL;

	in.numberofsegments = E.size() ? E.rows() : 0;
	in.segmentlist = (int*)calloc(E.size(), sizeof(int));

	for (size_t i = 0; i < E.rows(); i++)
	{
		in.segmentlist[2 * i] = E(i, 0);
		in.segmentlist[2 * i + 1] = E(i, 1);
	}

	in.segmentmarkerlist = (int*)calloc(E.rows(), sizeof(int));
	for (unsigned i = 0; i < E.rows(); ++i) in.segmentmarkerlist[i] = 1;

	in.numberofholes = hole.rows();
	if (hole.rows() > 0)
	{
		in.holelist = (double*)calloc(hole.size(), sizeof(double));
		for (size_t i = 0; i < hole.rows(); i++)
		{
			in.holelist[2 * i] = hole(i, 0);
			in.holelist[2 * i + 1] = hole(i, 1);
		}
	}
	in.numberofregions = 0;
	triangulateio out;
	out.pointlist = NULL;
	out.trianglelist = NULL;
	out.segmentlist = NULL;
	out.segmentmarkerlist = NULL;
	out.pointmarkerlist = NULL;
	std::string full_flags = "qYYQpz";
	triangulate(const_cast<char*>(full_flags.c_str()), &in, &out, 0);

	pts.resize(out.numberofpoints, 2);
	FV.resize(out.numberoftriangles, 3);

	for (size_t i = 0; i < out.numberofpoints; i++)
	{
		pts.row(i) << out.pointlist[2 * i], out.pointlist[2 * i + 1];
	}

	for (size_t i = 0; i < out.numberoftriangles; i++)
	{
		FV.row(i) << out.trianglelist[3 * i], out.trianglelist[3 * i + 1], out.trianglelist[3 * i + 2];
	}
	// Cleanup in
	free(in.pointlist);
	free(in.pointmarkerlist);
	free(in.segmentlist);
	free(in.segmentmarkerlist);
	if (hole.rows() > 0)
		free(in.holelist);
	// Cleanup out
	free(out.pointlist);
	free(out.trianglelist);
	free(out.segmentlist);
	free(out.segmentmarkerlist);
	free(out.pointmarkerlist);
}


void triangulate(const Eigen::MatrixXd &bnd_pts, const Eigen::MatrixXi &E, const Eigen::MatrixXd &hole, const double & area_threshold, Eigen::MatrixXd & pts, Eigen::MatrixXi & FV)
{
	int p_num = bnd_pts.rows();
	triangulateio in;

	in.numberofpoints = p_num;
	in.pointlist = (double *)calloc(bnd_pts.size(), sizeof(double));

	for (size_t i = 0; i < p_num; i++)
	{
		in.pointlist[2 * i] = bnd_pts(i, 0);
		in.pointlist[2 * i + 1] = bnd_pts(i, 1);
	}

	in.numberofpointattributes = 0;
	in.pointmarkerlist = (int *)calloc(p_num, sizeof(int));
	for (unsigned i = 0; i < p_num; ++i) in.pointmarkerlist[i] = 1;

	in.trianglelist = NULL;
	in.numberoftriangles = 0;
	in.numberofcorners = 0;
	in.numberoftriangleattributes = 0;
	in.triangleattributelist = NULL;

	in.numberofsegments = E.size() ? E.rows() : 0;
	in.segmentlist = (int*)calloc(E.size(), sizeof(int));

	for (size_t i = 0; i < E.rows(); i++)
	{
		in.segmentlist[2 * i] = E(i, 0);
		in.segmentlist[2 * i + 1] = E(i, 1);
	}

	in.segmentmarkerlist = (int*)calloc(E.rows(), sizeof(int));
	for (unsigned i = 0; i < E.rows(); ++i) in.segmentmarkerlist[i] = 1;

	in.numberofholes = hole.rows();
	if (hole.rows() > 0)
	{
		in.holelist = (double*)calloc(hole.size(), sizeof(double));
		for (size_t i = 0; i < hole.rows(); i++)
		{
			in.holelist[2 * i] = hole(i, 0);
			in.holelist[2 * i + 1] = hole(i, 1);
		}
	}
	in.numberofregions = 0;
	triangulateio out;
	out.pointlist = NULL;
	out.trianglelist = NULL;
	out.segmentlist = NULL;
	out.segmentmarkerlist = NULL;
	out.pointmarkerlist = NULL;
	std::string full_flags;
	if (area_threshold == 0.0)
	{
		full_flags = "qYYQpz";
		printf("NO area threshold constraints in the triangulation;\n");
	}
	else
	{
		full_flags = "qYYa" + std::to_string(area_threshold) + "Qpz";
	}

	triangulate(const_cast<char*>(full_flags.c_str()), &in, &out, 0);

	pts.resize(out.numberofpoints, 2);
	FV.resize(out.numberoftriangles, 3);

	for (size_t i = 0; i < out.numberofpoints; i++)
	{
		pts.row(i) << out.pointlist[2 * i], out.pointlist[2 * i + 1];
	}

	for (size_t i = 0; i < out.numberoftriangles; i++)
	{
		FV.row(i) << out.trianglelist[3 * i], out.trianglelist[3 * i + 1], out.trianglelist[3 * i + 2];
	}
	// Cleanup in
	free(in.pointlist);
	free(in.pointmarkerlist);
	free(in.segmentlist);
	free(in.segmentmarkerlist);
	if (hole.rows() > 0)
		free(in.holelist);
	// Cleanup out
	free(out.pointlist);
	free(out.trianglelist);
	free(out.segmentlist);
	free(out.segmentmarkerlist);
	free(out.pointmarkerlist);
}
