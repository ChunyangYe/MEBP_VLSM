#include "TriangleInterface.h"

/*bnd_pts����ʱ��˳��洢�ı߽��Ķ�ά���� n X 2��area_threshold���������ε�������Լ����
	pts���ɵ��������񶥵����꣬ע��ǰn����Ȼ��bnd_pts����һ�£��������¼���ĵ㣬
	FV�洢����ÿ�������������id��nf X 3*/
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


/*bnd_pts����ʱ��˳��洢�ı߽��Ķ�ά���� n X 2,ע��bnd_pts�ж����߽��ߣ������ش洢��
	E�洢���������߽��ϵıߵ����ӹ�ϵ������E��0,1;1,2;2,0����ʾ�������߷ֱ���0-1,1-2,2-0��
	hole�洢���Ƕ������꣨һ�����ö��ڲ�һ��������������ڲ��������������������������scaffold��
	pts���ɵ��������񶥵����꣬ע��ǰn����Ȼ��bnd_pts����һ�£��������¼���ĵ㣬
	FV�洢����ÿ�������������id��nf X 3*/

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
