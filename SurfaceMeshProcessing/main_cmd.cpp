#include"PolyMap\PolyMap.h"
#include<string>
#include<iostream>
#include<fstream>
#include<io.h>
#include<glog\logging.h>
#include<stack>
using namespace std;

class Solution {
public:
	struct edge {
		edge() {};
		edge(int id, int c) {
			toid = id;
			color = c;
		}
		int toid;
		int color;

	};
	struct node {
		node() {};
		node(int vid_, int color_, int dis_) :vid(vid_), color(color_), dis(dis_) {};
		int vid;
		int color;
		int dis;
		bool operator<(const node&rhs) const
		{
			return dis < rhs.dis;
		}
	};
	vector<int> shortestAlternatingPaths(int n, const vector<vector<int>>& red_edges, const vector<vector<int>>& blue_edges) {
		// int INT_MAX=red_edges.size()+blue_edges.size()+1;
		vector<int> shortpath(n, INT_MAX);
		vector<vector<edge>> edges;
		for (auto&var : red_edges)
		{
			edges[var[0]].emplace_back(var[1], 0);
		}
		for (auto&var : blue_edges)
		{
			edges[var[0]].emplace_back(var[1], 1);
		}
		shortpath[0] = 0;
		vector<int> isok(n, 0);
		priority_queue<node> q;
		q.push(node(0, -1, 0));

		while (!q.empty())
		{
			auto top = q.top();
			q.pop();
			isok[top.vid] = 1;
			for (auto&var : edges[top.vid])
			{
				int dis_tmp = top.dis + 1;
				if (isok[var.toid] == 0 && var.color != top.color&&shortpath[var.toid] > dis_tmp)
				{
					shortpath[var.toid] = dis_tmp;
					q.push(node(var.toid, var.color, dis_tmp));
				}
			}
		}

		for (auto&var : shortpath)
		{
			if (var == INT_MAX)
				var = -1;
		}
		return shortpath;
	}
};
void subset(vector<int>& s, int pos, vector<set<int>>& SS, set<int> cur)
{
	if (pos == s.size())
	{
		SS.push_back(cur);
		return;

	}
	subset(s, pos + 1, SS, cur);
	cur.insert(s[pos]);
	subset(s, pos + 1, SS, cur);
}

template<typename T>
struct cmp
{
	bool operator()(const T& lhs, const T& rhs)
	{
		return lhs > rhs;
	}
};
int main(int argc, char *argv[])
{
	if (argc < 4)
	{

		Solution ss;
		ss.shortestAlternatingPaths(1, {  }, {  });
		return 0;
		Mesh m1, m2;
		OpenMesh::IO::read_mesh(m1, "D1_00449.obj");
		OpenMesh::IO::read_mesh(m2, "D1_00449_comp_para_result.obj");
	
		priority_queue<int> q;
		
		int nv = m1.n_vertices();
		int nf = m1.n_faces();
		ofstream ofs("D1_00449_in_para.obj",ios::trunc);
		for (size_t i = 0; i < nv; i++)
		{
			auto vh = m1.vertex_handle(i);
			auto pos_ = m1.point(vh);
			ofs << "v " << pos_[0] << " " << pos_[1] << " " << pos_[2] << endl;
		}



		{
			OpenMesh::Vec3d para_max, para_min;
			para_max = OpenMesh::Vec3d(-DBL_MAX, -DBL_MAX, -DBL_MAX);
			para_min = OpenMesh::Vec3d(DBL_MAX, DBL_MAX, DBL_MAX);

			for (size_t i = 0; i < nv; i++)
			{
				auto vh = m2.vertex_handle(i);
				auto v0 = m2.point(vh);
				para_max.maximize(v0);
				para_min.minimize(v0);
			}
			double max_size = (para_max-para_min).max();


			//for (size_t i = 0; i < nv; i++)
			//{
			//	auto vh = m2.vertex_handle(i);
			//	auto pos_ = m2.point(vh);
			//	ofs << "v " << pos_[0] << " " << pos_[1] << " " << pos_[2] << endl;
			//}

			for (size_t i = 0; i < nv; i++)
			{
				auto vh = m2.vertex_handle(i);
				auto v_p = m2.point(vh);
				v_p = v_p - para_min;
				v_p /= max_size;
				ofs << "vt " << v_p[0] << " " << v_p[1] << endl;
			}

		}


		for (size_t i = 0; i < nf; i++)
		{
			auto fh = m1.face_handle(i);
			ofs << "f";
			for (auto itfv = m1.fv_begin(fh);  itfv!=m1.fv_end(fh) ; itfv++)
			{
				ofs << " " << itfv->idx() + 1 << "/" << itfv->idx() + 1;
			}
			ofs << endl;
		}
		ofs.close();
		cerr << "exe <mesh_file: *.obj> <method_type: our/outofcore/amgcl> <init_type: tutte/partition>" << endl;
	}
	else
	{
		string file_str = argv[1];
		string meth_str = argv[2];
		string init_str = argv[3];
		google::InitGoogleLogging(file_str.c_str());
		FLAGS_log_dir = ".";

		int meth_type = -1;
		int init_type = -1;
		if (meth_str == "our")
			meth_type = 0;
		else if (meth_str == "outofcore")
			meth_type = 1;
		else if (meth_str == "amgcl")
			meth_type = 2;

		if (init_str == "tutte")
			init_type = 0;
		else if (init_str == "partition")
			init_type = 1;

		PolyMap pm;
		pm.load_mesh(file_str);
		pm.run_methods(meth_type, init_type);
		google::ShutdownGoogleLogging();

	}
	return 0;



}
