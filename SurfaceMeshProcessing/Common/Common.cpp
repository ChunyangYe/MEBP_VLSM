#include "Common.h"

//使得cut成为有序点列
void orderTheBoundary(vector<list<int>>& order_boundary, const vector<int>& nextlocation)
{
	int n_es = nextlocation.size() / 2;
	set<int> isused;
	for (int i = 0; i < n_es; i++)
	{
		isused.insert(i);
	}

	int s_1 = 0;
	int s_2 = 1;

	//bool flag = true;
	while (!isused.empty())
	{
		s_1 = *(isused.begin()) * 2;
		s_2 = s_1 + 1;
		list<int> list_1, list_2;
		if (growFromP(s_1, isused, list_1, nextlocation))
		{
			growFromP(s_2, isused, list_2, nextlocation);
			for (auto it = list_2.begin(); it != list_2.end(); it++)
			{
				list_1.push_front(*it);
			}
		}
		order_boundary.emplace_back(list_1);
	}

}
/*加入判断边界线是否闭合,return true means open curve; return false means closed curve.
In closed case, the first and last item of order_boundary are same.*/
bool growFromP(int p, set<int>& isused, list<int>& order_boundary, const vector<int>& nextlocation)
{
	int loc = nextlocation[p];
	isused.erase(p / 2);
	if (loc == -1)
	{
		order_boundary.push_back(p);
		return true;
	}
	else
	{
		order_boundary.push_back(loc);
	}

	while (true)
	{
		isused.erase(loc / 2);
		if (loc % 2 == 1)
		{
			if ((loc - 1) == p)
			{
				order_boundary.push_back(nextlocation[p]);
				return false;
			}

			if (nextlocation[loc - 1] != -1)
			{
				order_boundary.push_back(nextlocation[loc - 1]);
				loc = nextlocation[loc - 1];
			}
			else
			{
				order_boundary.push_back(loc - 1);
				break;
			}
		}
		else
		{
			if ((loc + 1) == p)
			{
				order_boundary.push_back(nextlocation[p]);
				return false;
			}
			if (nextlocation[loc + 1] != -1)
			{
				order_boundary.push_back(nextlocation[loc + 1]);
				loc = nextlocation[loc + 1];
			}
			else
			{
				order_boundary.push_back(loc + 1);
				break;
			}
		}


	}

	return true;
}

void Tutte(int V_N, const Eigen::MatrixXi& F, const Eigen::VectorXi& bnd, const Eigen::MatrixXd& bnd_uv, Eigen::MatrixXd & uv_init)
{
	int F_N = F.rows();

	uv_init.resize(V_N, 2);
	std::set<int> bound_ids;
	for (size_t i = 0; i < bnd.rows(); i++)
	{
		bound_ids.insert(bnd(i));
		uv_init.row(bnd(i)) << bnd_uv(i, 0), bnd_uv(i, 1);
	}

	std::vector<std::set<int>> VV_tmp;
	VV_tmp.resize(V_N);
	for (size_t i = 0; i < F_N; i++)
	{
		int vid[3];

		for (size_t j = 0; j < F.cols(); j++)
		{
			vid[j] = F(i, j);
		}
		VV_tmp[vid[0]].insert(vid[1]);
		VV_tmp[vid[0]].insert(vid[2]);

		VV_tmp[vid[1]].insert(vid[0]);
		VV_tmp[vid[1]].insert(vid[2]);

		VV_tmp[vid[2]].insert(vid[0]);
		VV_tmp[vid[2]].insert(vid[1]);
	}

	PardisoSolver pardiso;
	vector<double> pardiso_tu;
	vector<double> pardiso_tv;

	pardiso.ia.reserve(V_N + 1);
	pardiso.ja.reserve(8 * V_N);
	pardiso.a.reserve(8 * V_N);
	pardiso_tu.resize(V_N, 0.0);
	pardiso_tv.resize(V_N, 0.0);

	for (size_t i = 0; i < V_N; i++)
	{
		pardiso.ia.push_back(pardiso.ja.size());

		if (bound_ids.count(i) > 0)
		{
			pardiso.ja.push_back(i);
			pardiso.a.push_back(1.0);

			pardiso_tu[i] = uv_init(i, 0);
			pardiso_tv[i] = uv_init(i, 1);

		}
		else
		{
			pardiso.ja.push_back(i);
			pardiso.a.push_back(VV_tmp[i].size());
			vector<int> row_id;
			row_id.reserve(VV_tmp[i].size());
			double bu = 0.0; double bv = 0.0;
			for (auto &vv_id : VV_tmp[i])
			{
				if (bound_ids.count(vv_id) > 0)
				{
					bu += uv_init(vv_id, 0);
					bv += uv_init(vv_id, 1);
				}
				else
				{
					if (vv_id > i)
					{
						row_id.push_back(vv_id);
					}
				}
			}
			sort(row_id.begin(), row_id.end(), less<int>());
			for (size_t j = 0; j < row_id.size(); j++)
			{
				pardiso.ja.push_back(row_id[j]);
				pardiso.a.push_back(-1.0);
			}
			pardiso_tu[i] = bu;
			pardiso_tv[i] = bv;
		}
	}
	pardiso.ia.push_back(pardiso.ja.size());

	pardiso.nnz = pardiso.ja.size();
	pardiso.num = V_N;
	pardiso.rhs = pardiso_tu;

	clock_t t1, t2;
	t1 = clock();
	pardiso.pardiso_init();
	t2 = clock();
	double t_init = (t2 - t1) / 1000.0;

	pardiso.factorize();
	t1 = clock();
	double t_factorize= (t1 - t2) / 1000.0;

	pardiso.pardiso_solver();
	t2 = clock();
	double t_solve = (t2 - t1) / 1000.0;
	printf("pardiso_init: %f ------//pardiso_factorize: %f ------//pardiso_solve: %f ------\n", t_init, t_factorize, t_solve);
	LOG(INFO) << "pardiso_init: "<< t_init <<" ------//pardiso_factorize: "<< t_factorize <<" ------//pardiso_solve:" << t_solve;
	for (size_t i = 0; i < V_N; i++)
	{
		uv_init(i, 0) = pardiso.result[i];
	}

	pardiso.rhs = pardiso_tv;
	pardiso.pardiso_solver();
	for (size_t i = 0; i < V_N; i++)
	{
		uv_init(i, 1) = pardiso.result[i];
	}
	pardiso.free_numerical_factorization_memory();
}

void preCalc_pardiso(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, PardisoSolver & pardiso)
{
	int V_N = V.rows();
	int F_N = F.rows();
	pardiso.ia.clear(); pardiso.ia.reserve(2 * V_N + 1);
	pardiso.ja.clear(); pardiso.ja.reserve(8 * V_N);

	std::vector<std::set<int>> VV_tmp;
	VV_tmp.resize(V_N);
	for (size_t i = 0; i < F_N; i++)
	{
		int vid[3];

		for (size_t j = 0; j < F.cols(); j++)
		{
			vid[j] = F(i, j);
		}
		VV_tmp[vid[0]].insert(vid[1]);
		VV_tmp[vid[0]].insert(vid[2]);

		VV_tmp[vid[1]].insert(vid[0]);
		VV_tmp[vid[1]].insert(vid[2]);

		VV_tmp[vid[2]].insert(vid[0]);
		VV_tmp[vid[2]].insert(vid[1]);
	}

	for (int i = 0; i < V_N; i++)
	{
		pardiso.ia.push_back(pardiso.ja.size());
		VV_tmp[i].insert(i);
		vector<int> row_id;
		for (auto&var : VV_tmp[i])
		{
			row_id.push_back(var);
		}

		vector<int>::iterator iter = std::find(row_id.begin(), row_id.end(), i);

		int dd = 0;
		for (int k = std::distance(row_id.begin(), iter); k < row_id.size(); k++)
		{
			pardiso.ja.push_back(row_id[k]);
			++dd;
		}
		for (int k = 0; k < row_id.size(); k++)
		{
			pardiso.ja.push_back(row_id[k] + V_N);
			++dd;
		}
	}
	for (int i = V_N; i < 2 * V_N; i++)
	{
		pardiso.ia.push_back(pardiso.ja.size());
		vector<int> row_id;
		for (auto&var : VV_tmp[i - V_N])
		{
			row_id.push_back(var);
		}
		vector<int>::iterator iter = std::find(row_id.begin(), row_id.end(), i - V_N);

		int dd = 0;
		for (int k = std::distance(row_id.begin(), iter); k < row_id.size(); k++)
		{
			pardiso.ja.push_back(row_id[k] + V_N);
			++dd;
		}
	}
	pardiso.ia.push_back(pardiso.ja.size());
	pardiso.nnz = pardiso.ja.size();
	pardiso.num = 2 * V_N;
}

void map_vertices_to_circle(
	const Eigen::MatrixXd& V,
	const Eigen::VectorXi& bnd,
	Eigen::MatrixXd& UV)
{
	// Get sorted list of boundary vertices
	std::vector<int> map_ij;
	map_ij.resize(V.rows());

	for (int i = 0; i < bnd.size(); i++)
	{
		map_ij[bnd[i]] = i;
	}

	// Map boundary to unit circle
	std::vector<double> len(bnd.size());
	len[0] = 0.;

	for (int i = 1; i < bnd.size(); i++)
	{
		len[i] = len[i - 1] + (V.row(bnd[i - 1]) - V.row(bnd[i])).norm();
	}
	double total_len = len[len.size() - 1] + (V.row(bnd[0]) - V.row(bnd[bnd.size() - 1])).norm();

	UV.resize(bnd.size(), 2);
	for (int i = 0; i < bnd.size(); i++)
	{
		double frac = len[i] * 2. * M_PI / total_len;
		UV.row(map_ij[bnd[i]]) << cos(frac), sin(frac);
	}

}

void writeObj(const Eigen::MatrixXd& V_in, const Eigen::MatrixXi& F_ref,const string& outfile)
{
	printf("write mesh------------------------begining\n");
	ofstream of_obj(outfile, ios::trunc);

	if (V_in.cols() == 3)
	{
		for (size_t vid = 0; vid < V_in.rows(); vid++)
		{
			of_obj  << "v " << fixed << setprecision(10) << V_in(vid, 0) << " " << V_in(vid, 1) << " " << V_in(vid, 2) << endl;
		}
	}
	else if(V_in.cols() == 2)
	{
		for (size_t vid = 0; vid < V_in.rows(); vid++)
		{
			of_obj << "v " << fixed << setprecision(10) << V_in(vid, 0) << " " << V_in(vid, 1) << " ";
			of_obj << 0 << endl;
		}
	}


	for (size_t fi = 0; fi < F_ref.rows(); fi++)
	{
		of_obj << "f " << F_ref(fi, 0) + 1 << " " << F_ref(fi, 1) + 1 << " " << F_ref(fi, 2) + 1 << endl;
	}
	of_obj.close();
	printf("write mesh------------------------finishing\n");
}



//void boundary_loop(const Eigen::MatrixXi &F_ref, std::vector<std::vector<int>>& boundaryloop)
//{
//	std::vector<std::vector<int>> boundaryEdges;
//	std::vector<std::vector<int>> edges;
//	int n_fvs = F_ref.cols();
//
//	for (int it = 0; it < F_ref.rows(); it++)
//	{
//		for (int i = 0; i < n_fvs; i++)
//		{
//			int var = F_ref(it, i);
//			int var_n = F_ref(it, (i + 1) % n_fvs);
//			if (var > var_n) std::swap(var, var_n);
//			std::vector<int> edge(4);
//			edge[0] = var;
//			edge[1] = var_n;
//			edge[2] = it;
//			edge[3] = i;
//			edges.emplace_back(edge);
//		}
//	}
//	std::sort(edges.begin(), edges.end());
//	int i = 1;
//	for (; i < edges.size();)
//	{
//		auto& r1 = edges[i - 1];
//		auto& r2 = edges[i];
//		if ((r1[0] == r2[0]) && (r1[1] == r2[1]))
//		{
//			i += 2;
//		}
//		else
//		{
//			boundaryEdges.emplace_back(edges[i - 1]);
//			i++;
//		}
//	}
//	if (i == edges.size())
//		boundaryEdges.emplace_back(edges.back());
//
//	for (auto&var:boundaryEdges)
//	{
//		var[0] = F_ref(var[2], var[3]);
//		var[1] = F_ref(var[2], (var[3] + 1) % n_fvs);
//	}
//	int ev0 = boundaryEdges.front()[0];
//	int ev1 = boundaryEdges.front()[1];
//
//	vector<int> visited;
//	visited.resize(boundaryEdges.size(), 0);
//	visited[0] = 1;
//	vector<int> loop0;
//	loop0.push_back(ev1);
//	while (ev1!=ev0)
//	{
//		for (int i = 1; i < boundaryEdges.size(); i++)
//		{
//			if (visited[i] == 1)
//				continue;
//			if (boundaryEdges[i][0] == ev1)
//			{
//				visited[i] = 1;
//				ev1 = boundaryEdges[i][1];
//				loop0.push_back(ev1);
//				break;
//			}
//		}
//	}
//	boundaryloop.emplace_back(loop0);
//}

//outside loop is saved in counterclockwise; inner loop is saved in clockwise;
void boundary_loop(const Eigen::MatrixXi &F_ref, std::vector<std::vector<int>>& boundaryloop)
{
	boundaryloop.clear();
	std::vector<std::vector<int>> boundaryEdges;
	std::vector<std::vector<int>> edges;
	int n_fvs = F_ref.cols();
	for (int it = 0; it < F_ref.rows(); it++)
	{
		for (int i = 0; i < n_fvs; i++)
		{
			int var = F_ref(it, i);
			int var_n = F_ref(it, (i + 1) % n_fvs);
			if (var > var_n) std::swap(var, var_n);
			std::vector<int> edge(4);
			edge[0] = var;
			edge[1] = var_n;
			edge[2] = it;
			edge[3] = i;
			edges.emplace_back(edge);
		}
	}
	std::sort(edges.begin(), edges.end());
	int i = 1;
	for (; i < edges.size();)
	{
		auto& r1 = edges[i - 1];
		auto& r2 = edges[i];
		if ((r1[0] == r2[0]) && (r1[1] == r2[1]))
		{
			i += 2;
		}
		else
		{
			boundaryEdges.emplace_back(edges[i - 1]);
			i++;
		}
	}
	if (i == edges.size())
		boundaryEdges.emplace_back(edges.back());

	set<pair<int,int>> boundary_es;
	for (auto&var : boundaryEdges)
	{
		var[0] = F_ref(var[2], var[3]);
		var[1] = F_ref(var[2], (var[3] + 1) % n_fvs);
		boundary_es.emplace(var[0], var[1]);
	}

	while (!boundary_es.empty())
	{
		pair<int,int> var = *(boundary_es.begin());
		int ev0 = var.first;
		int ev1 = var.second;
		boundary_es.erase(var);
		vector<int> loop0;
		loop0.push_back(ev1);
		while (ev1 != ev0)
		{
			for (auto itp = boundary_es.begin(); itp != boundary_es.end(); itp++)
			{
				if (itp->first == ev1)
				{
					ev1 = itp->second;
					loop0.push_back(ev1);
					boundary_es.erase(itp);
					break;
				}
			}
		}
		boundaryloop.emplace_back(loop0);
	}

}

void remove_unreferenced(int V_N, const Eigen::MatrixXi & Fsrc, Eigen::VectorXi& J, Eigen::MatrixXi& PF, Eigen::VectorXi& VI)
{
	Eigen::VectorXi mark;
	mark.setZero(V_N);
	for (size_t i = 0; i < Fsrc.rows(); i++)
	{
		for (size_t j = 0; j < Fsrc.cols(); j++)
		{
			mark(Fsrc(i, j)) = 1;
		}
	}
	int newsize = mark.sum();
	VI.resize(V_N);	
	J.resize(newsize);
	int counts = 0;

	for (size_t i = 0; i < mark.size(); i++)
	{
		if (mark(i) == 1)
		{
			VI(i) = counts;
			J(counts) = i;
			counts++;
		}
		else
		{
			VI(i) = -1;
		}
	}

	PF.resize(Fsrc.rows(), Fsrc.cols());
	for (size_t i = 0; i < Fsrc.rows(); i++)
	{
		for (size_t j = 0; j < Fsrc.cols(); j++)
		{
			PF(i, j) = VI(Fsrc(i, j));
		}
	}

}


double get_smallest_pos_quad_zero(double a, double b, double c)
{
	using namespace std;
	double t1, t2;
	if (std::abs(a) < 1.0e-10)
	{
		a *= 1e6;
		b *= 1e6;
		c *= 1e6;
	}
	if (std::abs(a) > 1.0e-10)
	{
		double delta_in = pow(b, 2) - 4 * a * c;
		if (delta_in <= 0)
		{
			return INFINITY;
		}

		double delta = sqrt(delta_in); // delta >= 0
		if (b >= 0) // avoid subtracting two similar numbers
		{
			double bd = -b - delta;
			t1 = 2 * c / bd;
			t2 = bd / (2 * a);
		}
		else
		{
			double bd = -b + delta;
			t1 = bd / (2 * a);
			t2 = (2 * c) / bd;
		}

		assert(std::isfinite(t1));
		assert(std::isfinite(t2));

		if (a < 0) std::swap(t1, t2); // make t1 > t2
		// return the smaller positive root if it exists, otherwise return infinity
		if (t1 > 0)
		{
			return t2 > 0 ? t2 : t1;
		}
		else
		{
			return INFINITY;
		}
	}
	else
	{
		if (b == 0) return INFINITY; // just to avoid divide-by-zero
		t1 = -c / b;
		return t1 > 0 ? t1 : INFINITY;
	}
}

double newton_equation(const double & a, const double & b, const double & K)
{
	double tt = 1;
	double E_d = pow(a, 2 * tt) + pow(b, 2 * tt) + pow(1 / a, 2 * tt) + pow(1 / b, 2 * tt) - K;
	while (abs(E_d) > 1e-5)
	{
		tt = tt - 1 / (2 * log(a)*pow(a, 2 * tt) + 2 * log(b)* pow(b, 2 * tt) + 2 * log(1 / a)* pow(1 / a, 2 * tt) + 2 * log(1 / b)*pow(1 / b, 2 * tt))*(pow(a, 2 * tt) + pow(b, 2 * tt) + pow(1 / a, 2 * tt) + pow(1 / b, 2 * tt) - K);
		E_d = pow(a, 2 * tt) + pow(b, 2 * tt) + pow(1 / a, 2 * tt) + pow(1 / b, 2 * tt) - K;
	}
	return tt;
}

void test_remove_unreferenced(const Eigen::MatrixXd & V, const Eigen::MatrixXi & Fsrc, Eigen::VectorXi & J, Eigen::MatrixXi & PF, Eigen::VectorXi & VI)
{
	Eigen::VectorXi mark;
	mark.setZero(V.rows());
	for (size_t i = 0; i < Fsrc.rows(); i++)
	{
		for (size_t j = 0; j < Fsrc.cols(); j++)
		{
			mark(Fsrc(i, j)) = 1;
		}
	}
	int newsize = mark.sum();
	VI.resize(V.rows());
	J.resize(newsize);
	int counts = 0;

	int marksize = mark.size() - 1;
	for (int i = marksize; i >=0; i--)
	{
		if (mark(i) == 1)
		{
			VI(i) = counts;
			J(counts) = i;
			counts++;
		}
		else
		{
			VI(i) = -1;
		}
	}

	PF.resize(Fsrc.rows(), Fsrc.cols());
	for (size_t i = 0; i < Fsrc.rows(); i++)
	{
		for (size_t j = 0; j < Fsrc.cols(); j++)
		{
			PF(i, j) = VI(Fsrc(i, j));
		}
	}

}

void graph_color(const Mesh & mesh, std::vector<std::vector<int>>& phase)
{
	int nsizes = mesh.n_vertices();
	std::vector<int> color(nsizes, -1);
	std::vector<int>possible_color;

	int ncolors = 0;
	for (int i = 0; i < nsizes; i++)
	{
		std::fill(possible_color.begin(), possible_color.end(), 1);
		auto vj_h = mesh.vertex_handle(i);
		for (auto itvv = mesh.cvv_begin(vj_h); itvv != mesh.cvv_end(vj_h); itvv++)
		{
			int c = color[itvv->idx()];
			if (c >= 0)
			{
				possible_color[c] = 0;
			}
		}
		int color_id = -1;
		for (auto j = 0; j < possible_color.size(); j++)
		{
			if (possible_color[j] != 0)
			{
				color_id = j;
				break;
			}
		}
		if (color_id < 0)
		{
			color_id = ncolors++;
			possible_color.resize(ncolors);
		}
		color[i] = color_id;
	}
	phase.clear();
	phase.resize(ncolors);
	for (int i = 0; i < nsizes; i++)
	{
		phase[color[i]].push_back(i);
	}

	printf("total colors_size: %d\n", ncolors);
	LOG(INFO) << "total colors: " << ncolors;
	for (size_t i = 0; i < ncolors; i++)
	{
		cout << phase[i].size() << endl;
		LOG(INFO) << "Color "<<i<<" has faces: " << phase[i].size();
	}

}

int preCalc_pardiso(const Eigen::MatrixXd & V, const Eigen::MatrixXi & F)
{
	int nnz = 0;
	int V_N = V.rows();
	int F_N = F.rows();

	std::vector<std::set<int>> VV_tmp;
	VV_tmp.resize(V_N);
	for (size_t i = 0; i < F_N; i++)
	{
		int vid[3];

		for (size_t j = 0; j < F.cols(); j++)
		{
			vid[j] = F(i, j);
		}
		VV_tmp[vid[0]].insert(vid[1]);
		VV_tmp[vid[0]].insert(vid[2]);

		VV_tmp[vid[1]].insert(vid[0]);
		VV_tmp[vid[1]].insert(vid[2]);

		VV_tmp[vid[2]].insert(vid[0]);
		VV_tmp[vid[2]].insert(vid[1]);
	}

	for (int i = 0; i < V_N; i++)
	{
		VV_tmp[i].insert(i);
		vector<int> row_id;
		for (auto&var : VV_tmp[i])
		{
			row_id.push_back(var);
		}
		vector<int>::iterator iter = std::find(row_id.begin(), row_id.end(), i);
		int k = std::distance(row_id.begin(), iter);
		nnz += (3*row_id.size() - 2*k);
	}
	return nnz;
}
