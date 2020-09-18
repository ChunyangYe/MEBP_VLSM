#include"SharedData.h"

SharedData::SharedData() {
	dim = 2;
	mv_num = 0;
	mf_num = 0;
	sf_num = 0;
	sv_num = 0;
	mesh_measure = 0.;
};

void SharedData::update_scaffold()
{
	mv_num = V.rows();
	mf_num = F.rows();

	V_num = w_uv.rows();
	sf_num = s_T.rows();

	sv_num = V_num - mv_num;
	F_num = sf_num + mf_num;

}
void SharedData::adjust_scaf_weight(double new_weight)
{
	scaf_weight = new_weight;
}

void SharedData::handle_mintri()
{
	double min_bnd_edge_len = numeric_limits<double>::infinity();
	int acc_bnd = 0;
	for (int i = 0; i < bnd_sizes.size(); i++)
	{
		int current_size = bnd_sizes[i];

		for (int e = acc_bnd; e < acc_bnd + current_size - 1; e++)
		{
			min_bnd_edge_len = (std::min)(min_bnd_edge_len,
				(w_uv.row(internal_bnd(e)) - w_uv.row(internal_bnd(e + 1))).squaredNorm());
		}
		min_bnd_edge_len = (std::min)(min_bnd_edge_len,
			(w_uv.row(internal_bnd(acc_bnd)) - w_uv.row(internal_bnd(acc_bnd + current_size - 1))).squaredNorm());

		acc_bnd += current_size;
	}
	area_threshold = min_bnd_edge_len / 4.0;
}

void SharedData::update_scaf_jacobian()
{
	handle_mintri();
	scaf_jac.resize(sf_num);
#pragma omp parallel for
	for (int i = 0; i < sf_num; i++)
	{
		int f0 = s_T(i, 0);
		int f1 = s_T(i, 1);
		int f2 = s_T(i, 2);

		Eigen::Vector2d x_(w_uv(f0, 0) - w_uv(f2, 0), w_uv(f0, 1) - w_uv(f2, 1));
		Eigen::Vector2d l_(w_uv(f1, 0) - w_uv(f2, 0), w_uv(f1, 1) - w_uv(f2, 1));

		double area_tri = abs(x_(0)*l_(1) - x_(1)*l_(0))/2.;
		double x1_0, x2_0, y2_0;
		if (area_tri > area_threshold)
		{
			x1_0 = x_.norm();
			x_ /= x1_0;
			Eigen::Vector2d y_(-x_(1), x_(0));
			x2_0 = l_.dot(x_);
			y2_0 = l_.dot(y_);
		}
		else
		{
			//cout << "area too small!!!!!!!!!!!!! " << endl;
			double h = sqrt((2 * area_threshold) / sqrt(3.0));
			x1_0 = h;
			x2_0 = h / 2.0;
			y2_0 = sqrt(3.0)*h / 2.0;
		}
		scaf_jac[i] = Eigen::Vector4d(1.0 / x1_0, -x2_0 / (x1_0*y2_0), 0.0, 1.0 / y2_0);
	}

}

void SharedData::prepare_data_with_partion(int loops_num)
{
	component_sizes.push_back(mesh.n_faces());
	int component_num = component_sizes.size();

	rescale_and_init_area();

	{
		double e_len_sum = 0.0;
		for (const auto& he_h_id : boundary_es_del)
		{
			auto he_h = mesh.halfedge_handle(he_h_id);
			e_len_sum += mesh.calc_edge_length(he_h);
		}
		e_len_avg = e_len_sum / boundary_es_del.size();
		min_tri_area = e_len_avg * e_len_avg*1e-10 / 4.0;
		cout << "e_len_avg: " << e_len_avg << " .min_tri_area: " << min_tri_area << endl;
	}

	init_V_F_VF();

	mv_num = mesh.n_vertices();
	mf_num = mesh.n_faces();
	printf("Vn: %d ; Fn: %d ;\n", mv_num, mf_num);

	init_loops(loops_num);

	if (component_num == 1)//single patch
	{	
		for (auto cur_bnd : boundary_loops)
		{
			internal_bnd.conservativeResize(internal_bnd.rows() + cur_bnd.size());
			internal_bnd.bottomRows(cur_bnd.size()) = Eigen::Map<Eigen::ArrayXi>(cur_bnd.data(), cur_bnd.size());
			bnd_sizes.push_back(cur_bnd.size());
		}
		mf_num = mesh.n_faces();
		mv_num = mesh.n_vertices();
		mesh_measure = 1.0;
	}
	else// multiple patches
	{

	}
	rect_frame_V.resize(0, 0);
	//mesh_improve(true);
	init_src_inv();

}

void SharedData::write_obj_src_uv(const string & outstr)
{

	ofstream of_obj(outstr, ios::trunc);
	of_obj << "mtllib " << model_name << ".mtl\n";

	for (size_t vid = 0; vid < V.rows(); vid++)
	{
		of_obj << "v " << fixed << setprecision(10) << V(vid, 0) << " " << V(vid, 1) << " " << V(vid, 2) << endl;
	}
	for (size_t vid = 0; vid < mv_num; vid++)
	{
		of_obj << "vt " << fixed << setprecision(10) << w_uv(vid, 0) << " " << w_uv(vid, 1) << endl;
	}

	for (size_t fi = 0; fi < mf_num; fi++)
	{
		of_obj << "f " << F(fi, 0) + 1 << "/" << F(fi, 0) + 1 << " " << F(fi, 1) + 1 << "/" << F(fi, 1) + 1 << " " << F(fi, 2) + 1 << "/" << F(fi, 2) + 1 << endl;
	}
	of_obj.close();

}

void SharedData::local_opt_onebyone(vector<std::vector<int>>& phase, int max_iter_num_, double convergence_rate_)
{
	double e_pre_;
	double cur_conv_rate = 1.0;

	printf("Enter the local refinement phase(onebyone_fixed_boundary) with init energy: %f.............\n", Energy_cur);
	LOG(INFO) << "Enter the local refinement phase(onebyone_fixed_boundary) with init energy: " << Energy_cur;
	clock_t t1 = clock();

	for (int loc_refine_ = 0; loc_refine_ < max_iter_num_; loc_refine_++)
	{
		e_pre_ = Energy_cur;
		for (const auto&c_vids : phase)
		{
#pragma omp parallel for
			for (int i = 0; i < c_vids.size(); i++)
			{
				if (boundary_vs.find(c_vids[i]) == boundary_vs.end())
					local_opt_1v(c_vids[i]);
			}
		}
		Energy_cur = calc_energy(w_uv);
		cur_conv_rate = (e_pre_ - Energy_cur) / e_pre_;

		if (cur_conv_rate < convergence_rate_)
		{
			break;
		}
	}
	clock_t t2 = clock();
	energy_record.push_back(Energy_cur);

	printf("Local-Refinement iterations Time: %f ,Energy: %f \n", (t2 - t1) / 1000.0, Energy_cur);
	LOG(INFO) << "Local-Refinement iterations Time: " << (t2 - t1) / 1000.0 << " Energy: " << Energy_cur;

}

double SharedData::calc_energy(const Eigen::MatrixXd & pos_cur)
{
	double e_sum = 0.0;
#pragma omp parallel for reduction(+:e_sum)
	for (int i = 0; i < mf_num; i++)
	{
		int f0 = F(i, 0);
		int f1 = F(i, 1);
		int f2 = F(i, 2);
		double q00 = pos_cur(f0, 0) - pos_cur(f2, 0);
		double q01 = pos_cur(f1, 0) - pos_cur(f2, 0);
		double q10 = pos_cur(f0, 1) - pos_cur(f2, 1);
		double q11 = pos_cur(f1, 1) - pos_cur(f2, 1);
		//p00 = src_inv00[i]; p01 = src_inv01[i]; p10 = src_inv10[i]; p11 = src_inv11[i];
		double p00 = upd_inv00[i];
		double p01 = upd_inv01[i];
		//			double p10 = upd_inv10[i]; 
		double p11 = upd_inv11[i];
		double j00 = p00 * q00;
		double j01 = p01 * q00 + p11 * q01;
		double j10 = p00 * q10;
		double j11 = p01 * q10 + p11 * q11;
		double det = j00 * j11 - j01 * j10;
		double tr = j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11;
		double E_d = src_weight[i] * (tr + tr / (det*det));
		e_sum += E_d;
	}
	return e_sum;
}


void SharedData::local_decrease_1v(int vid, Eigen::Vector2d & p_g_)
{
	p_g_.setZero();
	auto seed_h = mesh.vertex_handle(vid);
	for (auto itvf = mesh.vf_begin(seed_h); itvf != mesh.vf_end(seed_h); itvf++)
	{
		int fid_ = itvf->idx();
		int f0, f1, f2;
		double p00, p01, p10, p11, area_now;
		if (fid_ < mf_num)
		{
			f0 = F(fid_, 0);
			f1 = F(fid_, 1);
			f2 = F(fid_, 2);
			p00 = upd_inv00[fid_];
			p01 = upd_inv01[fid_];
			p10 = upd_inv10[fid_];
			p11 = upd_inv11[fid_];
			area_now = src_weight[fid_];
		}
		else
		{
			//		continue;
			int fid_i = fid_ - mf_num;
			f0 = s_T(fid_i, 0);
			f1 = s_T(fid_i, 1);
			f2 = s_T(fid_i, 2);
			p00 = scaf_jac[fid_i][0];
			p01 = scaf_jac[fid_i][1];
			p10 = scaf_jac[fid_i][2];
			p11 = scaf_jac[fid_i][3];
			area_now = scaf_weight;
		}
		double q00 = w_uv(f0, 0) - w_uv(f2, 0);
		double q01 = w_uv(f1, 0) - w_uv(f2, 0);
		double q10 = w_uv(f0, 1) - w_uv(f2, 1);
		double q11 = w_uv(f1, 1) - w_uv(f2, 1);

		double j00 = p00 * q00 + p10 * q01;
		double j01 = p01 * q00 + p11 * q01;
		double j10 = p00 * q10 + p10 * q11;
		double j11 = p01 * q10 + p11 * q11;
		double det = j00 * j11 - j01 * j10;
		double tr = j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11;

		double det_2plus1 = 2.0 / (det * det) + 2.0;
		double det_3tr2 = 2.0*tr / (det * det * det);
		double r0 = area_now * (j00 * det_2plus1 - det_3tr2 * j11);
		double r1 = area_now * (j01 * det_2plus1 + det_3tr2 * j10);
		double r2 = area_now * (j10 * det_2plus1 + det_3tr2 * j01);
		double r3 = area_now * (j11 * det_2plus1 - det_3tr2 * j00);

		if (vid == f0)
		{
			p_g_(0) -= (r0*p00 + r1 * p01);
			p_g_(1) -= (r2*p00 + r3 * p01);
		}
		else if (vid == f1)
		{
			p_g_(0) -= (r0*p10 + r1 * p11);
			p_g_(1) -= (r2*p10 + r3 * p11);
		}
		else
		{
			p_g_(0) += (r0*(p00 + p10) + r1 * (p01 + p11));
			p_g_(1) += (r2*(p00 + p10) + r3 * (p01 + p11));
		}
	}
}

void SharedData::local_opt_1v(int vid)
{
	Eigen::Vector2d p_g_;
	local_decrease_1v(vid, p_g_);
	auto seed_h = mesh.vertex_handle(vid);
	//max-step
	double max_step_ = numeric_limits<double>::infinity();
	for (auto itvoh = mesh.voh_begin(seed_h); itvoh != mesh.voh_end(seed_h); itvoh++)
	{
		int f1 = mesh.to_vertex_handle(*itvoh).idx();
		int f2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(*itvoh)).idx();
		double q00 = w_uv(vid, 0) - w_uv(f2, 0);
		double q01 = w_uv(f1, 0) - w_uv(f2, 0);
		double q10 = w_uv(vid, 1) - w_uv(f2, 1);
		double q11 = w_uv(f1, 1) - w_uv(f2, 1);

		double b = q00 * q11 - q10 * q01;
		double a = p_g_(0)*q11 - p_g_(1)*q01;
		double t_tmp = -b / a;
		if (t_tmp > 0 && max_step_ > t_tmp)
		{
			max_step_ = t_tmp;
		}
	}
	double step_ = 0.9*max_step_;
	//line-search

	double h_ = 0.5;
	double c_ = 0.2;
	double tt = -(p_g_.squaredNorm());
	Eigen::Vector2d p_1v_(w_uv(vid, 0), w_uv(vid, 1));
	double ex = local_energy_1v(vid, p_1v_);
	//	printf("ex==: %f\n", ex);

	Eigen::Vector2d pos_new = p_1v_ + step_ * p_g_;
	double e = local_energy_1v(vid, pos_new);
	//	printf("e==: %f\n", e);

	double e_cache = e;
	double step_cache = step_;
	int count_c = 0;
	while (e > ex + step_ * c_*tt&&count_c < 5)
	{
		step_ *= h_;
		pos_new = p_1v_ + step_ * p_g_;
		e = local_energy_1v(vid, pos_new);
		count_c++;
		if (e < e_cache)
		{
			e_cache = e;
			step_cache = step_;
		}
		//	printf("e==: %f\n", e);
	}
	if (e_cache < e)
	{
		step_ = step_cache;
	}

	w_uv(vid, 0) = w_uv(vid, 0) + step_ * p_g_(0);
	w_uv(vid, 1) = w_uv(vid, 1) + step_ * p_g_(1);

}

double SharedData::local_energy_1v(int vid, const Eigen::Vector2d& p_)
{
	double e_sum = 0.;
	auto seed_h = mesh.vertex_handle(vid);
	for (auto itvf = mesh.vf_begin(seed_h); itvf != mesh.vf_end(seed_h); itvf++)
	{
		int fid_ = itvf->idx();
		int f0, f1, f2;
		double p00, p01, p10, p11, area_now;
		if (fid_ < mf_num)
		{
			f0 = F(fid_, 0);
			f1 = F(fid_, 1);
			f2 = F(fid_, 2);
			p00 = upd_inv00[fid_];
			p01 = upd_inv01[fid_];
			p10 = upd_inv10[fid_];
			p11 = upd_inv11[fid_];
			area_now = src_weight[fid_];
		}
		else
		{
			//continue;
			int fid_i = fid_ - mf_num;
			f0 = s_T(fid_i, 0);
			f1 = s_T(fid_i, 1);
			f2 = s_T(fid_i, 2);
			p00 = scaf_jac[fid_i][0];
			p01 = scaf_jac[fid_i][1];
			p10 = scaf_jac[fid_i][2];
			p11 = scaf_jac[fid_i][3];
			area_now = scaf_weight;
		}

		double q00, q01, q10, q11;
		if (vid == f0)
		{
			q00 = p_(0) - w_uv(f2, 0);
			q01 = w_uv(f1, 0) - w_uv(f2, 0);
			q10 = p_(1) - w_uv(f2, 1);
			q11 = w_uv(f1, 1) - w_uv(f2, 1);
		}
		else if (vid == f1)
		{
			q00 = w_uv(f0, 0) - w_uv(f2, 0);
			q01 = p_(0) - w_uv(f2, 0);
			q10 = w_uv(f0, 1) - w_uv(f2, 1);
			q11 = p_(1) - w_uv(f2, 1);
		}
		else
		{
			q00 = w_uv(f0, 0) - p_(0);
			q01 = w_uv(f1, 0) - p_(0);
			q10 = w_uv(f0, 1) - p_(1);
			q11 = w_uv(f1, 1) - p_(1);
		}

		double j00 = p00 * q00 + p10 * q01;
		double j01 = p01 * q00 + p11 * q01;
		double j10 = p00 * q10 + p10 * q11;
		double j11 = p01 * q10 + p11 * q11;
		double det = j00 * j11 - j01 * j10;
		double tr = j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11;
		double e_ = area_now * (tr + tr / (det*det));

		e_sum += e_;
	}
	return e_sum;
}

void SharedData::mesh_improve(bool in_packing = false)
{
	using namespace Eigen;
	if (rect_frame_V.size() == 0) {
		MatrixXd m_uv = w_uv.topRows(mv_num);
		Matrix2d ob;// = rect_corners;
		VectorXd uv_max = m_uv.colwise().maxCoeff();
		VectorXd uv_min = m_uv.colwise().minCoeff();
		VectorXd uv_mid = (uv_max + uv_min) / 2.;
		double maxlen = (uv_max - uv_mid).maxCoeff();
		uv_min(0) = uv_mid(0) - maxlen;
		uv_min(1) = uv_mid(1) - maxlen;
		uv_max(0) = uv_mid(0) + maxlen;
		uv_max(1) = uv_mid(1) + maxlen;

		//double scaf_range = 3;
		Eigen::Array2d scaf_range(3, 3);
		ob.row(0) = uv_mid.array() + scaf_range * ((uv_min - uv_mid).array());
		ob.row(1) = uv_mid.array() + scaf_range * ((uv_max - uv_mid).array());

		ob.row(1) = ob.row(1) - ob.row(0);
		for (int i = 0; i < mv_num; ++i)
		{
			w_uv.row(i) = w_uv.row(i) - ob.row(0);
		}
		ob.row(0).setZero();

		Vector2d rect_len;
		rect_len << ob(1, 0) - ob(0, 0), ob(1, 1) - ob(0, 1);
		int frame_points = 20;
		interval = rect_len(0);

		rect_frame_V.resize(4 * frame_points, 2);
		for (int i = 0; i < frame_points; i++)
		{
			// 0,0;0,1
			rect_frame_V.row(i) << ob(0, 0), ob(0, 1) + i * rect_len(1) / frame_points;
			// 0,0;1,1
			rect_frame_V.row(i + frame_points)
				<< ob(0, 0) + i * rect_len(0) / frame_points, ob(1, 1);
			// 1,0;1,1
			rect_frame_V.row(i + 2 * frame_points) << ob(1, 0), ob(1, 1)
				- i * rect_len(1) /
				frame_points;
			// 1,0;0,1
			rect_frame_V.row(i + 3 * frame_points)
				<< ob(1, 0) - i * rect_len(0) / frame_points, ob(0, 1);
			// 0,0;0,1
		}
		frame_ids = Eigen::VectorXi::LinSpaced(rect_frame_V.rows(), mv_num,
			mv_num + rect_frame_V.rows());
	}

	MatrixXd V_bnd;
	V_bnd.resize(internal_bnd.size(), 2);
	for (int i = 0; i < internal_bnd.size(); i++) // redoing step 1.
	{
		V_bnd.row(i) = w_uv.row(internal_bnd(i));
	}

	// Concatenate Vert and Edge
	MatrixXd V_;
	MatrixXi E;

	{
		V_.resize(V_bnd.rows() + rect_frame_V.rows(), V_bnd.cols());
		V_ << V_bnd, rect_frame_V;
	}
	E.resize(V_.rows(), 2);
	for (int i = 0; i < E.rows(); i++)
		E.row(i) << i, i + 1;
	int acc_bs = 0;
	for (auto bs : bnd_sizes) {
		E(acc_bs + bs - 1, 1) = acc_bs;
		acc_bs += bs;
	}
	E(V_.rows() - 1, 1) = acc_bs;
	assert(acc_bs == internal_bnd.size());


	//ycy H 存储的是每个chart的第一个三角形的重心
	MatrixXd H = MatrixXd::Zero(component_sizes.size(), 2);
	{
		int hole_f = 0;
		int hole_i = 0;
		for (auto cs : component_sizes) {
			for (int i = 0; i < 3; i++)
				H.row(hole_i) += w_uv.row(F(hole_f, i)); // redoing step 2
			hole_f += cs;
			hole_i++;
		}
	}
	H /= 3.;

	MatrixXd uv2;
	triangulate(V_, E, H, uv2, s_T);
	auto bnd_n = internal_bnd.size();

	for (auto i = 0; i < s_T.rows(); i++)
	{
		for (auto j = 0; j < s_T.cols(); j++)
		{
			auto &x = s_T(i, j);
			if (x < bnd_n) x = internal_bnd(x);
			else x += mv_num - bnd_n;
		}
	}


	{
		//w_T.resize(m_T.rows() + s_T.rows(), 3);
		//w_T << m_T, s_T;
	}

	w_uv.conservativeResize(mv_num - bnd_n + uv2.rows(), 2);
	w_uv.bottomRows(uv2.rows() - bnd_n) = uv2.bottomRows(-bnd_n + uv2.rows());
	update_scaffold();

	//writeObj(w_uv, surface_F, "G:/a.obj");
}


void SharedData::add_new_patch(const Eigen::MatrixXd &V_ref,
	const Eigen::MatrixXi &F_ref,
	const Eigen::RowVectorXd &center/*, const Eigen::VectorXd& M*/) {

	Eigen::MatrixXd uv_init;
	Eigen::VectorXi bnd;
	Eigen::MatrixXd bnd_uv;

	int num_holes = boundary_loops.size() - 1;

	bnd = Eigen::Map<Eigen::VectorXi>(boundary_loops[0].data(),
		boundary_loops[0].size());

	map_vertices_to_circle(V_ref, bnd, bnd_uv);
	bnd_uv *= sqrt(1.0 / M_PI);
	bnd_uv.rowwise() += center;
	mesh_measure += 1.0;
	std::cout << "Mesh Measure " << mesh_measure << "; number holes " << num_holes << std::endl;

	if (num_holes == 0)
	{
		if (bnd.rows() == V_ref.rows())
		{
			std::cout << "All vert on boundary" << std::endl;
			uv_init.resize(V_ref.rows(), 2);
			for (int i = 0; i < bnd.rows(); i++) {
				uv_init.row(bnd(i)) = bnd_uv.row(i);
			}
		}
		else
		{
			Tutte(V_ref.rows(), F_ref, bnd, bnd_uv, uv_init);
		}
	}
	else
	{
		printf("A single patch without holes is needed;\n");
		exit(0);
	}

	//writeObj(uv_init, F_ref, "tutte.obj");

	if (mv_num == 0)
	{
		w_uv = uv_init;
	}
	else
	{
		Eigen::MatrixXd m_uv = w_uv.topRows(mv_num);
		w_uv.resize(m_uv.rows() + uv_init.rows(), 2);
		w_uv << m_uv, uv_init;
	}

	//m_M.conservativeResize(mf_num + F_ref.rows());
	//m_M.bottomRows(F_ref.rows()) = M;

	for (auto cur_bnd : boundary_loops)
	{
		internal_bnd.conservativeResize(internal_bnd.size() + cur_bnd.size());
		internal_bnd.bottomRows(cur_bnd.size()) = Eigen::Map<Eigen::ArrayXi>(cur_bnd.data(), cur_bnd.size()) + mv_num;
		bnd_sizes.push_back(cur_bnd.size());
	}

	//m_T.conservativeResize(mf_num + F_ref.rows(), 3);
	//m_T.bottomRows(F_ref.rows()) = F_ref.array() + mv_num;
	mf_num += F_ref.rows();

	//m_V.conservativeResize(mv_num + V_ref.rows(), 3);
	//m_V.bottomRows(V_ref.rows()) = V_ref;
	mv_num += V_ref.rows();
}


void SharedData::find_components(Eigen::VectorXi& F_flag)
{
	int mf_n = mesh.n_faces();
	F_flag.setConstant(mf_n, -1);
	component_sizes.clear();
	int color_count = 0;
	for (auto ifh = mesh.faces_begin(); ifh != mesh.faces_end(); ifh++)
	{
		if (F_flag[ifh->idx()] != -1)
			continue;
		queue<int> Q;
		Q.push(ifh->idx());
		component_sizes.push_back(0);
		while (!Q.empty())
		{
			int seed_ = Q.front();
			F_flag[seed_] = color_count;
			component_sizes[color_count]++;
			Q.pop();
			auto f_h = mesh.face_handle(seed_);
			for (auto iff = mesh.cff_begin(f_h); iff != mesh.cff_end(f_h); iff++)
			{
				if (F_flag[iff->idx()] == -1)
				{
					Q.push(iff->idx());
				}
			}
		}
		color_count++;
	}
	printf("find_components // total %d compenents;\n", color_count);
}

void SharedData::init_V_F_VF()
{
	int mv_n = mesh.n_vertices();
	int mf_n = mesh.n_faces();
	V.resize(mv_n, 3);
	F.resize(mf_n, 3);
//	VF.resize(mv_n);
	//assemble V,F;
	for (auto iv = mesh.vertices_begin(); iv != mesh.vertices_end(); iv++)
	{
		auto p = mesh.point(*iv);
		V(iv->idx(), 0) = p[0];
		V(iv->idx(), 1) = p[1];
		V(iv->idx(), 2) = p[2];
		//vector<int> fs_vec;
		//for (auto itvf = mesh.vf_begin(*iv); itvf != mesh.vf_end(*iv); itvf++)
		//{
		//	fs_vec.push_back(itvf->idx());
		//}
		//VF[iv->idx()] = fs_vec;
	}
	int fvi = 0;
	for (auto fi = mesh.faces_begin(); fi != mesh.faces_end(); fi++)
	{
		fvi = 0;
		for (auto fv = mesh.fv_begin(*fi); fv != mesh.fv_end(*fi); fv++)
		{
			F(fi->idx(), fvi) = fv->idx();
			fvi++;
		}
	}
}

void SharedData::rescale_and_init_area()
{
	int mf_n = mesh.n_faces();
	src_area.resize(mf_n);
	//resize the src mesh to unit-surface-area;
	double area_sum = 0.0;
	for (auto f_h = mesh.faces_begin(); f_h != mesh.faces_end(); f_h++)
	{
		auto he_h = mesh.halfedge_handle(*f_h);
		double area_ = mesh.calc_sector_area(he_h);
		area_sum += area_;
		src_area[f_h->idx()] = area_;
	}

	for (double &area_ : src_area)
		area_ /= area_sum;
	double originmesh_area_sqrt = sqrt(area_sum);
	area_same_factor = 1.0 / originmesh_area_sqrt;
	printf("total_area: %f ;area_same_factor: %f ;\n", area_sum, area_same_factor);
	for (auto it1 = mesh.vertices_begin(); it1 != mesh.vertices_end(); it1++)
	{
		mesh.set_point(*it1, area_same_factor*mesh.point(*it1));
	}

	switch (dis_w_type)
	{
	case SharedData::AREA:
		src_weight = src_area;
		break;
	case SharedData::UNIFORM:
	{
		double uniform_a = 1.0 / mf_n;
		src_weight.resize(mf_n, uniform_a);
		break;
	}
	default:
		break;
	}

	distortion.resize(mf_n);


	//OpenMesh::IO::write_mesh(mesh, "co_w.obj");
}

void SharedData::prepare_data(int loops_num)
{
	//Eigen::VectorXi F_flag;
	//find_components(F_flag);

	component_sizes.push_back(mesh.n_faces());
	int component_num = component_sizes.size();

	rescale_and_init_area();
	init_V_F_VF();

	mv_num = mesh.n_vertices();
	mf_num = mesh.n_faces();
	printf("Vn: %d ; Fn: %d ;\n", mv_num, mf_num);

	init_loops(loops_num);

	if (component_num == 1)//single patch
	{
		mf_num = 0;
		mv_num = 0;
//		Eigen::VectorXd M;
//		M.resize(src_area.size());
//#pragma omp parallel for
//		for (int i = 0; i < src_area.size(); i++)
//		{
//			M[i] = src_area[i];
//		}
		add_new_patch(V, F, Eigen::RowVector2d(0., 0.)/*,M*/);

	}
	else// multiple patches
	{

//#pragma omp parallel for
//		for (int i = 0; i < src_area.size(); i++)
//		{
//			src_area[i]=m_M[i];
//		}
	}

	rect_frame_V.resize(0, 0);
	//mesh_improve(true);
	init_src_inv();

}

void SharedData::init_src_inv()
{
	using namespace Eigen;
	upd_inv00.resize(mf_num);
	upd_inv01.resize(mf_num);
	upd_inv10.resize(mf_num);
	upd_inv11.resize(mf_num);

#pragma omp parallel for
	for (int i = 0; i < mf_num; i++)
	{
		int f0 = F(i, 0);
		int f1 = F(i, 1);
		int f2 = F(i, 2);

		//if (src_weight[i] > min_tri_area)
		{
			Vector3d x_(V(f0, 0) - V(f2, 0), V(f0, 1) - V(f2, 1), V(f0, 2) - V(f2, 2));
			double x1_0 = x_.norm();
			x_ /= x1_0;
			Vector3d l_(V(f1, 0) - V(f2, 0), V(f1, 1) - V(f2, 1), V(f1, 2) - V(f2, 2));

			Vector3d n_ = x_.cross(l_);
			n_.normalize();
			Vector3d y_ = n_.cross(x_);
			double x2_0 = l_.dot(x_);
			double y2_0 = l_.dot(y_);

			upd_inv00[i] = 1.0 / x1_0;
			upd_inv01[i] = -x2_0 / (x1_0*y2_0);
			upd_inv10[i] = 0.0;
			upd_inv11[i] = 1.0 / y2_0;
		}

		if (src_weight[i] < min_tri_area)
		{
			cout << " area: " << src_weight[i] << endl;
		}
	}
	restore_upd2src();

}

void SharedData::restore_upd2src()
{
	is_interp = false;
}

void SharedData::init_loops(int loops_num)
{
//	find_boundary_loops();

	boundary_vs.clear();
	for (auto&var : boundary_loops)
		boundary_vs.insert(var.begin(), var.end());
	grow_from_seeds(order1_fs_b, boundary_vs, loops_num);
	boundary_loop_for_part(order1_fs_b, boundary_inner_loop);

}

void SharedData::find_boundary_loops()
{
	vector<int> heh_visit(mesh.n_halfedges(), 0);
	for (auto heh = mesh.halfedges_begin(); heh != mesh.halfedges_end(); heh++)
	{
		if (heh_visit[heh->idx()] == 1)
			continue;
		else
		{
			heh_visit[heh->idx()] = 1;
			if (mesh.is_boundary(*heh))
			{
				auto he_start = *heh;
				vector<int> bvs_;
				int from_v = mesh.from_vertex_handle(he_start).idx();
				int to_v_copy = mesh.to_vertex_handle(he_start).idx();
				do
				{
					auto right_hh = mesh.opposite_halfedge_handle(he_start);
					heh_visit[right_hh.idx()] = 1;
					from_v = mesh.from_vertex_handle(he_start).idx();
					boundary_es_del.push_back(right_hh.idx());
					bvs_.push_back(from_v);
					he_start = mesh.prev_halfedge_handle(he_start);
					heh_visit[he_start.idx()] = 1;
				} while (from_v != to_v_copy);
				boundary_loops.emplace_back(bvs_);
			}
		}
	}
}

void SharedData::grow_from_seeds(set<int>& fs_set, const set<int>& seeds_v, int loops_num)
{
	set<int> contain_vs;
	std::queue<pair<int, int>> Q;
	for (auto&vid : seeds_v)
	{
		Q.emplace(vid, 0);
		contain_vs.insert(vid);
		auto v_h = mesh.vertex_handle(vid);
		for (auto itvf = mesh.vf_begin(v_h); itvf != mesh.vf_end(v_h); itvf++)
		{
			fs_set.insert(itvf->idx());
		}
	}
	while (!Q.empty())
	{
		auto var = Q.front();
		Q.pop();
		auto v_h = mesh.vertex_handle(var.first);
		if (var.second < loops_num - 1)
			for (auto itvv = mesh.vv_begin(v_h); itvv != mesh.vv_end(v_h); itvv++)
			{
				if (contain_vs.count(itvv->idx()) == 0)
				{
					Q.emplace(itvv->idx(), var.second + 1);
					contain_vs.insert(itvv->idx());
					
					for (auto itvf = mesh.vf_begin(*itvv); itvf != mesh.vf_end(*itvv); itvf++)
					{
						int vf_ = itvf->idx();
						fs_set.insert(vf_);
					}

				}
			}
	}
	//printf("There %d faces in the region of %d loops-vertex around the boundary;\n", fs_set.size(), loops);

	cout << "===============check========================" << " " << fs_set.size() << endl;
	clock_t t0 = clock();
	vector<int> f_flag(mf_num, 0);
	for (auto&var : fs_set)
		f_flag[var] = 1;

	for (auto&var : order1_fs_b)
		f_flag[var] = 1;

	set<int> v_check_wait;
	for (auto&var : fs_set)
	{
		auto f_h = mesh.face_handle(var);
		for (auto itfv = mesh.fv_begin(f_h); itfv != mesh.fv_end(f_h); itfv++)
		{
			v_check_wait.insert(itfv->idx());
		}
	}

	queue<int> Qv;
	for (auto&var : v_check_wait)
	{
		Qv.push(var);
	}

	int h_c_ = 0;
	while (!Qv.empty())
	{
		int var = Qv.front();
		Qv.pop();
		auto v_h = mesh.vertex_handle(var);
		auto itvf = mesh.vf_begin(v_h);
		int score_ = 0;
		int flag_f_pre = f_flag[itvf->idx()];
		for (; itvf != mesh.vf_end(v_h); itvf++)
		{
			int flag_f_cur = f_flag[itvf->idx()];
			if (flag_f_cur != flag_f_pre)
				score_++;
			flag_f_pre = flag_f_cur;
		}

		if (score_ > 2)
		{
			h_c_++;
			for (auto itvoh = mesh.voh_begin(v_h); itvoh != mesh.voh_end(v_h); itvoh++)
			{
				int f_id_ = mesh.face_handle(*itvoh).idx();
				fs_set.insert(f_id_);
				f_flag[f_id_] = 1;
				int v_id_ = mesh.to_vertex_handle(*itvoh).idx();
				Qv.push(v_id_);
			}
		}
	}

	clock_t t1 = clock();
	cout << "===============check========================" << (t1 - t0) / 1000.0 << " " << fs_set.size() << " hc " << h_c_ << endl;

}

void SharedData::boundary_loop_for_part(set<int>& part, std::vector<int>& boundaryloop)
{
	typedef tuple<int, int, int> Triplet3;
	boundaryloop.clear();

	//set<Triplet3> boundary_es;
	//std::vector<Triplet3> boundaryEdges;


	std::vector<Triplet3> edges;
	set<int> b_es_clean;
	for (auto& var : part)
	{
		auto f_h = mesh.face_handle(var);
		for (auto ifh = mesh.fh_begin(f_h); ifh != mesh.fh_end(f_h); ifh++)
		{
			int v0 = mesh.from_vertex_handle(*ifh).idx();
			int v1 = mesh.to_vertex_handle(*ifh).idx();
			if (v0 > v1) std::swap(v0, v1);
			edges.emplace_back(v0, v1, ifh->idx());
		}
	}

	std::sort(edges.begin(), edges.end());
	int i = 1;
	for (; i < edges.size();)
	{
		auto& r1 = edges[i - 1];
		auto& r2 = edges[i];
		if ((get<0>(r1) == get<0>(r2)) && (get<1>(r1) == get<1>(r2)))
		{
			i += 2;
		}
		else
		{
			b_es_clean.insert(get<2>(edges[i - 1]));
			i++;
		}
	}
	if (i == edges.size())
		b_es_clean.insert(get<2>(edges.back()));

	for (auto& var : boundary_es_del)
		b_es_clean.erase(var);

	std::vector<std::vector<int>> boundaryloop_hh;
	while (!b_es_clean.empty())
	{
		int h_h_id = *(b_es_clean.begin());
		int next_h_h_id = h_h_id;
		vector<int> loop0;
		do
		{
			loop0.push_back(next_h_h_id);
			auto h_h = mesh.opposite_halfedge_handle(mesh.halfedge_handle(next_h_h_id));
			OpenMesh::HalfedgeHandle next_h_h = h_h;
			do
			{
				next_h_h = mesh.opposite_halfedge_handle(mesh.prev_halfedge_handle(next_h_h));
				if (b_es_clean.count(next_h_h.idx()) > 0)
					break;
			} while (next_h_h != h_h);
			next_h_h_id = next_h_h.idx();
			b_es_clean.erase(next_h_h_id);
		} while (next_h_h_id != h_h_id);
		boundaryloop_hh.emplace_back(loop0);
	}


	{
		set<int> part_add;
		int maxlong_id;
		int maxlong_nums = 0;
		for (size_t i = 0; i < boundaryloop_hh.size(); i++)
		{
			if (maxlong_nums < boundaryloop_hh[i].size())
			{
				maxlong_id = i;
				maxlong_nums = boundaryloop_hh[i].size();
			}
		}
		//vector<int> maxlong = boundaryloop_hh[maxlong_id];

		vector<int> visited_fs(mf_num, 0);
		for (size_t i = 0; i < boundaryloop_hh.size(); i++)
		{
			if (i == maxlong_id)
			{
				continue;
			}
			set<int> forbid_f_set;
			for (auto&var : boundaryloop_hh[i])
			{
				forbid_f_set.insert(mesh.face_handle(mesh.halfedge_handle(var)).idx());
			}
			queue<int> q;
			int seed_ = mesh.face_handle(mesh.opposite_halfedge_handle(mesh.halfedge_handle(boundaryloop_hh[i].front()))).idx();
			q.push(seed_);
			visited_fs[seed_] = 1;
			while (!q.empty())
			{
				auto fh_ = mesh.face_handle(q.front());
				part_add.insert(q.front());
				q.pop();
				for (auto itff = mesh.ff_begin(fh_); itff != mesh.ff_end(fh_); itff++)
				{
					if (visited_fs[itff->idx()] == 0 && forbid_f_set.count(itff->idx()) == 0)
					{
						q.push(itff->idx());
						visited_fs[itff->idx()] = 1;
					}
				}
			}
		}
		for (auto&var : part_add)
			part.insert(var);
		printf("There are totally %d faces set fixed at boundary region;\n", part.size());


		{
			first_inner_bl= boundaryloop_hh[maxlong_id];
			auto& var = boundaryloop_hh[maxlong_id];
			vector<int> loop0;
			int begin_v = mesh.from_vertex_handle(mesh.halfedge_handle(var.front())).idx();
			int end_v = mesh.to_vertex_handle(mesh.halfedge_handle(var.back())).idx();
			if (begin_v == end_v);
			cout << "First inner-boundary vertex number: " << var.size() << "; Success!" << endl;

			for (auto it = var.begin(); it != var.end(); it++)
			{
				loop0.push_back(mesh.from_vertex_handle(mesh.halfedge_handle(*it)).idx());
			}
			boundaryloop = loop0;
		}

	}

}

void SharedData::update_source_jacobian(const Eigen::VectorXd& pos_)
{

}

void SharedData::update_uv_from_outersrc(const Eigen::VectorXd & pos_)
{
	w_uv.block(0, 0, mv_num, 1)= pos_.block(0, 0, mv_num, 1);
	w_uv.block(0, 1, mv_num, 1)= pos_.block(mv_num, 0, mv_num, 1);
}

void SharedData::find_order1_region_near_boundary_face(set<int>& fs_set, int loops)
{
	fs_set.clear();
	std::vector<int> boundary_f;
	for (auto itf = mesh.faces_begin(); itf != mesh.faces_end(); itf++)
	{
		if (mesh.is_boundary(*itf))
			boundary_f.push_back(itf->idx());
	}
	vector<int> contain_fs(mf_num, 0);

	std::queue<pair<int, int>> Q;
	for (size_t i = 0; i < boundary_f.size(); i++)
	{
		Q.emplace(boundary_f[i], 0);
		contain_fs[boundary_f[i]] = 1;
		fs_set.insert(boundary_f[i]);
	}

	while (!Q.empty())
	{
		auto var = Q.front();
		Q.pop();
		auto f_h = mesh.face_handle(var.first);
		if (var.second < loops - 1)
			for (auto itff = mesh.ff_begin(f_h); itff != mesh.ff_end(f_h); itff++)
			{
				if (contain_fs[itff->idx()] == 0)
				{
					Q.emplace(itff->idx(), var.second + 1);
					contain_fs[itff->idx()] = 1;
					fs_set.insert(itff->idx());
				}
			}
	}
	printf("There %d faces in the region of %d loops-face around the boundary;\n", fs_set.size(), loops);
}

// find loops boundary faces and store them into fs_set
void SharedData::find_order1_region_near_boundary_vertex(set<int>& fs_set, int loops)
{
	fs_set.clear();
	std::set<int> boundary_v;
	for (auto&lp : boundary_loops)
		boundary_v.insert(lp.begin(), lp.end());
	grow_from_seeds(fs_set, boundary_v, loops);
}

void SharedData::find_order1_region_near_boundary_by_radius(set<int>& fs_set, const set<int>& seeds_v, int loops_num)
{
	for (auto&vid : seeds_v)
	{
		std::queue<pair<int, int>> Q;

		set<int> contain_vs, contain_fs;
		Q.emplace(vid, 0);
		contain_vs.insert(vid);
		auto v_h = mesh.vertex_handle(vid);

		double max_e = 0.;
		for (auto itvhe = mesh.vv_begin(v_h); itvhe !=mesh.vv_end(v_h); itvhe++)
		{		
			double e_l_ = (w_uv.row(itvhe->idx()) - w_uv.row(vid)).norm();
			if (max_e < e_l_)
				max_e = e_l_;
		}
		max_e *= loops_num;
		for (auto itvf = mesh.vf_begin(v_h); itvf != mesh.vf_end(v_h); itvf++)
		{
			contain_fs.insert(itvf->idx());
			fs_set.insert(itvf->idx());
		}

		while (!Q.empty())
		{
			auto var = Q.front();
			Q.pop();
			auto v_h = mesh.vertex_handle(var.first);
			double e_l_ = (w_uv.row(var.first) - w_uv.row(vid)).norm();
			if (/*var.second < 2*loops_num&&*/e_l_<max_e)
				for (auto itvv = mesh.vv_begin(v_h); itvv != mesh.vv_end(v_h); itvv++)
				{
					if (contain_vs.count(itvv->idx()) == 0)
					{
						Q.emplace(itvv->idx(), var.second + 1);
						contain_vs.insert(itvv->idx());

						for (auto itvf = mesh.vf_begin(*itvv); itvf != mesh.vf_end(*itvv); itvf++)
						{
							int vf_ = itvf->idx();
							if (contain_fs.count(vf_) == 0)
							{
								contain_fs.insert(vf_);
								fs_set.insert(vf_);
							}
						}
					}
				}
		}

	}

	//printf("There %d faces in the region of %d loops-vertex around the boundary;\n", fs_set.size(), loops);

	//check
	cout << "===============check========================" << " " << fs_set.size() << endl;
	clock_t t0 = clock();
	vector<int> f_flag(mf_num, 0);
	for (auto&var : fs_set)
		f_flag[var] = 1;

	for (auto&var : order1_fs_b)
		f_flag[var] = 1;

	set<int> v_check_wait;
	for (auto&var : fs_set)
	{
		auto f_h = mesh.face_handle(var);
		for (auto itfv = mesh.fv_begin(f_h); itfv !=mesh.fv_end(f_h) ; itfv++)
		{
			v_check_wait.insert(itfv->idx());
		}
	}

	queue<int> Qv;
	for (auto&var : v_check_wait)
	{
		Qv.push(var);
	}
 
	int h_c_ = 0;
	while (!Qv.empty())
	{
		int var = Qv.front();
		Qv.pop();
		auto v_h = mesh.vertex_handle(var);
		auto itvf = mesh.vf_begin(v_h);
		int score_ = 0;
		int flag_f_pre = f_flag[itvf->idx()];
		for (; itvf != mesh.vf_end(v_h); itvf++)
		{
			int flag_f_cur = f_flag[itvf->idx()];
			if (flag_f_cur != flag_f_pre)
				score_++;
			flag_f_pre = flag_f_cur;
		}

		if (score_ > 2)
		{
			h_c_++;
			for (auto itvoh = mesh.voh_begin(v_h); itvoh !=mesh.voh_end(v_h) ; itvoh++)
			{
				int f_id_ = mesh.face_handle(*itvoh).idx();
				fs_set.insert(f_id_);
				f_flag[f_id_] = 1;
				int v_id_=mesh.to_vertex_handle(*itvoh).idx();
				Qv.push(v_id_);
			}
		}
	}

	clock_t t1 = clock();
	cout << "===============check========================"<<(t1-t0)/1000.0 << " " << fs_set.size() << " hc " << h_c_<< endl;

}

