#include "Parafun.h"
#define eps 1e-16
using namespace Eigen;
using namespace std;
Parafun::~Parafun()
{
	delete spline_pardiso;
	spline_pardiso = NULL;
}

Parafun::Parafun(SharedData & data) :d_(data)
{
	iter_general_cout = 0;
	spline_pardiso = NULL;
	spline_surface = NULL;
};

void Parafun::init(int grid_num_)
{
	if (spline_pardiso != NULL)
	{
		delete spline_pardiso;
		spline_pardiso = NULL;
	}
	if (spline_surface != NULL)
	{
		delete spline_surface;
		spline_surface = NULL;
	}
	spline_surface = new BSplineSurface();
	spline_surface->SetDegree(3, 3);
	spline_surface->SetNumberControlPoints(grid_num_, grid_num_);
	//	spline_surface->SetIdentity();

	int m = spline_surface->ControlPoints().size();
	int n = spline_surface->ControlPoints()[0].size();
	total_cntrl = m * n;
	Control_N = (m - 2)*(n - 2);
	Interval_N = (m - 3)*(n - 3);

	Assist_N = 0;

	cntrl_vs_index();
	Pre_calculate_spline();

	int mf_n = d_.mf_num;

	u_x_m.resize(mf_n);
	u_y_m.resize(mf_n);
	v_x_m.resize(mf_n);
	v_y_m.resize(mf_n);
	//	energy_per_face.resize(mf_n);
	pos_of_control.resize(total_cntrl * 2);
}

void Parafun::cntrl_vs_index()
{
	const auto & ctrlpts = spline_surface->ControlPoints();
	int m = ctrlpts.size();
	int n = ctrlpts[0].size();
	int bnd_n = 2 * (m + n) - 4;
	//	fix_cntrls.resize(bnd_n, d_.dim);
	var_cntrls.clear();
	cntrl2index.clear();
	var_cntrls.reserve(total_cntrl);
	cntrl2index.resize(total_cntrl, -1);

	for (size_t i = 1; i < m - 1; i++)
	{
		for (size_t j = 1; j < n - 1; j++)
		{
			var_cntrls.push_back(i*n + j);
		}
	}

	for (size_t i = 0; i < var_cntrls.size(); i++)
	{
		cntrl2index[var_cntrls[i]] = i;
	}
}
void Parafun::Pre_calculate_spline()
{
	spline_pardiso = new PardisoSolver();

	typedef Triplet<int> T;
	std::vector<T> tripletlist;

	spline_pardiso->ia.clear(); spline_pardiso->ia.reserve(2 * Control_N + 1);
	spline_pardiso->ja.clear(); spline_pardiso->ja.reserve(49 * Control_N);

	int m = spline_surface->ControlPoints().size();
	int n = spline_surface->ControlPoints()[0].size();

	int row_con;
	int  index = 0;
	for (int i = 1; i < m - 1; ++i)
	{
		for (int j = 1; j < n - 1; ++j)
		{
			int dd = 0;
			row_con = (i - 1)*(n - 2) + j - 1;
			spline_pardiso->ia.push_back(spline_pardiso->ja.size());
			for (int j1 = j; j1 < min(n - 1, j + 4); ++j1)
			{
				index = (i - 1)*(n - 2) + j1 - 1;
				spline_pardiso->ja.push_back(index);
				tripletlist.push_back(T(row_con, index, dd));
				dd++;
			}
			for (int i1 = i + 1; i1 < min(m - 1, i + 4); ++i1)
			{
				for (int j1 = max(1, j - 3); j1 < min(n - 1, j + 4); ++j1)
				{
					index = (i1 - 1)*(n - 2) + j1 - 1;
					spline_pardiso->ja.push_back(index);
					tripletlist.push_back(T(row_con, index, dd));
					dd++;
				}
			}
			for (int i1 = max(1, i - 3); i1 < min(m - 1, i + 4); ++i1)
			{
				for (int j1 = max(1, j - 3); j1 < min(n - 1, j + 4); ++j1)
				{
					index = (i1 - 1)*(n - 2) + j1 - 1;
					spline_pardiso->ja.push_back(index + Control_N);
					tripletlist.push_back(T(row_con, index + Control_N, dd));
					dd++;
				}
			}
		}
	}

	for (int i = 1; i < m - 1; ++i)
	{
		for (int j = 1; j < n - 1; ++j)
		{
			int dd = 0;
			row_con = (i - 1)*(n - 2) + j - 1 + Control_N;
			spline_pardiso->ia.push_back(spline_pardiso->ja.size());
			for (int j1 = j; j1 < min(n - 1, j + 4); ++j1)
			{
				index = (i - 1)*(n - 2) + j1 - 1 + Control_N;
				spline_pardiso->ja.push_back(index);
				dd++;
			}
			for (int i1 = i + 1; i1 < min(m - 1, i + 4); ++i1)
			{
				for (int j1 = max(1, j - 3); j1 < min(n - 1, j + 4); ++j1)
				{
					index = (i1 - 1)*(n - 2) + j1 - 1 + Control_N;
					spline_pardiso->ja.push_back(index);
					dd++;
				}
			}
		}
	}

	find_cntrl_in_rows.resize(Control_N, 2 * Control_N);
	find_cntrl_in_rows.setFromTriplets(tripletlist.begin(), tripletlist.end());

	spline_pardiso->ia.push_back(spline_pardiso->ja.size());


	spline_pardiso->nnz = spline_pardiso->ja.size();
	spline_pardiso->a.resize(spline_pardiso->nnz);
	spline_pardiso->num = 2 * Control_N;
	spline_pardiso->pardiso_init();
}

void Parafun::after_mesh_improve()
{
	d_.adjust_scaf_weight(d_.Energy_cur*d_.mesh_measure / d_.sf_num / 100.0);
	d_.update_scaf_jacobian();

	total_num = d_.V_num;
	V_N = d_.V_num - d_.frame_ids.size();
	F_N = d_.F_num;
}

void Parafun::spline_init()
{
	spline_surface->UpdateKnots(d_.interval);
	spline_surface->SetIdentity();

	int m = spline_surface->ControlPoints().size();
	int n = spline_surface->ControlPoints()[0].size();
	int degx = static_cast<int>(spline_surface->DegreeX());
	int degy = static_cast<int>(spline_surface->DegreeY());
	double inter_len_ = spline_surface->get_inter_len();
	int mf_n = d_.mf_num;
	int sf_n = d_.sf_num;
	int f_n = mf_n + sf_n;

	//sample

	samplepointsx.clear();
	samplepointsy.clear();
	sample_at_face.clear();
	samplepointsx.reserve(mf_n * Sample_per_face + sf_n * Sample_per_scaffold);
	samplepointsy.reserve(mf_n * Sample_per_face + sf_n * Sample_per_scaffold);
	sample_at_face.reserve(mf_n * Sample_per_face + sf_n);

	//	const double *pos = position_of_mesh.data();
	double cx, cy;
	double x0, y0, x1, y1, x2, y2;
	int id = 0;
	int id2 = 0;

	for (int i = 0; i < d_.mf_num; ++i)
	{
		int f0 = d_.F(i, 0);
		int f1 = d_.F(i, 1);
		int f2 = d_.F(i, 2);

		x0 = d_.w_uv(f0, 0);
		y0 = d_.w_uv(f0, 1);
		x1 = d_.w_uv(f1, 0);
		y1 = d_.w_uv(f1, 1);
		x2 = d_.w_uv(f2, 0);
		y2 = d_.w_uv(f2, 1);

		cx = (x0 + x1 + x2) / 3.0;
		cy = (y0 + y1 + y2) / 3.0;

		samplepointsx.push_back(cx);
		samplepointsy.push_back(cy);
		sample_at_face.push_back(i);
		if (d_.distortion[i] > 0.01*max_energy)
		{
			id++;
			samplepointsx.push_back((cx + x0) / 2.0);
			samplepointsy.push_back((cy + y0) / 2.0);
			sample_at_face.push_back(i);
			samplepointsx.push_back((cx + x1) / 2.0);
			samplepointsy.push_back((cy + y1) / 2.0);
			sample_at_face.push_back(i);
			samplepointsx.push_back((cx + x2) / 2.0);
			samplepointsy.push_back((cy + y2) / 2.0);
			sample_at_face.push_back(i);

		}
	}
	Sample_N_m = sample_at_face.size();
//	std::cout << max_energy << "	" << id << "	" << id2 << "	" << Sample_N_m << std::endl;

	for (int i = 0; i < d_.sf_num; ++i)
	{
		int f0 = d_.s_T(i, 0);
		int f1 = d_.s_T(i, 1);
		int f2 = d_.s_T(i, 2);

		x0 = d_.w_uv(f0, 0);
		y0 = d_.w_uv(f0, 1);
		x1 = d_.w_uv(f1, 0);
		y1 = d_.w_uv(f1, 1);
		x2 = d_.w_uv(f2, 0);
		y2 = d_.w_uv(f2, 1);

		cx = (x0 + x1 + x2) / 3.0;
		cy = (y0 + y1 + y2) / 3.0;

		samplepointsx.push_back(cx);
		samplepointsy.push_back(cy);
		sample_at_face.push_back(i);
	}

	Sample_N = samplepointsx.size();

	vector<int> is_assist;
	is_assist.resize(Interval_N, -1);
	for (int i = 0; i < Sample_N; i++)
	{
		int rx = (int)(samplepointsx[i] / inter_len_) + degx;
		int ry = (int)(samplepointsy[i] / inter_len_) + degy;
		int index = (rx - degx) * (n - 3) + (ry - degy);
		is_assist[index] = 1;
	}

	for (size_t i = 0; i < m; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			pos_of_control[i*n + j] = spline_surface->ControlPoints()[i][j][0];
			pos_of_control[i*n + j + total_cntrl] = spline_surface->ControlPoints()[i][j][1];
		}
	}

	Assist_N = 0;
	for (int k = 0; k < Interval_N; ++k)
	{
		if (is_assist[k] < 0)
		{
			int i = k / (n - 3);
			int j = k % (n - 3);
			double x = (i + 0.5) * inter_len_;
			double y = (j + 0.5) * inter_len_;
			samplepointsx.push_back(x);
			samplepointsy.push_back(y);
			Assist_N++;
		}
	}

	samplepointsx.shrink_to_fit();
	samplepointsy.shrink_to_fit();
	sample_at_face.shrink_to_fit();
}

double Parafun::perform_1iteration_spline(bool use_CM)
{
	long time_beg, time_end;
	time_beg = clock();
	if (use_CM)
		CM_spline();
	else
		SLIM_spline();


	time_end = clock();
	time_consumption = (time_end - time_beg) / 1000.0;

	return compute_energy(d_.w_uv, false);

}
void Parafun::CM_spline()
{
	spline_init();
	typedef Eigen::Triplet<double> Tri;

	int m = spline_surface->ControlPoints().size();
	int n = spline_surface->ControlPoints()[0].size();
	int degx = static_cast<int>(spline_surface->DegreeX());
	int degy = static_cast<int>(spline_surface->DegreeY());
	int l_upper_ = (degx + 1)*(degy + 1);

//	double *pos = pos_of_control.data();
	spline_pardiso->a.clear();
	spline_pardiso->a.resize(spline_pardiso->nnz, 0.);
	spline_pardiso->rhs.clear();
	spline_pardiso->rhs.resize(2 * Control_N, 0.0);

	clock_t t11 = clock();

#pragma omp parallel for
	for (int i = 0; i < Sample_N_m; i++)
	{
		vector<pair<int, double>> hess_compents_record;
		vector<pair<int, double>> grad_compents_record;

		int var_num_per_sample = 2 * (degx + 1)*(degy + 1);
		hess_compents_record.reserve(var_num_per_sample*(var_num_per_sample + 1) / 2);
		grad_compents_record.reserve(var_num_per_sample);

		int f_id = sample_at_face[i];
		double area_now = d_.src_weight[f_id];
		double j00 = u_x_m[f_id];
		double j01 = u_y_m[f_id];
		double j10 = v_x_m[f_id];
		double j11 = v_y_m[f_id];

		double tr = j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11;
		double det = j00 * j11 - j01 * j10;

		double alpha_0 = j00 + j11;
		double alpha_1 = j10 - j01;
		double beta_0 = j00 - j11;
		double beta_1 = j10 + j01;

		double alpha_norm = sqrt(alpha_0*alpha_0 + alpha_1 * alpha_1);
		double beta_norm = sqrt(beta_0*beta_0 + beta_1 * beta_1);

		double sig0 = 0.5*(alpha_norm + beta_norm);
		double sig1 = 0.5*(alpha_norm - beta_norm);

		std::vector<double> coeff_u, coeff_v;
		int grid_x, grid_y;
		coeff_u.resize(l_upper_);
		coeff_v.resize(l_upper_);
		std::vector<double> du, dv;
		spline_surface->Derivative(samplepointsx[i], samplepointsy[i], grid_x, grid_y, du, dv);

		for (int l_idx_ = 0; l_idx_ < l_upper_; l_idx_++)
		{
			coeff_u[l_idx_] = du[l_idx_] * u_x_m[f_id] + dv[l_idx_] * v_x_m[f_id];
			coeff_v[l_idx_] = du[l_idx_] * u_y_m[f_id] + dv[l_idx_] * v_y_m[f_id];
		}

		double det_1 = 1.0 / det;
		double det_2plus1 = det_1 * det_1;
		double tr_det_4_3 = 3.*tr*det_2plus1 * det_2plus1;
		double det_3_2 = 2.*det_2plus1 * det_1;

		double det_3tr2 = det_2plus1 * tr*det_1;
		det_2plus1 += 1.0;
		double r0 = area_now * (j00 * det_2plus1 - det_3tr2 * j11);
		double r1 = area_now * (j01 * det_2plus1 + det_3tr2 * j10);
		double r2 = area_now * (j10 * det_2plus1 + det_3tr2 * j01);
		double r3 = area_now * (j11 * det_2plus1 - det_3tr2 * j00);


		double uu0 = (det_2plus1 + tr_det_4_3 * j11*j11 - 2.*det_3_2*j00*j11)*area_now;
		double uu12 = (det_3_2*(j00*j10 - j01 * j11) - tr_det_4_3 * j11*j10)*area_now;
		double uu3 = (det_2plus1 + tr_det_4_3 * j10*j10 + 2.*det_3_2*j01*j10)*area_now;

		double vv0 = (det_2plus1 + tr_det_4_3 * j01*j01 + 2.*det_3_2*j01*j10)*area_now;
		double vv12 = (det_3_2*(j01 * j11 - j00 * j10) - tr_det_4_3 * j00*j01)*area_now;
		double vv3 = (det_2plus1 + tr_det_4_3 * j00*j00 - 2.*det_3_2*j00*j11)*area_now;

		double uv0 = (det_3_2*(j00*j01 - j10 * j11) - tr_det_4_3 * j01*j11)*area_now;
		double uv1 = (tr_det_4_3*j00*j11 - det_3_2 * (j00*j00 + j11 * j11 + tr / 2.))*area_now;
		double uv2 = (tr_det_4_3*j01*j10 + det_3_2 * (j01*j01 + j10 * j10 + tr / 2.))*area_now;
		double uv3 = (det_3_2*(j10 * j11 - j00 * j01) - tr_det_4_3 * j00*j10)*area_now;

		double sig_flag = (sig1 + sig0 - 1.0 / (sig1 * sig1 * sig1) - 1.0 / (sig0 * sig0 * sig0));

		double sig_u0, sig_u1, sig_u2;
		bool is_sig_flag = (sig_flag < 0);
		if (is_sig_flag)
		{
			double alpha_1 = 1.0 / alpha_norm;
			double alpha_3 = 0.5*alpha_1 * alpha_1*alpha_1;
			sig_u0 = (j01 - j10)*(j01 - j10)*alpha_3;
			sig_u1 = (j00 + j11)*(j01 - j10)*alpha_3;
			sig_u2 = (j00 + j11)*(j00 + j11)*alpha_3;
			sig_flag *= area_now;
			sig_u0 *= sig_flag;
			sig_u1 *= sig_flag;
			sig_u2 *= sig_flag;
			//sig_u0//-sig_u1//..//sig_u2
			//sig_u2//sig_u1//..//sig_u0
			//sig_u1//sig_u0//-sig_u2//-sig_u1
		}
		for (int i1 = 0; i1 < degx + 1; i1++)
		{
			for (int j1 = 0; j1 < degy + 1; j1++)
			{
				int id = (i1 + grid_x - degx) * n + (j1 + grid_y - degy);
				int var = cntrl2index[id];
				if (var == -1)
					continue;

				int l_var_n1 = i1 * (degy + 1) + j1;
				double g_u_ = r0 * coeff_u[l_var_n1] + r1 * coeff_v[l_var_n1];
				double g_v_ = r2 * coeff_u[l_var_n1] + r3 * coeff_v[l_var_n1];
				grad_compents_record.emplace_back(var, g_u_);
				grad_compents_record.emplace_back(var + Control_N, g_v_);

				for (int i2 = 0; i2 < degx + 1; i2++)
				{
					for (int j2 = 0; j2 < degy + 1; j2++)
					{
						int id2 = (i2 + grid_x - degx) * n + (j2 + grid_y - degy);
						int var2 = cntrl2index[id2];
						if (var2 == -1)
							continue;
						int l_var_n2 = i2 * (degy + 1) + j2;

						double uij11 = coeff_u[l_var_n1] * coeff_u[l_var_n2];
						double uij12 = coeff_u[l_var_n1] * coeff_v[l_var_n2];
						double uij21 = coeff_v[l_var_n1] * coeff_u[l_var_n2];
						double uij22 = coeff_v[l_var_n1] * coeff_v[l_var_n2];

						double h_ui_vj_plus = uv0 * uij11 + uv1 * uij12 + uv2 * uij21 + uv3 * uij22;
						if (is_sig_flag)
						{
							h_ui_vj_plus -= (sig_u1 * uij11 + sig_u0 * uij12 - sig_u2 * uij21 - sig_u1 * uij22);
						}
						hess_compents_record.emplace_back(spline_pardiso->ia[var] + find_cntrl_in_rows.coeff(var, var2 + Control_N) - 1, h_ui_vj_plus);
						if (var2 >= var)
						{
							double h_ui_uj_plus = uu0 * uij11 + uu12 * uij12 + uu12 * uij21 + uu3 * uij22;
							double h_vi_vj_plus = vv0 * uij11 + vv12 * uij12 + vv12 * uij21 + vv3 * uij22;
							if (is_sig_flag)
							{
								h_ui_uj_plus -= (sig_u0 * uij11 - sig_u1 * uij12 - sig_u1 * uij21 + sig_u2 * uij22);
								h_vi_vj_plus -= (sig_u2 * uij11 + sig_u1 * uij12 + sig_u1 * uij21 + sig_u0 * uij22);
							}
							int d = find_cntrl_in_rows.coeff(var, var2);
							hess_compents_record.emplace_back(spline_pardiso->ia[var] + d - 1, h_ui_uj_plus);
							hess_compents_record.emplace_back(spline_pardiso->ia[var + Control_N] + d - 1, h_vi_vj_plus);
						}
					}
				}
			}
		}
#pragma omp critical
		{
			for (const auto&var : hess_compents_record)
				spline_pardiso->a[var.first] += var.second;
			for (const auto&var : grad_compents_record)
				spline_pardiso->rhs[var.first] -= var.second;
		}
	}

#pragma omp parallel for
	for (int i = Sample_N_m; i < Sample_N + Assist_N; i++)
	{
		vector<pair<int, double>> hess_compents_record;
		int var_num_per_sample = 2 * (degx + 1)*(degy + 1);
		hess_compents_record.reserve(var_num_per_sample*(var_num_per_sample + 1) / 2);

		double area_now = d_.scaf_weight;
		double j00 = 1.;
		double j01 = 0.;

		if (i < Sample_N)
		{
			int f0 = d_.s_T(sample_at_face[i], 0);
			int f1 = d_.s_T(sample_at_face[i], 1);
			int f2 = d_.s_T(sample_at_face[i], 2);
			double q00 = d_.w_uv(f0, 0) - d_.w_uv(f2, 0);
			double q01 = d_.w_uv(f1, 0) - d_.w_uv(f2, 0);
			j00 = d_.scaf_jac[sample_at_face[i]][0] * q00;
			j01 = d_.scaf_jac[sample_at_face[i]][1] * q00 + d_.scaf_jac[sample_at_face[i]][3] * q01;
		}

		area_now *= 2.;

		double a2p1 = (1. + j00 * j00)*area_now;
		double ab = j00 * j01*area_now;
		double b2p1 = (1. + j01 * j01)*area_now;
		double uv1 = (j00*j00 - 1.)*area_now;
		double uv2 = (1. - j01 * j01)*area_now;

		std::vector<double> coeff_u, coeff_v;
		int grid_x, grid_y;
		spline_surface->Derivative(samplepointsx[i], samplepointsy[i], grid_x, grid_y, coeff_u, coeff_v);

		for (int i1 = 0; i1 < degx + 1; i1++)
		{
			for (int j1 = 0; j1 < degy + 1; j1++)
			{
				int id = (i1 + grid_x - degx) * n + (j1 + grid_y - degy);
				int var = cntrl2index[id];
				if (var == -1)
					continue;
				int l_var_n1 = i1 * (degy + 1) + j1;
				for (int i2 = 0; i2 < degx + 1; i2++)
				{
					for (int j2 = 0; j2 < degy + 1; j2++)
					{
						int id2 = (i2 + grid_x - degx) * n + (j2 + grid_y - degy);
						int var2 = cntrl2index[id2];
						if (var2 == -1)
							continue;
						int l_var_n2 = i2 * (degy + 1) + j2;

						double uij11 = coeff_u[l_var_n1] * coeff_u[l_var_n2];
						double uij12 = coeff_u[l_var_n1] * coeff_v[l_var_n2];
						double uij21 = coeff_v[l_var_n1] * coeff_u[l_var_n2];
						double uij22 = coeff_v[l_var_n1] * coeff_v[l_var_n2];

						double h_ui_vj_plus = -ab * uij11 + uv1 * uij12 + uv2 * uij21 + ab * uij22;

						hess_compents_record.emplace_back(spline_pardiso->ia[var] + find_cntrl_in_rows.coeff(var, var2 + Control_N) - 1, h_ui_vj_plus);
						if (var2 >= var)
						{
							double h_ui_uj_plus = a2p1 * uij11 + ab * uij12 + ab * uij21 + b2p1 * uij22;
							double h_vi_vj_plus = b2p1 * uij11 - ab * uij12 - ab * uij21 + a2p1 * uij22;

							int d = find_cntrl_in_rows.coeff(var, var2);
							hess_compents_record.emplace_back(spline_pardiso->ia[var] + d - 1, h_ui_uj_plus);
							hess_compents_record.emplace_back(spline_pardiso->ia[var + Control_N] + d - 1, h_vi_vj_plus);
						}
					}
				}
			}
		}

#pragma omp critical
		{
			for (const auto&var : hess_compents_record)
				spline_pardiso->a[var.first] += var.second;
		}

	}

	clock_t t22 = clock();
	printf("SamplePoints: %d ; spline-assmbly: %f ;", Sample_N + Assist_N, (t22 - t11) / 1000.0);
	LOG(INFO) << "SamplePoints: " << Sample_N + Assist_N << " ; spline-assmbly: " << (t22 - t11) / 1000.0;

	clear_samples();
	spline_pardiso->factorize();
	spline_pardiso->pardiso_solver();

	VectorXd d(2 * total_cntrl);
	d.setZero();
#pragma omp parallel for
	for (int i = 0; i < Control_N; i++)
	{
		int id = var_cntrls[i];
		d(id) = spline_pardiso->result[i];
		d(id + total_cntrl) = spline_pardiso->result[i + Control_N];
	}

	double temp_t;
	double alpha;
	pos_of_control += d;

	spline_pardiso->free_numerical_factorization_memory();

	Eigen::MatrixXd direction;
	Reconstruct(direction);
	Eigen::MatrixXd negative_grad_norm;
	calc_gradient_norm(negative_grad_norm);
	double neg_grad_dot_d = -(negative_grad_norm.col(0).dot(direction.col(0)));
	neg_grad_dot_d -= (negative_grad_norm.col(1).dot(direction.col(1)));
	if (neg_grad_dot_d < 0)
	{
		double temp_t1;
		max_step(d_.w_uv, direction, temp_t1);
		alpha = 0.95 * temp_t1;
		backtracking_line_search(d_.w_uv, direction, neg_grad_dot_d, alpha);
		d_.w_uv += alpha * direction;
	}
	iter_general_cout++;
}
void Parafun::SLIM_spline()
{
	spline_init();
	typedef Eigen::Triplet<double> Tri;

	int m = spline_surface->ControlPoints().size();
	int n = spline_surface->ControlPoints()[0].size();
	int degx = static_cast<int>(spline_surface->DegreeX());
	int degy = static_cast<int>(spline_surface->DegreeY());
	int l_upper_ = (degx + 1)*(degy + 1);


	spline_pardiso->a.clear();
	spline_pardiso->a.resize(spline_pardiso->ja.size(), 0.);
	spline_pardiso->rhs.clear();
	spline_pardiso->rhs.resize(2 * Control_N, 0.0);

	clock_t t11 = clock();

#pragma omp parallel for
	for (int i = 0; i < Sample_N_m; i++)
	{
		vector<pair<int, double>> hess_compents_record;
		vector<pair<int, double>> grad_compents_record;

		int var_num_per_sample = 2 * (degx + 1)*(degy + 1);
		hess_compents_record.reserve(var_num_per_sample*(var_num_per_sample + 1) / 2);
		grad_compents_record.reserve(var_num_per_sample);

		int f_id = sample_at_face[i];
		double area_now = d_.src_weight[f_id];
		double j00 = u_x_m[f_id];
		double j01 = u_y_m[f_id];
		double j10 = v_x_m[f_id];
		double j11 = v_y_m[f_id];

		double alpha_0 = j00 + j11;
		double alpha_1 = j10 - j01;
		double beta_0 = j00 - j11;
		double beta_1 = j10 + j01;

		double alpha_norm = 0.5*sqrt(alpha_0*alpha_0 + alpha_1 * alpha_1);
		double beta_norm = 0.5*sqrt(beta_0*beta_0 + beta_1 * beta_1);

		double sig0 = alpha_norm + beta_norm;
		double sig1 = alpha_norm - beta_norm;

		double sig0_1 = 1.0 / sig0;
		double sig1_1 = 1.0 / sig1;

		double new_sig0 = sqrt(1.0 + sig0_1 + sig0_1 * sig0_1 + sig0_1 * sig0_1*sig0_1);
		double new_sig1 = sqrt(1.0 + sig1_1 + sig1_1 * sig1_1 + sig1_1 * sig1_1*sig1_1);

		double temp = 0.;
		if (beta_norm > 1e-10)
		{
			temp = (new_sig1 - new_sig0) / (sig1*sig1 - sig0 * sig0);
		}
		double sig_square_sum = 0.5*(sig0 * sig0 + sig1 * sig1);
		double sig_new_sum = 0.5*(new_sig0 + new_sig1);

		double w00 = temp * (j00*j00 + j01 * j01 - sig_square_sum) + sig_new_sum;
		double w01 = temp * (j00*j10 + j01 * j11);
		double w11 = temp * (j10*j10 + j11 * j11 - sig_square_sum) + sig_new_sum;

		double w1 = area_now * (w00 * w00 + w01 * w01);
		double w2 = area_now * (w00 + w11)* w01;
		double w3 = area_now * (w01 * w01 + w11 * w11);

		double det = j00 * j11 - j01 * j10;
		double tr = (j00*j00 + j01 * j01 + j10 * j10 + j11 * j11);
		double det_1 = 1.0 / det;
		double det_2plus1 = det_1 * det_1;
		double det_3tr2 = det_2plus1 * tr*det_1;
		det_2plus1 += 1.0;
		double r0 = area_now * (j00 * det_2plus1 - det_3tr2 * j11);
		double r1 = area_now * (j01 * det_2plus1 + det_3tr2 * j10);
		double r2 = area_now * (j10 * det_2plus1 + det_3tr2 * j01);
		double r3 = area_now * (j11 * det_2plus1 - det_3tr2 * j00);

		std::vector<double> coeff_u, coeff_v;
		int grid_x, grid_y;
		coeff_u.resize(l_upper_);
		coeff_v.resize(l_upper_);
		std::vector<double> du, dv;
		spline_surface->Derivative(samplepointsx[i], samplepointsy[i], grid_x, grid_y, du, dv);

		for (int l_idx_ = 0; l_idx_ < l_upper_; l_idx_++)
		{
			coeff_u[l_idx_] = du[l_idx_] * u_x_m[f_id] + dv[l_idx_] * v_x_m[f_id];
			coeff_v[l_idx_] = du[l_idx_] * u_y_m[f_id] + dv[l_idx_] * v_y_m[f_id];
		}

		for (int i1 = 0; i1 < degx + 1; i1++)
		{
			for (int j1 = 0; j1 < degy + 1; j1++)
			{
				int id = (i1 + grid_x - degx) * n + (j1 + grid_y - degy);
				int var = cntrl2index[id];
				if (var == -1)
					continue;
				int l_var_n1 = i1 * (degy + 1) + j1;
				grad_compents_record.emplace_back(var, r0* coeff_u[l_var_n1] + r1 * coeff_v[l_var_n1]);
				grad_compents_record.emplace_back(var + Control_N, r2* coeff_u[l_var_n1] + r3 * coeff_v[l_var_n1]);

				for (int i2 = 0; i2 < degx + 1; i2++)
				{
					for (int j2 = 0; j2 < degy + 1; j2++)
					{
						int id2 = (i2 + grid_x - degx) * n + (j2 + grid_y - degy);
						int var2 = cntrl2index[id2];
						if (var2 == -1)
							continue;
						int l_var_n2 = i2 * (degy + 1) + j2;
						double dx2dy2 = coeff_u[l_var_n1] * coeff_u[l_var_n2] + coeff_v[l_var_n1] * coeff_v[l_var_n2];
						hess_compents_record.emplace_back(spline_pardiso->ia[var] + find_cntrl_in_rows.coeff(var, var2 + Control_N) - 1, w2 * dx2dy2);
						//spline_pardiso_a[spline_pardiso->ia[var] + find_cntrl_in_rows.coeff(var, var2 + Control_N)] += w2 * (dxdx + dydy);
						if (var2 >= var)
						{
							int d = find_cntrl_in_rows.coeff(var, var2);

							hess_compents_record.emplace_back(spline_pardiso->ia[var] + d - 1, w1 * dx2dy2);
							hess_compents_record.emplace_back(spline_pardiso->ia[var + Control_N] + d - 1, w3 * dx2dy2);
							//spline_pardiso_a[spline_pardiso->ia[var] + d] += w1 * (dxdx + dydy);
							//spline_pardiso_a[spline_pardiso->ia[var + Control_N] + d] += w3 * (dxdx + dydy);
						}
					}
				}
			}
		}

#pragma omp critical
		{
			for (const auto&var : hess_compents_record)
				spline_pardiso->a[var.first] += var.second;
			for (const auto&var : grad_compents_record)
				spline_pardiso->rhs[var.first] -= var.second;
		}

	}

#pragma omp parallel for
	for (int i = Sample_N_m; i < Sample_N + Assist_N; i++)
	{

		vector<pair<int, double>> hess_compents_record;
		int var_num_per_sample = 2 * (degx + 1)*(degy + 1);
		hess_compents_record.reserve(var_num_per_sample*(var_num_per_sample + 1) / 2);

		double area_now = d_.scaf_weight;
		double w1 = 4.*area_now;
		double w2 = 0.;
		double w3 = 4.*area_now;

		std::vector<double> coeff_u, coeff_v;
		int grid_x, grid_y;
		spline_surface->Derivative(samplepointsx[i], samplepointsy[i], grid_x, grid_y, coeff_u, coeff_v);

		for (int i1 = 0; i1 < degx + 1; i1++)
		{
			for (int j1 = 0; j1 < degy + 1; j1++)
			{
				int id = (i1 + grid_x - degx) * n + (j1 + grid_y - degy);
				int var = cntrl2index[id];
				if (var == -1)
					continue;
				int l_var_n1 = i1 * (degy + 1) + j1;

				for (int i2 = 0; i2 < degx + 1; i2++)
				{
					for (int j2 = 0; j2 < degy + 1; j2++)
					{
						int id2 = (i2 + grid_x - degx) * n + (j2 + grid_y - degy);
						int var2 = cntrl2index[id2];
						if (var2 == -1)
							continue;

						if (var2 >= var)
						{
							int l_var_n2 = i2 * (degy + 1) + j2;
							double dx2dy2 = coeff_u[l_var_n1] * coeff_u[l_var_n2] + coeff_v[l_var_n1] * coeff_v[l_var_n2];
							int d = find_cntrl_in_rows.coeff(var, var2);
							hess_compents_record.emplace_back(spline_pardiso->ia[var] + d - 1, w1 * dx2dy2);
							hess_compents_record.emplace_back(spline_pardiso->ia[var + Control_N] + d - 1, w3 * dx2dy2);
						}
					}
				}
			}
		}
#pragma omp critical
		{
			for (const auto&var : hess_compents_record)
				spline_pardiso->a[var.first] += var.second;
		}
	}

	clock_t t22 = clock();
	printf("SamplePoints: %d ; spline-assmbly: %f ;", Sample_N + Assist_N, (t22 - t11) / 1000.0);
	LOG(INFO)<< "SamplePoints: "<< Sample_N + Assist_N <<" ; spline-assmbly: "<< (t22 - t11) / 1000.0;
	clear_samples();
	spline_pardiso->factorize();
	spline_pardiso->pardiso_solver();

	VectorXd d(2 * total_cntrl);
	d.setZero();
#pragma omp parallel for
	for (int i = 0; i < Control_N; i++)
	{
		int id = var_cntrls[i];
		d(id) = spline_pardiso->result[i];
		d(id + total_cntrl) = spline_pardiso->result[i + Control_N];
	}

	double temp_t;
	double alpha;
	pos_of_control += d;

	spline_pardiso->free_numerical_factorization_memory();

	Eigen::MatrixXd direction;
	Reconstruct(direction);
	Eigen::MatrixXd negative_grad_norm;
	calc_gradient_norm(negative_grad_norm);
	double neg_grad_dot_d = -(negative_grad_norm.col(0).dot(direction.col(0)));
	neg_grad_dot_d -= (negative_grad_norm.col(1).dot(direction.col(1)));
	//cout << neg_grad_dot_d << endl;
	if (neg_grad_dot_d < 0)
	{
		double temp_t1;
		max_step(d_.w_uv, direction, temp_t1);
		alpha = 0.95 * temp_t1;
		backtracking_line_search(d_.w_uv, direction, neg_grad_dot_d, alpha);
		d_.w_uv += alpha * direction;
	}

	//Energysource();
	iter_general_cout++;
}
void Parafun::BPE_spline(bool is_slim)
{
	enter_this_solver();
	total_num = d_.w_uv.rows();

	compute_energy(d_.w_uv, false);
	double energy_pre = d_.Energy_cur;
	double energy_cur = energy_pre;

	int slim_iter_num = 0;
	int cm_iter_num = 0;

	double conv_percent = 1;

	long time_beg, time_end;
	time_beg = clock();
	if (is_slim)
	{
		while (iter_general_cout < MAX_ITER_NUM && conv_percent > convgence_con_rate)
		{
			printf("---------------------SplineMap-----iter: %d ------------------------\n", iter_general_cout);
			energy_pre = energy_cur;
			d_.mesh_improve(true);
			after_mesh_improve();
			clock_t t11 = clock();
			SLIM_spline();
			clock_t t22 = clock();
			printf("spline-t: %f s\n", (t22 - t11) / 1000.0);
			LOG(INFO) << "spline-1itertime: " << (t22 - t11) / 1000.0;
			slim_iter_num++;
			energy_cur = compute_energy(d_.w_uv, false);
			d_.energy_record.push_back(d_.Energy_cur);

			d_.energy_record3.emplace_back(d_.Energy_cur, (t22 - t11) / 1000.0,0);

			cout << "spline-energy: " << energy_cur << " energy decrease: " << energy_pre - energy_cur << endl;
			LOG(INFO) << "spline-energy: " << energy_cur << " energy decrease: " << energy_pre - energy_cur;
			conv_percent = abs(energy_cur - energy_pre) / energy_pre;
		}
	}
	else
	{
		conv_percent = 1;
		while (iter_general_cout < MAX_ITER_NUM && conv_percent > convgence_con_rate)
		{
			printf("---------------------SplineMap-----iter: %d --------------------------\n", iter_general_cout);
			energy_pre = energy_cur;
			d_.mesh_improve(true);
			after_mesh_improve();

			clock_t t11 = clock();
			CM_spline();
			clock_t t22 = clock();
			printf("spline-t: %f s\n", (t22 - t11) / 1000.0);
			LOG(INFO) << "spline-1itertime: " << (t22 - t11) / 1000.0;
			cm_iter_num++;
			energy_cur = compute_energy(d_.w_uv, false);
			d_.energy_record.push_back(d_.Energy_cur);
			d_.energy_record3.emplace_back(d_.Energy_cur, (t22 - t11) / 1000.0, 0);

			cout << "spline-energy: " << energy_cur << " energy decrease: " << energy_pre - energy_cur << endl;
			LOG(INFO) << "spline-energy: " << energy_cur << " energy decrease: " << energy_pre - energy_cur;
			conv_percent = abs(energy_cur - energy_pre) / energy_pre;
		}
	}

	time_end = clock();
	time_consumption = (time_end - time_beg) / 1000.0;

	cout << "time_consumption: " << time_consumption << " s; sum_iter: " << iter_general_cout << endl;
	LOG(INFO) << "time_consumption: " << time_consumption << " s; sum_iter: " << iter_general_cout << endl;
	return;

	delete spline_pardiso;
	spline_pardiso = NULL;

	leave_this_solver();
}

void Parafun::iteratively_BPE_spline()
{
	if (set_ctrl_num == 1200)
	{
		//100-3
		init(100);
		convgence_con_rate = 1e-2;
		BPE_spline();
		//200-3
		init(200);
		convgence_con_rate = 1e-3;
		BPE_spline();

		BPE_spline(false);
	}
}

void Parafun::max_step(const Eigen::MatrixXd &x, const Eigen::MatrixXd &d, double &step)
{
	double temp_t = numeric_limits<double>::infinity();
	double a, b, c, b1, b2, tt, tt1, tt2;
	double x0, x1, x2, x3, x4, x5, d0, d1, d2, d3, d4, d5;

	for (int i = 0; i < d_.mf_num; ++i)
	{
		int f0 = d_.F(i, 0);
		int f1 = d_.F(i, 1);
		int f2 = d_.F(i, 2);

		x0 = x(f0); x1 = x(f1); x2 = x(f2); x3 = x(f0 + total_num); x4 = x(f1 + total_num); x5 = x(f2 + total_num);
		d0 = d(f0); d1 = d(f1); d2 = d(f2); d3 = d(f0 + total_num); d4 = d(f1 + total_num); d5 = d(f2 + total_num);

		//x0 = x[f0]; x1 = x[f1]; x2 = x[f2]; x3 = x[f0 + total_num]; x4 = x[f1 + total_num]; x5 = x[f2 + total_num];
		//d0 = d[f0]; d1 = d[f1]; d2 = d[f2]; d3 = d[f0 + total_num]; d4 = d[f1 + total_num]; d5 = d[f2 + total_num];

		a = (d1 - d0) * (d5 - d3) - (d4 - d3) * (d2 - d0);
		b1 = (d1 - d0) * (x5 - x3) + (x1 - x0) * (d5 - d3);
		b2 = (x4 - x3) * (d2 - d0) + (x2 - x0) * (d4 - d3);
		b = b1 - b2;
		c = (x1 - x0) * (x5 - x3) - (x4 - x3) * (x2 - x0);
		tt = get_smallest_pos_quad_zero(a, b, c);

		if (temp_t > tt)
		{
			temp_t = tt;
		}

	}

	for (int i = 0; i < d_.sf_num; ++i)
	{
		int f0 = d_.s_T(i, 0);
		int f1 = d_.s_T(i, 1);
		int f2 = d_.s_T(i, 2);

		x0 = x(f0); x1 = x(f1); x2 = x(f2); x3 = x(f0 + total_num); x4 = x(f1 + total_num); x5 = x(f2 + total_num);
		d0 = d(f0); d1 = d(f1); d2 = d(f2); d3 = d(f0 + total_num); d4 = d(f1 + total_num); d5 = d(f2 + total_num);

		//x0 = x[f0]; x1 = x[f1]; x2 = x[f2]; x3 = x[f0 + total_num]; x4 = x[f1 + total_num]; x5 = x[f2 + total_num];
		//d0 = d[f0]; d1 = d[f1]; d2 = d[f2]; d3 = d[f0 + total_num]; d4 = d[f1 + total_num]; d5 = d[f2 + total_num];

		a = (d1 - d0) * (d5 - d3) - (d4 - d3) * (d2 - d0);
		b1 = (d1 - d0) * (x5 - x3) + (x1 - x0) * (d5 - d3);
		b2 = (x4 - x3) * (d2 - d0) + (x2 - x0) * (d4 - d3);
		b = b1 - b2;
		c = (x1 - x0) * (x5 - x3) - (x4 - x3) * (x2 - x0);
		tt = get_smallest_pos_quad_zero(a, b, c);

		if (temp_t > tt)
		{
			temp_t = tt;
		}

	}
	step = temp_t;
}

void Parafun::backtracking_line_search(const Eigen::MatrixXd &x, const Eigen::MatrixXd &d, const double &tt, double &alpha)
{
	double h = 0.5;
	double c = 0.2;

	double ex = compute_energy(x, true);
	Eigen::MatrixXd x_new = x + alpha * d;
	double e = compute_energy(x_new, true);
	while (e > ex + alpha * c * tt)
	{
		alpha = h * alpha;
		x_new = x + alpha * d;
		e = compute_energy(x_new, true);
	}
}

double Parafun::compute_energy(const Eigen::MatrixXd & x, bool whole)
{
	double total_e_area = 0;
	const double *pos_cur = x.data();

	//cout << "compute energy " << end_e_area << endl;
	if (whole)
	{
#pragma omp parallel for
		for (int i = 0; i < d_.mf_num; ++i)
		{
			int f0 = d_.F(i, 0);
			int f1 = d_.F(i, 1);
			int f2 = d_.F(i, 2);

			double q00 = pos_cur[f0] - pos_cur[f2];
			double q01 = pos_cur[f1] - pos_cur[f2];
			double q10 = pos_cur[f0 + total_num] - pos_cur[f2 + total_num];
			double q11 = pos_cur[f1 + total_num] - pos_cur[f2 + total_num];

			double j00 = d_.upd_inv00[i] * q00 + d_.upd_inv10[i] * q01;
			double j01 = d_.upd_inv01[i] * q00 + d_.upd_inv11[i] * q01;
			double j10 = d_.upd_inv00[i] * q10 + d_.upd_inv10[i] * q11;
			double j11 = d_.upd_inv01[i] * q10 + d_.upd_inv11[i] * q11;

			double det = j00 * j11 - j01 * j10;
			double E_1 = (j00*j00 + j01 * j01 + j10 * j10 + j11 * j11);
			double E_2 = E_1 / (det*det);
#pragma omp critical 
			{
				total_e_area += ((E_1 + E_2)*d_.src_weight[i]);
			}

		}
#pragma omp parallel for
		for (int i = 0; i < d_.sf_num; ++i)
		{
			int f0 = d_.s_T(i, 0);
			int f1 = d_.s_T(i, 1);
			int f2 = d_.s_T(i, 2);

			double q00 = pos_cur[f0] - pos_cur[f2];
			double q01 = pos_cur[f1] - pos_cur[f2];
			double q10 = pos_cur[f0 + total_num] - pos_cur[f2 + total_num];
			double q11 = pos_cur[f1 + total_num] - pos_cur[f2 + total_num];

			double j00 = d_.scaf_jac[i][0] * q00;
			double j01 = d_.scaf_jac[i][1] * q00 + d_.scaf_jac[i][3] * q01;
			double j10 = d_.scaf_jac[i][0] * q10;
			double j11 = d_.scaf_jac[i][1] * q10 + d_.scaf_jac[i][3] * q11;

			double det = j00 * j11 - j01 * j10;

			double E_1 = (j00*j00 + j01 * j01 + j10 * j10 + j11 * j11);
			double E_2 = E_1 / (det*det);
#pragma omp critical 
			{
				total_e_area += ((E_1 + E_2)*d_.scaf_weight);
			}
		}
	}
	else
	{
		//max_energy = 0;
		for (int i = 0; i < d_.mf_num; ++i)
		{
			int f0 = d_.F(i, 0);
			int f1 = d_.F(i, 1);
			int f2 = d_.F(i, 2);

			double q00 = pos_cur[f0] - pos_cur[f2];
			double q01 = pos_cur[f1] - pos_cur[f2];
			double q10 = pos_cur[f0 + total_num] - pos_cur[f2 + total_num];
			double q11 = pos_cur[f1 + total_num] - pos_cur[f2 + total_num];

			double j00 = d_.upd_inv00[i] * q00 + d_.upd_inv10[i] * q01;
			double j01 = d_.upd_inv01[i] * q00 + d_.upd_inv11[i] * q01;
			double j10 = d_.upd_inv00[i] * q10 + d_.upd_inv10[i] * q11;
			double j11 = d_.upd_inv01[i] * q10 + d_.upd_inv11[i] * q11;

			u_x_m[i] = j00;	u_y_m[i] = j01; v_x_m[i] = j10;	v_y_m[i] = j11;

			double det = j00 * j11 - j01 * j10;
			double E_1 = (j00*j00 + j01 * j01 + j10 * j10 + j11 * j11);
			double E_2 = E_1 / (det*det);

			total_e_area += ((E_1 + E_2)*d_.src_weight[i]);

			d_.distortion[i] = (E_1 + E_2);

		}
		d_.Energy_cur = total_e_area / d_.mesh_measure;
	}

	return total_e_area / d_.mesh_measure;
}

void Parafun::Reconstruct(Eigen::MatrixXd& direction)
{
	auto &ctrlpts = spline_surface->ControlPoints();
	auto m = ctrlpts.size();
	auto n = ctrlpts[0].size();

#pragma omp parallel for
	for (int idx_ = 0; idx_ < total_cntrl; idx_++)
	{
		int i = idx_ / n;
		int j = idx_ - i * n;
		ctrlpts[i][j][0] = pos_of_control(idx_);
		ctrlpts[i][j][1] = pos_of_control(idx_ + total_cntrl);
	}
	//for (size_t i = 0; i < m; i++)
	//{
	//	for (size_t j = 0; j < n; j++)
	//	{
	//		ctrlpts[i][j][0] = pos_of_control(i * n + j);
	//		ctrlpts[i][j][1] = pos_of_control(i * n + j + total_cntrl);
	//	}
	//}

	Eigen::MatrixXd position_of_mesh_spline;
	position_of_mesh_spline.resize(total_num, 2);
#pragma omp parallel for
	for (int i = 0; i < total_num; ++i)
	{
		Vector3d p = spline_surface->getPoint(d_.w_uv(i, 0), d_.w_uv(i, 1));
		position_of_mesh_spline(i, 0) = p[0];
		position_of_mesh_spline(i, 1) = p[1];
	}

	direction = position_of_mesh_spline - d_.w_uv;
}

void Parafun::calc_gradient_norm(Eigen::MatrixXd& negative_grad_norm)
{
	negative_grad_norm.resize(total_num, 2);
	negative_grad_norm.setZero();

#pragma omp parallel for
	for (int i = 0; i < d_.mf_num; ++i)
	{
		double area_now = d_.src_weight[i];

		int f0 = d_.F(i, 0);
		int f1 = d_.F(i, 1);
		int f2 = d_.F(i, 2);

		double p00 = d_.upd_inv00[i];
		double p01 = d_.upd_inv01[i];
		double p10 = d_.upd_inv10[i];
		double p11 = d_.upd_inv11[i];

		double j00 = u_x_m[i];
		double j01 = u_y_m[i];
		double j10 = v_x_m[i];
		double j11 = v_y_m[i];
		double det = j00 * j11 - j01 * j10;
		double tr = (j00*j00 + j01 * j01 + j10 * j10 + j11 * j11);

		double det_2plus1 = 2.0 / det / det + 2.0;
		double det_3tr2 = 2.0*tr / det / det / det;
		double r0 = area_now * (j00 * det_2plus1 - det_3tr2 * j11);
		double r1 = area_now * (j01 * det_2plus1 + det_3tr2 * j10);
		double r2 = area_now * (j10 * det_2plus1 + det_3tr2 * j01);
		double r3 = area_now * (j11 * det_2plus1 - det_3tr2 * j00);

#pragma omp critical
		{
			negative_grad_norm(f0) -= (r0*p00 + r1 * p01);
			negative_grad_norm(f1) -= (r0*p10 + r1 * p11);
			negative_grad_norm(f2) += (r0*(p00 + p10) + r1 * (p01 + p11));
			negative_grad_norm(f0 + total_num) -= (r2*p00 + r3 * p01);
			negative_grad_norm(f1 + total_num) -= (r2*p10 + r3 * p11);
			negative_grad_norm(f2 + total_num) += (r2*(p00 + p10) + r3 * (p01 + p11));

			//negative_grad_norm[f0] -= (r0*p00 + r1 * p01);
			//negative_grad_norm[f1] -= (r0*p10 + r1 * p11);
			//negative_grad_norm[f2] += (r0*(p00 + p10) + r1 * (p01 + p11));
			//negative_grad_norm[f0 + total_num] -= (r2*p00 + r3 * p01);
			//negative_grad_norm[f1 + total_num] -= (r2*p10 + r3 * p11);
			//negative_grad_norm[f2 + total_num] += (r2*(p00 + p10) + r3 * (p01 + p11));
		}

	}
}

void Parafun::leave_this_solver()
{
	clear_samples();
	spline_pardiso->a.clear();
	spline_pardiso->rhs.clear();
	u_x_m.clear();
	u_y_m.clear();
	v_x_m.clear();
	v_y_m.clear();
	pos_of_control.resize(0);
}

void Parafun::enter_this_solver()
{
	int mf_n = d_.mf_num;
	u_x_m.resize(mf_n);
	u_y_m.resize(mf_n);
	v_x_m.resize(mf_n);
	v_y_m.resize(mf_n);
	pos_of_control.resize(total_cntrl * 2);

}

void Parafun::clear_samples()
{
	samplepointsx.clear();
	samplepointsy.clear();
	sample_at_face.clear();
}
