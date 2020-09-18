#include "SimplifyMesh.h"



SimplifyMesh::SimplifyMesh(SharedData& sharedd_):sD(sharedd_),mesh(sharedd_.mesh)
{
	V_N = mesh.n_vertices();
	F_N = mesh.n_faces();
}

SimplifyMesh::~SimplifyMesh()
{
}

void SimplifyMesh::prepare_for_dijkstra()
{
	e_len.resize(mesh.n_edges(),1);
}

void SimplifyMesh::color_src_mesh()
{
	prepare_for_dijkstra();
	vector<AbstractEdge> abes;
	vector<int> f_colors;
	f_colors.resize(F_N, 0);

	//color from 1 to 40
	vector<int> f_color_num(100, 0);
	int count_ = 0;

	int f_hope_num = min(F_N / 16,1000000);
	cout << "-----------------f_hope_num------------------- " << f_hope_num << endl;
	color_disk(sD.boundary_es_del, f_colors, count_, 0, F_N, f_hope_num);

	f_color_num.assign(count_ + 1, 0);
	for (int i = 0; i < f_colors.size(); i++)
	{
		f_color_num[f_colors[i]]++;
	}

	for (size_t i = 0; i < f_color_num.size(); i++)
	{
		cout << i << " : " << f_color_num[i] << endl;
	}

	vector<AFace> Afs;
	bool is_handle_ = true;
	set<int> b_nodes;
	while (is_handle_)
	{
		cout << "---------------outer---------------------" << endl;
		find_inter_line(abes, b_nodes, f_colors, sD.boundary_es_del);
		is_handle_ = handle_encircle_case(Afs, abes, f_colors, 0, f_color_num, f_hope_num);
	}
	cout << "---------------curve_to_line---------------------" << endl;

	curve_to_line(Afs, abes, f_colors, f_color_num);

	is_handle_ = true;
	while (is_handle_)
	{
		cout << "---------------outer---------------------" << endl;
		find_inter_line(abes, b_nodes, f_colors, sD.boundary_es_del);
		is_handle_ = handle_encircle_case(Afs, abes, f_colors, 0, f_color_num, f_hope_num);
	}

	split_for_tutte(abes, f_colors, f_color_num);

	tutte_the_parts(Afs, abes, b_nodes, f_colors, f_color_num);
	e_len.clear();
	cout << "-----------------------------finish------- " << endl;


}

void SimplifyMesh::color_disk(const vector<int>& b_hes, vector<int>& f_cs, int & color_start, int color_split, int f_total_num, int f_part_hope, pair<int, int> b_s_e_idx)
{
	if (b_hes.size() < 3)
	{
		cout << "Wrong input!!!" << b_hes.size() << endl;
		return;
	}

	int part_num = ceil((double)f_total_num / f_part_hope);
	if (part_num <= 1)
		return;


	int color_count_ = color_start;
	if (color_split == 0)
	{
		part_num = max(24, min(part_num, 48));
	}
	else
	{
		part_num = max(2, min(part_num, 8));
	}
	vector<int> f_c_num(color_count_ + part_num + 2, 0);

	queue<pair<int, int>> Q;

	{
		int b_he_num = b_hes.size();
		int part_he_num = ceil((double)b_he_num / part_num);
		int p_c_ = 0;
		int backid = b_hes.back();

		for (int j = 0; j < b_he_num; j++)
		{
			if (j > 0)
				backid = b_hes[j - 1];
			auto he_h = mesh.halfedge_handle(b_hes[j]);
			if (j >= p_c_ * part_he_num)
			{
				color_count_++;
				p_c_++;
			}
			auto f_h = mesh.face_handle(he_h);
			if (f_cs[f_h.idx()] == color_split)
			{
				f_cs[f_h.idx()] = color_count_;
				f_c_num[color_count_]++;
				Q.emplace(f_h.idx(), color_count_);
			}
		}
		cout << "p1_c_: " << p_c_ << endl;
	}



	while (!Q.empty())
	{
		auto node_ = Q.front();
		Q.pop();
		auto seed_f_h = mesh.face_handle(node_.first);
		for (auto itvf = mesh.ff_begin(seed_f_h); itvf != mesh.ff_end(seed_f_h); itvf++)
		{
			if (f_cs[itvf->idx()] == color_split)
			{
				f_cs[itvf->idx()] = node_.second;
				f_c_num[node_.second]++;
				Q.emplace(itvf->idx(), node_.second);
			}
		}
	}


	vector<AFace> Afs;
	vector<AbstractEdge> abes_p;
	bool is_handle_ = true;
	set<int> b_nodes_;
	while (is_handle_)
	{
		cout << "------------------------------------" << endl;
		find_inter_line(abes_p, b_nodes_, f_cs, b_hes);
		is_handle_ = handle_encircle_case(Afs, abes_p, f_cs, color_start, f_c_num, f_part_hope);
	}
	cout << "----------------handle_over--------------------" << endl;

	int f_split_num = 0;
	for (int i = 0; i < f_cs.size(); i++)
	{
		if (f_cs[i] == color_split)
			f_split_num++;
	}
	if (f_split_num != 0)
	{
		cout << "Split faces num Wrong!!!" << endl;
	}
	int new_color_beginning_idx = color_start + 1;
	color_start = color_count_;
	//	return;
	for (int c_i = new_color_beginning_idx; c_i < f_c_num.size(); c_i++)
	{
		if (f_c_num[c_i] > f_part_hope)
		{
			// recursively
			int af_i = 0;
			for (; af_i < Afs.size(); af_i++)
			{
				if (c_i == Afs[af_i].f_color)
				{
					break;
				}
			}
			//vector<vector<int>> c_b_he_lines;
			vector<int> c_b_hes;

			auto & ae_ids = Afs[af_i].ae_ids;
			vector<bool> ae_is_visited_(ae_ids.size(), false);

			int idx_ae_ids = -1;
			for (int i = 0; i < ae_ids.size(); i++)
			{
				if (ae_ids[i].first)
				{
					auto v_h_ = mesh.vertex_handle(abes_p[ae_ids[i].second].start_id);
					if (!mesh.is_boundary(v_h_))
					{
						idx_ae_ids = i;
						break;
					}
				}
				else
				{
					auto v_h_ = mesh.vertex_handle(abes_p[ae_ids[i].second].end_id);
					if (!mesh.is_boundary(v_h_))
					{
						idx_ae_ids = i;
						break;
					}
				}
			}

			int link_node_id = -1;
			while (idx_ae_ids != -1)
			{
				if (ae_ids[idx_ae_ids].first)
				{
					//c_b_he_lines.push_back(abes_p[ae_ids[idx_ae_ids].second].he_line);

					auto &he_line_ = abes_p[ae_ids[idx_ae_ids].second].he_line;
					for (int i = 0; i < he_line_.size(); i++)
					{
						c_b_hes.push_back(he_line_[i]);
					}
					link_node_id = abes_p[ae_ids[idx_ae_ids].second].end_id;
				}
				else
				{
					auto & he_line_origin = abes_p[ae_ids[idx_ae_ids].second].he_line;
					//vector<int> he_line_opp(he_line_origin.size());
					//for (int i = he_line_origin.size() - 1; i >= 0; i--)
					//{
					//	he_line_opp[he_line_origin.size() - 1 - i] = mesh.opposite_halfedge_handle(mesh.halfedge_handle(he_line_origin[i])).idx();
					//}
					//c_b_he_lines.emplace_back(he_line_opp);

					for (int i = he_line_origin.size() - 1; i >= 0; i--)
					{
						c_b_hes.push_back(mesh.opposite_halfedge_handle(mesh.halfedge_handle(he_line_origin[i])).idx());
					}

					link_node_id = abes_p[ae_ids[idx_ae_ids].second].start_id;
				}
				ae_is_visited_[idx_ae_ids] = true;
				idx_ae_ids = -1;
				for (int i = 0; i < ae_ids.size(); i++)
				{
					if (!ae_is_visited_[i])
					{
						if (ae_ids[i].first)
						{
							if (link_node_id == abes_p[ae_ids[i].second].start_id)
							{
								idx_ae_ids = i;
								break;
							}
						}
						else
						{
							if (link_node_id == abes_p[ae_ids[i].second].end_id)
							{
								idx_ae_ids = i;
								break;
							}
						}
					}
				}
			}
			cout << "enter color: " << c_i << " with " << c_b_hes.size() << " boundary parts " << f_c_num[c_i] << endl;

			pair<int, int> b_s_e_idx_ = pair<int, int>(0, 0);
			bool is_find_first_ = false;

			for (int i = 0; i < c_b_hes.size(); i++)
			{
				auto he_h_ = mesh.halfedge_handle(c_b_hes[i]);
				int f0 = mesh.opposite_face_handle(he_h_).idx();
				if (f0 == -1)
				{
					if (!is_find_first_)
					{
						b_s_e_idx_.first = i;
						is_find_first_ = true;
					}
				}
				else
				{
					if (is_find_first_)
					{
						b_s_e_idx_.second = i;
						break;
					}
				}
			}

			color_disk(c_b_hes, f_cs, color_start, c_i, f_c_num[c_i], f_part_hope, b_s_e_idx_);
		}
	}

}

bool SimplifyMesh::handle_encircle_case(vector<AFace>& Afs, vector<AbstractEdge>& abes_p, vector<int>& f_cs, int color_start_, vector<int>& f_color_num, int f_part_hope)
{
	bool res_ = false;
	Afs.clear();
	for (int k = 0; k < abes_p.size(); k++)
	{
		int j_l = -1;
		int j_r = -1;
		for (int i = 0; i < Afs.size(); i++)
		{
			if (abes_p[k].color_left == Afs[i].f_color)
			{
				Afs[i].ae_ids.emplace_back(true, k);
				j_l = i;
				break;
			}
		}
		if (j_l == -1)
		{
			vector<pair<bool, int>> vec_tmp = { pair<bool,int>(true,k) };
			Afs.emplace_back(abes_p[k].color_left, vec_tmp);
		}

		for (int i = 0; i < Afs.size(); i++)
		{
			if (abes_p[k].color_right == Afs[i].f_color)
			{
				Afs[i].ae_ids.emplace_back(false, k);
				j_r = i;
				break;
			}
		}
		if (j_r == -1)
		{
			vector<pair<bool, int>> vec_tmp = { pair<bool,int>(false,k) };
			Afs.emplace_back(abes_p[k].color_right, vec_tmp);
		}
	}

	for (size_t i = 0; i < Afs.size(); i++)
	{
		cout << "color: " << Afs[i].f_color << " es: " << Afs[i].ae_ids.size();
		if (Afs[i].f_color != -1)
			cout << " fs: " << f_color_num[Afs[i].f_color];
		cout << endl;
	}

	//sort the ae_ids

	vector<vector<vector<pair<bool, int>>>> ae_ids_sort_all(Afs.size());

	for (int af_i = 0; af_i < Afs.size(); af_i++)
	{
		if (Afs[af_i].f_color == -1)
			continue;

		auto & ae_ids = Afs[af_i].ae_ids;
		vector<bool> ae_is_visited_(ae_ids.size(), false);
		vector<vector<pair<bool, int>>>& ae_ids_sort = ae_ids_sort_all[af_i];

		int idx_ae_ids = -1;
		for (int i = 0; i < ae_ids.size(); i++)
		{
			if (ae_ids[i].first)
			{
				auto v_h_ = mesh.vertex_handle(abes_p[ae_ids[i].second].start_id);
				if (!mesh.is_boundary(v_h_))
				{
					idx_ae_ids = i;
					ae_ids_sort.emplace_back(vector<pair<bool, int>>{});
					break;
				}
			}
			else
			{
				auto v_h_ = mesh.vertex_handle(abes_p[ae_ids[i].second].end_id);
				if (!mesh.is_boundary(v_h_))
				{
					idx_ae_ids = i;
					ae_ids_sort.emplace_back(vector<pair<bool, int>>{});
					break;
				}
			}
		}

		//	int link_node_id = -1;
		while (idx_ae_ids != -1)
		{
			OpenMesh::HalfedgeHandle he_guide;
			if (ae_ids[idx_ae_ids].first)
			{
				auto &he_line_ = abes_p[ae_ids[idx_ae_ids].second].he_line;
				he_guide = mesh.opposite_halfedge_handle(mesh.halfedge_handle(he_line_.back()));
				//	link_node_id = abes_p[ae_ids[idx_ae_ids].second].end_id;
			}
			else
			{
				auto & he_line_origin = abes_p[ae_ids[idx_ae_ids].second].he_line;
				he_guide = mesh.halfedge_handle(he_line_origin.front());
				//	link_node_id = abes_p[ae_ids[idx_ae_ids].second].start_id;
			}
			ae_is_visited_[idx_ae_ids] = true;
			ae_ids_sort.back().push_back(ae_ids[idx_ae_ids]);
			int backid = he_guide.idx();
			do
			{
				auto he_opp = mesh.opposite_halfedge_handle(he_guide);
				int f0 = mesh.face_handle(he_opp).idx();
				if (f0 == -1 || f_cs[f0] != Afs[af_i].f_color)
				{
					break;
				}
				he_guide = mesh.next_halfedge_handle(he_opp);
			} while (he_guide.idx() != backid);

			idx_ae_ids = -1;
			for (int i = 0; i < ae_ids.size(); i++)
			{
				if (!ae_is_visited_[i])
				{
					if (ae_ids[i].first)
					{
						if (he_guide.idx() == abes_p[ae_ids[i].second].he_line.front())
						{
							idx_ae_ids = i;
							break;
						}
					}
					else
					{
						int he_check_ = mesh.opposite_halfedge_handle(mesh.halfedge_handle(abes_p[ae_ids[i].second].he_line.back())).idx();
						if (he_guide.idx() == he_check_)
						{
							idx_ae_ids = i;
							break;
						}
					}
				}
			}

			if (idx_ae_ids == -1)
			{
				for (int i = 0; i < ae_ids.size(); i++)
				{
					if (!ae_is_visited_[i])
					{
						if (ae_ids[i].first)
						{
							auto v_h_ = mesh.vertex_handle(abes_p[ae_ids[i].second].start_id);
							if (!mesh.is_boundary(v_h_))
							{
								idx_ae_ids = i;
								ae_ids_sort.emplace_back(vector<pair<bool, int>>{});
								break;
							}
						}
						else
						{
							auto v_h_ = mesh.vertex_handle(abes_p[ae_ids[i].second].end_id);
							if (!mesh.is_boundary(v_h_))
							{
								idx_ae_ids = i;
								ae_ids_sort.emplace_back(vector<pair<bool, int>>{});
								break;
							}
						}
					}
				}
			}

		}

	}


	vector<bool> af_is_set_(Afs.size(), false);

	for (int af_i = 0; af_i < Afs.size(); af_i++)
	{
		if (Afs[af_i].f_color == -1)
			continue;

		if (af_is_set_[af_i])
			continue;

		vector<vector<pair<bool, int>>>& ae_ids_sort = ae_ids_sort_all[af_i];

		if (ae_ids_sort.size() > 1)
		{
			int max_long_id = -1;
			int max_long = 0;
			for (int i = 0; i < ae_ids_sort.size(); i++)
			{
				int long_size_ = 0;
				for (const auto&var_ : ae_ids_sort[i])
				{
					long_size_ += abes_p[var_.second].he_line.size();
				}
				if (long_size_ > max_long)
				{
					max_long = long_size_;
					max_long_id = i;
				}
			}
			//Afs[af_i].ae_ids = ae_ids_sort[max_long_id];
			for (int i = 0; i < ae_ids_sort.size(); i++)
			{
				if (ae_ids_sort[i].size() < 2)
				{
					cout << "impossible: " << ae_ids_sort[i].size() << endl;
				}
				if (i == max_long_id)
					continue;

				int long_side_id = -1;
				int long_side = 0;
				int c_to = -1;
				for (int j = 0; j < ae_ids_sort[i].size(); j++)
				{
					if (abes_p[ae_ids_sort[i][j].second].color_left != -1 && abes_p[ae_ids_sort[i][j].second].color_right != -1)
					{
						if (abes_p[ae_ids_sort[i][j].second].color_left != abes_p[ae_ids_sort[i][j].second].color_right)
							if (long_side < abes_p[ae_ids_sort[i][j].second].he_line.size())
							{
								long_side = abes_p[ae_ids_sort[i][j].second].he_line.size();
								long_side_id = j;
								if (ae_ids_sort[i][j].first)
								{
									c_to = abes_p[ae_ids_sort[i][j].second].color_right;
								}
								else
								{
									c_to = abes_p[ae_ids_sort[i][j].second].color_left;
								}
							}
					}
				}
				if (-1 != c_to)
				{
					queue<int> Q;
					int c_to_add_num = 0;
					for (const auto&var_ : ae_ids_sort[i])
					{
						if (var_.first)
						{
							for (const auto&heid_ : abes_p[var_.second].he_line)
							{
								auto he_h_ = mesh.halfedge_handle(heid_);
								int f0 = mesh.face_handle(he_h_).idx();
								if (f_cs[f0] == Afs[af_i].f_color)
								{
									f_cs[f0] = c_to;
									c_to_add_num++;
									Q.push(f0);
								}
							}
						}
						else
						{
							for (const auto&heid_ : abes_p[var_.second].he_line)
							{
								auto he_h_ = mesh.halfedge_handle(heid_);
								int f0 = mesh.opposite_face_handle(he_h_).idx();
								if (f_cs[f0] == Afs[af_i].f_color)
								{
									f_cs[f0] = c_to;
									c_to_add_num++;
									Q.push(f0);
								}
							}
						}
					}

					while (!Q.empty())
					{
						auto node_ = Q.front();
						Q.pop();
						auto seed_f_h = mesh.face_handle(node_);
						for (auto itvf = mesh.ff_begin(seed_f_h); itvf != mesh.ff_end(seed_f_h); itvf++)
						{
							if (f_cs[itvf->idx()] == Afs[af_i].f_color)
							{
								f_cs[itvf->idx()] = c_to;
								c_to_add_num++;

								Q.emplace(itvf->idx());
							}
						}
					}
					if (c_to_add_num > 0)
					{
						f_color_num[c_to] += c_to_add_num;
						f_color_num[Afs[af_i].f_color] -= c_to_add_num;
						abes_p[ae_ids_sort[i][long_side_id].second].color_left = c_to;
						abes_p[ae_ids_sort[i][long_side_id].second].color_right = c_to;
					}
					res_ = true;
					cout << "Handle_NODesire_case--- c: " << Afs[af_i].f_color << " num: " << c_to_add_num << endl;
				}
				//int inter_he_id_ = ae_ids_sort[i][long_side_id].second;
				//for (size_t j = 0; j < Afs.size(); j++)
				//{
				//	if (Afs[j].f_color == c_to)
				//	{
				//		break;
				//	}
				//}

			}
			cout << "Handle_NODesire_case--- c: " << Afs[af_i].f_color << "over." << endl;
			//res_ = true;
		}
		af_is_set_[af_i] = true;
	}



	for (int i = 0; i < Afs.size(); i++)
	{
		if (Afs[i].f_color <= color_start_)
			continue;

		if (Afs[i].f_color == -1)
			continue;
		if (f_color_num[Afs[i].f_color] == 0)
			continue;

		if (f_color_num[Afs[i].f_color] < 0.1*f_part_hope)
		{
			cout << "Handle_small_case--- c: " << Afs[i].f_color << " c_num: " << f_color_num[Afs[i].f_color] << endl;

			int max_line_num_ = 0;
			int max_idx_ = -1;
			int max_l_hx_ = 0;
			int max_idx_hx_ = -1;

			for (int j = 0; j < Afs[i].ae_ids.size(); j++)
			{
				const auto& var_ = Afs[i].ae_ids[j];
				if (max_line_num_ < abes_p[var_.second].he_line.size())
				{
					if (abes_p[var_.second].color_left != -1 && abes_p[var_.second].color_right != -1)
					{
						max_l_hx_= abes_p[var_.second].he_line.size();
						max_idx_hx_ = j;

						int two_color_sum_ = f_color_num[abes_p[var_.second].color_left] + f_color_num[abes_p[var_.second].color_right];
						if (two_color_sum_ < f_part_hope)
						{
							max_line_num_ = abes_p[var_.second].he_line.size();
							max_idx_ = j;
						}
					}
				}
			}

			if (max_idx_ == -1)
			{
				if (f_color_num[Afs[i].f_color] < 0.01*f_part_hope)
				{
					max_idx_ = max_idx_hx_;
				}
				else
				{
					continue;
				}
			}

			queue<int> Q;
			int c_to = -1;
			int c_change_num = 0;
			if (Afs[i].ae_ids[max_idx_].first)
			{
				c_to = abes_p[Afs[i].ae_ids[max_idx_].second].color_right;
				const auto&var_ = abes_p[Afs[i].ae_ids[max_idx_].second].he_line;
				for (int j = 0; j < var_.size(); j++)
				{
					auto he_h = mesh.halfedge_handle(var_[j]);
					int f1 = mesh.face_handle(he_h).idx();
					if (f_cs[f1] == Afs[i].f_color)
					{
						f_cs[f1] = c_to;
						Q.push(f1);
						c_change_num++;
					}
				}
			}
			else
			{
				c_to = abes_p[Afs[i].ae_ids[max_idx_].second].color_left;
				const auto&var_ = abes_p[Afs[i].ae_ids[max_idx_].second].he_line;
				for (int j = 0; j < var_.size(); j++)
				{
					auto he_h = mesh.halfedge_handle(var_[j]);
					int f2 = mesh.opposite_face_handle(he_h).idx();
					if (f_cs[f2] == Afs[i].f_color)
					{
						f_cs[f2] = c_to;
						Q.push(f2);
						c_change_num++;
					}
				}
			}

			while (!Q.empty())
			{
				auto node_ = Q.front();
				Q.pop();
				auto seed_f_h = mesh.face_handle(node_);
				for (auto itvf = mesh.ff_begin(seed_f_h); itvf != mesh.ff_end(seed_f_h); itvf++)
				{
					if (f_cs[itvf->idx()] == Afs[i].f_color)//== delete_color
					{
						f_cs[itvf->idx()] = c_to;
						Q.push(itvf->idx());
						c_change_num++;
					}
				}
			}

			if (c_change_num != f_color_num[Afs[i].f_color])
			{
				cout << "delete small part Wrong !!! c_change_num: " << c_change_num << " src_delete_color num: " << f_color_num[Afs[i].f_color] << endl;
			}

			f_color_num[c_to] += f_color_num[Afs[i].f_color];
			f_color_num[Afs[i].f_color] = 0;

			res_ = true;

		}
		else if (Afs[i].ae_ids.size() == 2)
		{
			if (color_start_ > 0)
			{
				set<int> cs_outer_;
				bool is_skip_flag_ = false;
				for (const auto&var_ : Afs[i].ae_ids)
				{
					if (var_.first)
					{
						for (const auto&heid_ : abes_p[var_.second].he_line)
						{
							auto he_h_ = mesh.halfedge_handle(heid_);
							int fid_ = mesh.opposite_face_handle(he_h_).idx();
							if (fid_ == -1)
								cs_outer_.insert(-1);
							else
								cs_outer_.insert(f_cs[fid_]);
						}
					}
					else
					{
						for (const auto&heid_ : abes_p[var_.second].he_line)
						{
							auto he_h_ = mesh.halfedge_handle(heid_);
							int fid_ = mesh.face_handle(he_h_).idx();
							if (fid_ == -1)
								cs_outer_.insert(-1);
							else
								cs_outer_.insert(f_cs[fid_]);
						}
					}

					if (cs_outer_.size() > 2)
					{
						is_skip_flag_ = true;
						break;
					}
				}
				if (is_skip_flag_)
					continue;
			}

			int delete_color = Afs[i].f_color;
			int de1 = Afs[i].ae_ids.front().second;
			int de2 = Afs[i].ae_ids.back().second;

			int c1 = Afs[i].ae_ids.front().first ? abes_p[de1].color_right : abes_p[de1].color_left;
			int c2 = Afs[i].ae_ids.back().first ? abes_p[de2].color_right : abes_p[de2].color_left;

			cout << "Handle_encircle_case--- c: " << Afs[i].f_color << " c_num: " << f_color_num[Afs[i].f_color] << " ;neighbor ( " << c1 << " , " << c2 << " )" << endl;

			vector<int> s_path_;
			int start_id = abes_p[de1].start_id;
			int end_id = abes_p[de1].end_id;

			if (Afs[i].ae_ids.front().first)
			{
				start_id = abes_p[de1].end_id;
				end_id = abes_p[de1].start_id;
			}

			if (-1 == c1||-1==c2)
			{
				int c_to_;
				if (-1 == c1)
				{
					c_to_ = c2;
					cout << "Handle_encircle_case need process!!! add ALL the faces in color: " << delete_color << " to color: " << c2 << endl;
				}
				else
				{
					c_to_ = c1;
					cout << "Handle_encircle_case need process!!! add ALL the faces in color: " << delete_color << " to color: " << c1 << endl;
				}

				queue<int> Q;
				int c_add_to_num_ = 0;
				if (Afs[i].ae_ids.front().first)
				{
					for(auto&heid_: abes_p[de1].he_line)
					{
						int f0 = mesh.face_handle(mesh.halfedge_handle(heid_)).idx();
						if (f_cs[f0] == delete_color)
						{
							Q.push(f0);
							f_cs[f0] = c_to_;
							c_add_to_num_++;
						}
					}
				}
				else
				{
					for (auto&heid_ : abes_p[de1].he_line)
					{
						int f0 = mesh.opposite_face_handle(mesh.halfedge_handle(heid_)).idx();
						if (f_cs[f0] == delete_color)
						{
							Q.push(f0);
							f_cs[f0] = c_to_;
							c_add_to_num_++;
						}
					}
				}

				while (!Q.empty())
				{
					auto node_ = Q.front();
					Q.pop();
					auto seed_f_h = mesh.face_handle(node_);
					for (auto itvf = mesh.ff_begin(seed_f_h); itvf != mesh.ff_end(seed_f_h); itvf++)
					{
						if (f_cs[itvf->idx()] == delete_color)//== delete_color
						{
							f_cs[itvf->idx()] = c_to_;
							Q.push(itvf->idx());
							c_add_to_num_++;
						}
					}
				}
				f_color_num[c_to_] += c_add_to_num_;
				f_color_num[delete_color] -= c_add_to_num_;
				if (c_add_to_num_ != f_color_num[delete_color])
				{
					cout << "Handle_encircle_case process wrong!!!" << endl;
				}
				else
				{
					cout << "Handle_encircle_case process OK." << endl;
				}

				res_ = true;
			}
			else
			{
				//dijkstra find s_path_
				shortest_path(end_id, set<int>{start_id}, f_cs, vector<int>{ c1, delete_color,c2}, s_path_);

				set<int> block_fids_1, block_fids_2;
				if (Afs[i].ae_ids.front().first)
				{
					for (const auto&heid_ : abes_p[de1].he_line)
					{
						auto he_h = mesh.halfedge_handle(heid_);
						int f_id_ = mesh.opposite_face_handle(he_h).idx();
						block_fids_1.insert(f_id_);
					}
				}
				else
				{
					for (const auto&heid_ : abes_p[de1].he_line)
					{
						auto he_h = mesh.halfedge_handle(heid_);
						int f_id_ = mesh.face_handle(he_h).idx();
						block_fids_1.insert(f_id_);
					}
				}

				if (Afs[i].ae_ids.back().first)
				{
					for (const auto&heid_ : abes_p[de2].he_line)
					{
						auto he_h = mesh.halfedge_handle(heid_);
						int f_id_ = mesh.opposite_face_handle(he_h).idx();
						block_fids_2.insert(f_id_);
					}
				}
				else
				{
					for (const auto&heid_ : abes_p[de2].he_line)
					{
						auto he_h = mesh.halfedge_handle(heid_);
						int f_id_ = mesh.face_handle(he_h).idx();
						block_fids_2.insert(f_id_);
					}
				}



				//change f_cs
				queue<int> Q1, Q2;
				int c1_num = 0;
				int c2_num = 0;
				for (int j = 0; j < s_path_.size(); j++)
				{
					auto he_h = mesh.halfedge_handle(s_path_[j]);
					int f1 = mesh.face_handle(he_h).idx();
					int f2 = mesh.opposite_face_handle(he_h).idx();
					block_fids_2.insert(f1);
					block_fids_1.insert(f2);

					if (-1 != f1 && f_cs[f1] != c1)//== delete_color
					{
						if (f_cs[f1] == c2)
							c2_num--;
						f_cs[f1] = c1;
						Q1.push(f1);
						c1_num++;
					}
					if (-1 != f2 && f_cs[f2] != c2)//== delete_color
					{
						if (f_cs[f2] == c1)
							c1_num--;
						f_cs[f2] = c2;
						Q2.push(f2);
						c2_num++;
					}
				}

				while (!Q1.empty())
				{
					auto node_ = Q1.front();
					Q1.pop();
					auto seed_f_h = mesh.face_handle(node_);
					for (auto itvf = mesh.ff_begin(seed_f_h); itvf != mesh.ff_end(seed_f_h); itvf++)
					{
						if (block_fids_1.find(itvf->idx()) == block_fids_1.end() && f_cs[itvf->idx()] != c1)//== delete_color
						{
							if (f_cs[itvf->idx()] == c2)
								c2_num--;
							f_cs[itvf->idx()] = c1;
							c1_num++;
							Q1.push(itvf->idx());
						}
					}
				}

				while (!Q2.empty())
				{
					auto node_ = Q2.front();
					Q2.pop();
					auto seed_f_h = mesh.face_handle(node_);
					for (auto itvf = mesh.ff_begin(seed_f_h); itvf != mesh.ff_end(seed_f_h); itvf++)
					{
						if (block_fids_2.find(itvf->idx()) == block_fids_2.end() && f_cs[itvf->idx()] != c2)//== delete_color
						{
							if (f_cs[itvf->idx()] == c1)
								c1_num--;
							f_cs[itvf->idx()] = c2;
							c2_num++;
							Q2.push(itvf->idx());
						}
					}
				}

				if ((c1_num + c2_num) != f_color_num[delete_color])
				{
					cout << "s_path_-size: " << s_path_.size() << " start_id: " << start_id << " end_id: " << end_id << endl;
					cout << "Handle_encircle_case Wrong!!! c1: " << c1_num << " c2: " << c2_num << " c12-ori: " << f_color_num[delete_color] << endl;
				}
				else
				{
					cout << "Handle_encircle_case OK... c1: " << c1_num << " c2: " << c2_num << " c12-ori: " << f_color_num[delete_color] << endl;
				}
				if (-1 != c1)
					f_color_num[c1] += c1_num;
				if (-1 != c2)
					f_color_num[c2] += c2_num;
				f_color_num[delete_color] = 0;

				for (int j = 0; j < Afs.size(); j++)
				{
					if (Afs[j].f_color == c1)
					{
						vector<int> new_line_;

						pair<bool, int> k1(false, -1);
						pair<bool, int> k2(false, -1);

						auto& var = Afs[j].ae_ids;

						for (auto it2 = var.begin(); it2 != var.end();)
						{
							if (it2->second == de1)
							{
								it2 = var.erase(it2);
								continue;
							}
							if (it2->first)
							{
								if (abes_p[it2->second].end_id == start_id)
								{
									k1 = *it2;
									it2 = var.erase(it2);
								}
								else if (abes_p[it2->second].start_id == end_id)
								{
									k2 = *it2;
									it2 = var.erase(it2);
								}
								else
								{
									it2++;
								}
							}
							else
							{
								if (abes_p[it2->second].end_id == end_id)
								{
									k2 = *it2;
									it2 = var.erase(it2);
								}
								else if (abes_p[it2->second].start_id == start_id)
								{
									k1 = *it2;
									it2 = var.erase(it2);
								}
								else
								{
									it2++;
								}
							}

						}

						if (k1.second != -1)
						{
							if (k1.first)
							{
								for (const auto&he_ : abes_p[k1.second].he_line)
								{
									new_line_.push_back(he_);
								}
							}
							else
							{
								for (auto itre = abes_p[k1.second].he_line.rbegin(); itre != abes_p[k1.second].he_line.rend(); itre++)
								{
									new_line_.push_back(mesh.opposite_halfedge_handle(mesh.halfedge_handle(*itre)).idx());
								}

							}
						}

						for (const auto&he_ : s_path_)
						{
							new_line_.push_back(he_);
						}

						if (k2.second != -1)
						{
							if (k2.first)
							{
								for (const auto&he_ : abes_p[k2.second].he_line)
								{
									new_line_.push_back(he_);
								}
							}
							else
							{
								for (auto itre = abes_p[k2.second].he_line.rbegin(); itre != abes_p[k2.second].he_line.rend(); itre++)
								{
									new_line_.push_back(mesh.opposite_halfedge_handle(mesh.halfedge_handle(*itre)).idx());
								}
							}
						}

						abes_p.emplace_back(c1, c2, new_line_);
						abes_p.back().set_start_end(mesh.from_vertex_handle(mesh.halfedge_handle(new_line_.front())).idx(),
							mesh.to_vertex_handle(mesh.halfedge_handle(new_line_.back())).idx());

						var.emplace_back(true, abes_p.size());

						for (int j2 = 0; j2 < Afs.size(); j2++)
						{
							if (Afs[j2].f_color == c2)
							{
								auto& var2 = Afs[j2].ae_ids;

								for (auto it2 = var2.begin(); it2 != var2.end();)
								{
									if (it2->second == de2)
									{
										it2 = var2.erase(it2);
									}
									else if (it2->second == k1.second)
									{
										it2 = var2.erase(it2);
									}
									else if (it2->second == k2.second)
									{
										it2 = var2.erase(it2);
									}
									else
									{
										it2++;
									}
								}
								var2.emplace_back(false, abes_p.size());
								break;
							}
						}

						break;
					}
				}
				res_ = true;
			}

		}
	}

	return res_;
}

void SimplifyMesh::find_inter_line(vector<AbstractEdge>& abes_p, set<int>& b_nodes, const vector<int>& f_cs, const vector<int>& b_hes)
{
	//f_color: f_cs ; boundary: b_hes ;======> abes_p ; b_nodes ;
	set<int> e_visited;
	set<int> end_vids;
	vector<pair<int, int>> node_direction;
	abes_p.clear();
	b_nodes.clear();
	{
		int b_hes_num = b_hes.size();
		int start_he_id = -1;
		int backid = b_hes.back();
		bool is_found_ = false;
		int i = 0;
		for (; i < b_hes_num; i++)
		{
			if (i > 0)
				backid = b_hes[i - 1];
			auto he_h = mesh.halfedge_handle(b_hes[i]);
			auto he_h_prev = mesh.prev_halfedge_handle(he_h);
			while (he_h_prev.idx() != backid)
			{
				int f0 = mesh.face_handle(he_h_prev).idx();
				auto he_h_opp = mesh.opposite_halfedge_handle(he_h_prev);
				int f1 = mesh.face_handle(he_h_opp).idx();
				if (f_cs[f0] != f_cs[f1])
				{
					start_he_id = he_h_opp.idx();
					is_found_ = true;
					break;
				}
				he_h_prev = mesh.prev_halfedge_handle(he_h_opp);
			}
			if (is_found_)
				break;
		}
		if (!is_found_)
		{
			cout << "Only one color, abes size: " << abes_p.size() << endl;
			return;
		}

		for (int j1 = 0; j1 < b_hes_num; j1++)
		{
			int j = (i + j1) % b_hes_num;
			backid = b_hes[(j - 1 + b_hes_num) % b_hes_num];
			is_found_ = false;
			auto he_h = mesh.halfedge_handle(b_hes[j]);
			auto he_h_prev = mesh.prev_halfedge_handle(he_h);
			while (he_h_prev.idx() != backid)
			{
				int f0 = mesh.face_handle(he_h_prev).idx();
				auto he_h_opp = mesh.opposite_halfedge_handle(he_h_prev);
				int f1 = mesh.face_handle(he_h_opp).idx();
				if (f_cs[f0] != f_cs[f1])
				{
					int startid_ = mesh.from_vertex_handle(he_h_opp).idx();
					node_direction.emplace_back(startid_, he_h_opp.idx());
					is_found_ = true;
					//break;
				}
				he_h_prev = mesh.prev_halfedge_handle(he_h_opp);
			}

			if (is_found_)
			{
				int f_opp_id_ = mesh.opposite_face_handle(he_h).idx();
				abes_p.emplace_back(f_cs[mesh.face_handle(he_h).idx()], -1, vector<int>{});

				//if (f_opp_id_ == -1)
				//{
				//	abes_p.emplace_back(f_cs[mesh.face_handle(he_h).idx()], -1, vector<int>{});
				//}
				//else
				//{
				//	abes_p.emplace_back(f_cs[mesh.face_handle(he_h).idx()], f_cs[f_opp_id_], vector<int>{});
				//}

				b_nodes.insert(mesh.from_vertex_handle(he_h).idx());
			}

			abes_p.back().he_line.push_back(b_hes[j]);

		}

		for (auto&var : abes_p)
		{
			int s_id_ = mesh.from_vertex_handle(mesh.halfedge_handle(var.he_line.front())).idx();
			int e_id_ = mesh.to_vertex_handle(mesh.halfedge_handle(var.he_line.back())).idx();
			var.set_start_end(s_id_, e_id_);
		}

		for (int i2 = 0; i2 < b_hes.size(); i2++)
		{
			auto he_h = mesh.halfedge_handle(b_hes[i2]);
			end_vids.insert(mesh.from_vertex_handle(he_h).idx());
		}

		//auto he_h = mesh.halfedge_handle(start_he_id);
		//int startid_ = mesh.from_vertex_handle(he_h).idx();
		//node_direction.emplace_back(startid_, he_h.idx());
	}

	for (int i = 0; i < node_direction.size(); i++)
	{
		if (e_visited.find(node_direction[i].second / 2) != e_visited.end())
			continue;

		vector<int> he_line_;
		e_visited.insert(node_direction[i].second / 2);
		auto he_h = mesh.halfedge_handle(node_direction[i].second);
		int start_vid = mesh.from_vertex_handle(he_h).idx();
		int color_left = f_cs[mesh.face_handle(he_h).idx()];
		int color_right = f_cs[mesh.opposite_face_handle(he_h).idx()];
		he_line_.push_back(he_h.idx());

		auto v_h = mesh.to_vertex_handle(he_h);

		while (end_vids.find(v_h.idx()) == end_vids.end())
		{
			int branch_num = 0;
			auto origin_he = mesh.opposite_halfedge_handle(he_h);
			auto he_out = mesh.next_halfedge_handle(he_h);
			while (he_out.idx() != origin_he.idx())
			{
				int f0 = mesh.face_handle(he_out).idx();
				int f1 = mesh.opposite_face_handle(he_out).idx();
				if (f_cs[f0] != f_cs[f1])
				{
					branch_num++;
					he_h = he_out;
					v_h = mesh.to_vertex_handle(he_h);
				}
				he_out = mesh.next_halfedge_handle(mesh.opposite_halfedge_handle(he_out));
			}
			if (branch_num > 1)
			{
				e_visited.insert(origin_he.idx() / 2);
				v_h = mesh.from_vertex_handle(origin_he);
				he_out = mesh.next_halfedge_handle(mesh.opposite_halfedge_handle(origin_he));
				while (he_out.idx() != origin_he.idx())
				{
					int f0 = mesh.face_handle(he_out).idx();
					int f1 = mesh.opposite_face_handle(he_out).idx();
					if (f_cs[f0] != f_cs[f1])
					{
						node_direction.emplace_back(v_h.idx(), he_out.idx());
					}
					he_out = mesh.next_halfedge_handle(mesh.opposite_halfedge_handle(he_out));
				}
				break;
			}
			else
			{
				he_line_.push_back(he_h.idx());
			}
		}
		if (end_vids.find(v_h.idx()) != end_vids.end())
		{
			e_visited.insert(he_h.idx() / 2);
		}
		abes_p.emplace_back(start_vid, v_h.idx(), color_left, color_right, he_line_);
	}

	cout << "abes size: " << abes_p.size() << endl;

}

void SimplifyMesh::find_first_ring_with_b(const vector<int>& b_he_lines, const int color_split, vector<int>& f_cs, set<int>& ring1_vids)
{
	//color the first ring to -1
	{
		int backid = b_he_lines.back();
		for (int i = 0; i < b_he_lines.size(); i++)
		{
			if (i > 0)
				backid = b_he_lines[i - 1];
			auto he_h_prev = mesh.opposite_halfedge_handle(mesh.halfedge_handle(b_he_lines[i]));
			do
			{
				auto he_h_opp = mesh.opposite_halfedge_handle(he_h_prev);
				he_h_prev = mesh.prev_halfedge_handle(he_h_opp);
				auto f_h = mesh.face_handle(he_h_prev);
				if (f_cs[f_h.idx()] == color_split)
				{
					f_cs[f_h.idx()] = -1;
				}
			} while (he_h_prev.idx() != backid);
		}
	}

	int ring1_start_he_id = -1;
	bool is_found_ring1_ = false;
	//find first he for ring1-loop
	{
		int backid = b_he_lines.back();
		for (int i = 0; i < b_he_lines.size(); i++)
		{
			auto he_h = mesh.halfedge_handle(b_he_lines[i]);
			if (i > 0)
				backid = b_he_lines[i - 1];

			auto he_h_prev = mesh.prev_halfedge_handle(he_h);
			while (he_h_prev.idx() != backid)
			{
				auto he_check = mesh.prev_halfedge_handle(he_h_prev);
				int f0 = mesh.face_handle(he_check).idx();
				int f1 = mesh.opposite_face_handle(he_check).idx();

				if (f0 != -1 && f1 != -1 && f_cs[f0] != f_cs[f1])
				{
					ring1_start_he_id = mesh.opposite_halfedge_handle(he_check).idx();
					is_found_ring1_ = true;
					break;
				}

				auto he_h_opp = mesh.opposite_halfedge_handle(he_h_prev);
				he_h_prev = mesh.prev_halfedge_handle(he_h_opp);
			}
			if (is_found_ring1_)
				break;
		}
	}

	ring1_vids.clear();
	if (is_found_ring1_)
	{
		auto he_h = mesh.halfedge_handle(ring1_start_he_id);
		int start_vid_ = mesh.from_vertex_handle(he_h).idx();
		auto v_h = mesh.to_vertex_handle(he_h);
		ring1_vids.insert(v_h.idx());
		while (v_h.idx() != start_vid_)
		{
			int branch_num = 0;
			auto origin_he = mesh.opposite_halfedge_handle(he_h);
			auto he_out = mesh.next_halfedge_handle(he_h);
			while (he_out.idx() != origin_he.idx())
			{
				int f0 = mesh.face_handle(he_out).idx();
				int f1 = mesh.opposite_face_handle(he_out).idx();
				if (f_cs[f0] != f_cs[f1])
				{
					branch_num++;
					he_h = he_out;
					v_h = mesh.to_vertex_handle(he_h);
				}
				he_out = mesh.next_halfedge_handle(mesh.opposite_halfedge_handle(he_out));
			}
			if (branch_num > 1)
			{
				break;
			}
			else
			{
				ring1_vids.insert(v_h.idx());
			}
		}
	}
	else
	{
		cout << " Can't do this operation" << endl;
		return;
	}

	//color the first ring backto color_split
	{
		int backid = b_he_lines.back();
		for (int i = 0; i < b_he_lines.size(); i++)
		{
			if (i > 0)
				backid = b_he_lines[i - 1];
			auto he_h_prev = mesh.opposite_halfedge_handle(mesh.halfedge_handle(b_he_lines[i]));
			do
			{
				auto he_h_opp = mesh.opposite_halfedge_handle(he_h_prev);
				he_h_prev = mesh.prev_halfedge_handle(he_h_opp);
				auto f_h = mesh.face_handle(he_h_prev);
				if (f_cs[f_h.idx()] == -1)
				{
					f_cs[f_h.idx()] = color_split;
				}
			} while (he_h_prev.idx() != backid);
		}
	}


}

void SimplifyMesh::shortest_path(int start_id, const set<int>& end_set, const vector<int>& f_cs, const vector<int>& c_vec, vector<int>& path, const set<int>& es_ok_)
{

	vector<pair<int, int>> color_pairs;
	if (c_vec.size() == 3)
	{
		color_pairs.emplace_back(c_vec[1], c_vec[1]);
	}
	else if (c_vec.size() == 2)
	{
		for (size_t i = 0; i < c_vec.size(); i++)
		{
			color_pairs.emplace_back(c_vec[i], c_vec[i]);
		}
	}

	for (size_t i = 1; i < c_vec.size(); i++)
	{
		color_pairs.emplace_back(c_vec[i-1], c_vec[i]);
		color_pairs.emplace_back(c_vec[i], c_vec[i-1]);
	}


	auto is_edge_ok = [&](OpenMesh::HalfedgeHandle& he_h)->bool { 
		int e_id_ = he_h.idx() / 2;
		if (es_ok_.size() != 0 &&es_ok_.find(e_id_) == es_ok_.end())
			return false;
		int f0 = mesh.face_handle(he_h).idx();
		int f1 = mesh.opposite_face_handle(he_h).idx();
		if (f0 == -1 || f1 == -1)
			return false;
		for (const auto&var : color_pairs)
		{
			if (f_cs[f0] == var.first&&f_cs[f1] == var.second)
				return true;
		}
		return false;
	};

	std::priority_queue<pair<double, int>, vector<pair<double, int>>, std::greater<pair<double, int>>> Q;
	vector<double> r_dis(V_N, numeric_limits<double>::infinity());
	Q.emplace(0., start_id);
	r_dis[start_id] = 0.;
	vector<int> trace_(V_N, -1);

	int end_p_ = -1;
	while (!Q.empty())
	{
		auto node_ = Q.top();
		Q.pop();
		if (end_set.find(node_.second) != end_set.end())
		{
			end_p_ = node_.second;
			break;
		}
		auto seed_h = mesh.vertex_handle(node_.second);
		for (auto itvv = mesh.voh_begin(seed_h); itvv != mesh.voh_end(seed_h); itvv++)
		{
			if (is_edge_ok(*itvv))
			{
				int vid_ = mesh.to_vertex_handle(*itvv).idx();
				double v_dis_tmp = r_dis[node_.second] + e_len[itvv->idx() / 2];
				if (v_dis_tmp < r_dis[vid_])
				{
					r_dis[vid_] = v_dis_tmp;
					Q.emplace(v_dis_tmp, vid_);
					trace_[vid_] = node_.second;
				}
			}
		}
	}
	vector<int> v_path_;
	while (end_p_ != -1)
	{
		v_path_.push_back(end_p_);
		end_p_ = trace_[end_p_];
	}

	path.clear();
	for (int j = 0; j < v_path_.size(); j++)
	{
		auto f_vh = mesh.vertex_handle(v_path_[j]);
		int n_id_ = v_path_[(j + 1) % v_path_.size()];
		for (auto itvoh = mesh.voh_begin(f_vh); itvoh != mesh.voh_end(f_vh); itvoh++)
		{
			int t_vh_id = mesh.to_vertex_handle(*itvoh).idx();
			if (t_vh_id == n_id_)
			{
				path.push_back(itvoh->idx());
				break;
			}
		}
	}
}

void SimplifyMesh::split_for_tutte(vector<AbstractEdge>& abes_p, vector<int>& f_cs, vector<int>& f_color_num)
{
//	ofstream ofs_es("split_es.txt", ios::trunc);
	int add_fs = 0;
	vector<int> add_fs_color;
	int new_f_num = F_N;
	for (const auto&var : abes_p)
	{
		if (var.he_line.size() < 2)
			continue;
		if (var.color_left == -1 || var.color_right == -1)
			continue;

		set<int> p_b_vs_set;
		p_b_vs_set.insert(var.start_id);
		p_b_vs_set.insert(var.end_id);
		for (int i = 0; i < var.he_line.size(); i++)
		{
			auto he_h = mesh.halfedge_handle(var.he_line[i]);
			int vid_ = mesh.to_vertex_handle(he_h).idx();
			p_b_vs_set.insert(vid_);
		}


		int backid = mesh.to_vertex_handle(mesh.halfedge_handle(var.he_line.front())).idx();
		auto v_h_s_ = mesh.vertex_handle(var.start_id);
		for (auto itvoh = mesh.voh_begin(v_h_s_); itvoh != mesh.voh_end(v_h_s_); itvoh++)
		{
			auto v_t_ = mesh.to_vertex_handle(*itvoh);
			if (v_t_.idx() != backid && p_b_vs_set.find(v_t_.idx()) != p_b_vs_set.end())
			{
				auto e_h_ = mesh.edge_handle(*itvoh);
				auto v0_ = mesh.point(v_h_s_);
				auto v1_ = mesh.point(v_t_);
				auto v_mid = (v0_ + v1_) / 2.0;
				int c_to_ = f_cs[mesh.face_handle(*itvoh).idx()];
				mesh.split(e_h_, v_mid);
//				ofs_es << e_h_ << endl;
				add_fs += 2;
				int fs_tmp = mesh.n_faces() - new_f_num;
				for (size_t j = 0; j < fs_tmp; j++)
				{
					add_fs_color.push_back(c_to_);
					f_color_num[c_to_]++;
				}
				new_f_num = mesh.n_faces();
				break;
			}
		}

		for (int i = 1; i < var.he_line.size(); i++)
		{
			backid = var.he_line[i - 1];
			auto he_h = mesh.halfedge_handle(var.he_line[i]);
			auto he_h_prev = mesh.prev_halfedge_handle(he_h);
			while (he_h_prev.idx() != backid)
			{
				int v_f_ = mesh.from_vertex_handle(he_h_prev).idx();
				auto he_h_opp = mesh.opposite_halfedge_handle(he_h_prev);
				auto e_h_ = mesh.edge_handle(he_h_opp);
				he_h_prev = mesh.prev_halfedge_handle(he_h_opp);

				if (p_b_vs_set.find(v_f_) != p_b_vs_set.end())
				{
					auto v0_ = mesh.point(mesh.from_vertex_handle(he_h_opp));
					auto v1_ = mesh.point(mesh.to_vertex_handle(he_h_opp));
					auto v_mid = (v0_ + v1_) / 2.0;
					mesh.split(e_h_, v_mid);
//					ofs_es << e_h_ << endl;

					add_fs += 2;
					int fs_tmp = mesh.n_faces() - new_f_num;
					for (size_t j = 0; j < fs_tmp; j++)
					{
						add_fs_color.push_back(var.color_left);
						f_color_num[var.color_left]++;
					}
					new_f_num = mesh.n_faces();
				}
			}
		}

		for (int i = var.he_line.size() - 2; i >= 0; i--)
		{
			auto he_h1 = mesh.halfedge_handle(var.he_line[i]);
			auto he_h2 = mesh.halfedge_handle(var.he_line[i + 1]);
			auto he_h1_opp = mesh.opposite_halfedge_handle(he_h1);
			backid = mesh.opposite_halfedge_handle(he_h2).idx();

			auto he_h_prev = mesh.prev_halfedge_handle(he_h1_opp);
			while (he_h_prev.idx() != backid)
			{
				int v_f_ = mesh.from_vertex_handle(he_h_prev).idx();
				auto he_h_opp = mesh.opposite_halfedge_handle(he_h_prev);
				auto e_h_ = mesh.edge_handle(he_h_opp);
				he_h_prev = mesh.prev_halfedge_handle(he_h_opp);

				if (p_b_vs_set.find(v_f_) != p_b_vs_set.end())
				{
					auto v0_ = mesh.point(mesh.from_vertex_handle(he_h_opp));
					auto v1_ = mesh.point(mesh.to_vertex_handle(he_h_opp));
					auto v_mid = (v0_ + v1_) / 2.0;
					mesh.split(e_h_, v_mid);
//					ofs_es << e_h_ << endl;

					add_fs += 2;
					int fs_tmp = mesh.n_faces() - new_f_num;
					for (size_t j = 0; j < fs_tmp; j++)
					{
						add_fs_color.push_back(var.color_right);
						f_color_num[var.color_right]++;
					}
					new_f_num = mesh.n_faces();
				}
			}
		}
	}
//	ofs_es.close();
	int new_F_N = mesh.n_faces();
	f_cs.resize(new_F_N);
	for (int i = F_N; i < new_F_N; i++)
	{
		f_cs[i] = add_fs_color[i - F_N];

	}

	cout << "F_N: " << F_N << " add_fs: " << add_fs << " after split, F_N: " << mesh.n_faces() << endl;


}

void SimplifyMesh::tutte_the_parts(vector<AFace>& Afs, vector<AbstractEdge>& abes_p, set<int>& b_nodes_, vector<int>& f_cs, vector<int>& f_color_num)
{
	//tutte the boundary
	sD.w_uv = Eigen::MatrixXd::Zero(mesh.n_vertices(), 2);

	{
		double factor_ = sqrt(1.0 / M_PI);
		int b_vs_num = sD.boundary_es_del.size();
		for (int i = 0; i < b_vs_num; i++)
		{
			double frac = i * 2. * M_PI / b_vs_num;
			int bvid_ = mesh.from_vertex_handle(mesh.halfedge_handle(sD.boundary_es_del[i])).idx();
			sD.w_uv.row(bvid_) << factor_ * cos(frac), factor_*sin(frac);
		}
	}


	std::map<int, int> vid2idx_map;
	vector<int> idx2vid_vec;

	int corner_count = 0;
	for (auto&var : abes_p)
	{
		if (vid2idx_map.find(var.start_id) == vid2idx_map.end())
		{
			vid2idx_map[var.start_id] = corner_count;
			//var.start_id = corner_count;
			idx2vid_vec.push_back(var.start_id);
			corner_count++;
		}
		if (vid2idx_map.find(var.end_id) == vid2idx_map.end())
		{
			vid2idx_map[var.end_id] = corner_count;
			//var.end_id = corner_count;
			idx2vid_vec.push_back(var.end_id);
			corner_count++;
		}
	}

	//cout << "------------------------------------ " << abes_p.size() << endl;
	//for (int i = 0; i < abes_p.size(); i++)
	//{
	//	auto&var = abes_p[i];
	//	cout << i << " " << var.color_left << " " << var.color_right << " ( " << var.start_id << ", " << var.end_id << " )= " << " ( " << vid2idx_map[var.start_id] << ", " << vid2idx_map[var.end_id] << " ) " << var.he_line.size() << endl;
	//}

	int b_nodes_size = b_nodes_.size();

	Eigen::MatrixXd adjj = Eigen::MatrixXd::Zero(corner_count, corner_count);

	for (auto&var : abes_p)
	{
		int e_sum_ = var.he_line.size();
		adjj(vid2idx_map[var.start_id], vid2idx_map[var.end_id]) = 1.0 / e_sum_;
		adjj(vid2idx_map[var.end_id], vid2idx_map[var.start_id]) = 1.0 / e_sum_;
	}

	Eigen::MatrixXd buv;
	buv.resize(corner_count - b_nodes_size, 2);

	for (size_t i = 0; i < corner_count; i++)
	{
		double i_sum = adjj.row(i).sum();
		for (size_t j = 0; j < corner_count; j++)
		{
			adjj(i, j) = adjj(i, j) / i_sum;
		}
		if (i >= b_nodes_size)
		{
			Eigen::RowVector2d buv_i(0., 0.);
			for (size_t j = 0; j < b_nodes_size; j++)
			{
				buv_i -= adjj(i, j)*sD.w_uv.row(idx2vid_vec[j]);
			}
			buv.row(i - b_nodes_size) = buv_i;
		}
		adjj(i, i) = -1.0;
	}

	cout << "----------------------------------------- " << endl;

	cout << b_nodes_size << "  vs  " << corner_count << endl;

	cout << "----------------------------------------- " << endl;

	Eigen::MatrixXd AA = adjj.block(b_nodes_size, b_nodes_size, corner_count - b_nodes_size, corner_count - b_nodes_size);

	Eigen::FullPivLU<Eigen::MatrixXd>lu_;
	lu_.compute(AA);
	Eigen::MatrixXd node_res;
	node_res.resize(corner_count - b_nodes_size, 2);

	node_res.col(0) = lu_.solve(buv.col(0));
	node_res.col(1) = lu_.solve(buv.col(1));

	cout << "----------------------------------------- " << endl;

	cout << node_res << endl;

	cout << "----------------------------------------- " << endl;


	for (size_t j = 0; j < corner_count - b_nodes_size; j++)
	{
		sD.w_uv.row(idx2vid_vec[j + b_nodes_size]) = node_res.row(j);
	}

	for (const auto&var_ : abes_p)
	{
		if (var_.color_right == -1)
			continue;

		int hes_size = var_.he_line.size();
		Eigen::RowVector2d s_pos_ = sD.w_uv.row(var_.start_id);
		Eigen::RowVector2d e_pos_ = sD.w_uv.row(var_.end_id);
		Eigen::RowVector2d per_ = (e_pos_ - s_pos_) / hes_size;

		for (int j = 1; j < hes_size; j++)
		{
			int vid_ = mesh.from_vertex_handle(mesh.halfedge_handle(var_.he_line[j])).idx();
			sD.w_uv.row(vid_) = s_pos_ + j * per_;
		}
	}

	//ofstream ofs("./data/" + sD.model_name + "adj.txt", ios::trunc);
	//for (size_t i = 0; i < corner_count; i++)
	//{
	//	for (size_t j = 0; j < corner_count; j++)
	//	{
	//		ofs << adjj(i, j) << " ";
	//	}
	//	ofs << endl;
	//}
	//ofs.close();
	//ofs.open("./data/" + sD.model_name + "bxy.txt", ios::trunc);
	//for (const auto&var : b_nodes_)
	//{
	//	//auto vpos_ = mesh.point(mesh.vertex_handle(var));
	//	ofs << vid2idx_map[var] << " " << sD.w_uv(var, 0) << " " << sD.w_uv(var, 1) << endl;
	//}
	//ofs.close();

	for (size_t i = 0; i < Afs.size(); i++)
	{
		if (Afs[i].f_color == -1)
			continue;

		cout << "tutte-the--chart===================: " << i << " color: " << Afs[i].f_color << " fs_num: " << f_color_num[Afs[i].f_color] << endl;
		LOG(INFO) << "tutte-the--chart===================: " << i << " color: " << Afs[i].f_color << " fs_num: " << f_color_num[Afs[i].f_color] << endl;

		set<int> fids_set;
		queue<int> Q;
		for (size_t j = 0; j < Afs[i].ae_ids.size(); j++)
		{
			if (Afs[i].ae_ids[j].first)
			{
				auto &var_ = abes_p[Afs[i].ae_ids[j].second].he_line;
				for (size_t j = 0; j < var_.size(); j++)
				{
					auto he_h = mesh.halfedge_handle(var_[j]);
					int f1 = mesh.face_handle(he_h).idx();
					if (fids_set.find(f1) == fids_set.end() && f_cs[f1] == Afs[i].f_color)
					{
						Q.push(f1);
						fids_set.insert(f1);
					}
				}
			}
			else
			{
				auto &var_ = abes_p[Afs[i].ae_ids[j].second].he_line;
				for (size_t j = 0; j < var_.size(); j++)
				{
					auto he_h = mesh.halfedge_handle(var_[j]);
					int f1 = mesh.opposite_face_handle(he_h).idx();
					if (fids_set.find(f1) == fids_set.end() && f_cs[f1] == Afs[i].f_color)
					{
						Q.push(f1);
						fids_set.insert(f1);
					}
				}
			}
		}
		while (!Q.empty())
		{
			auto node_ = Q.front();
			Q.pop();
			auto seed_f_h = mesh.face_handle(node_);
			for (auto itvf = mesh.ff_begin(seed_f_h); itvf != mesh.ff_end(seed_f_h); itvf++)
			{
				if (fids_set.find(itvf->idx()) == fids_set.end() && f_cs[itvf->idx()] == Afs[i].f_color)//== delete_color
				{
					Q.push(itvf->idx());
					fids_set.insert(itvf->idx());
				}
			}
		}
		if (fids_set.size() != f_color_num[Afs[i].f_color])
		{
			cout << "Wrong!!! " << fids_set.size() << " vs " << f_color_num[Afs[i].f_color] << endl;
		}

	//	sD.patches_fids.push_back(fids_set);

		tutte_for_one_part(Afs[i], abes_p, fids_set, f_cs);
	}

}

void SimplifyMesh::tutte_for_one_part(AFace & af, vector<AbstractEdge>& abes_p, set<int>& fids_set, vector<int>& f_cs)
{

	int mv_num = mesh.n_vertices();

	Eigen::MatrixXi PFsrc;
	PFsrc.resize(fids_set.size(), 3);

	int c_ = 0;
	for (auto&var : fids_set)
	{
		//PFsrc.row(c_) = F.row(var);
		//c_++;

		auto f_h_ = mesh.face_handle(var);
		auto itfv = mesh.fv_begin(f_h_);
		PFsrc(c_, 0) = itfv->idx();
		itfv++;
		PFsrc(c_, 1) = itfv->idx();
		itfv++;
		PFsrc(c_, 2) = itfv->idx();
		c_++;

	}
	Eigen::VectorXi v_src2s;// size: mv_num, others -1;
	Eigen::VectorXi v_s2src;// size: fix_V.rows();

	Eigen::MatrixXi p_F;
	remove_unreferenced(mv_num, PFsrc, v_s2src, p_F, v_src2s);

	int p_b_num = 0;
	for (size_t i = 0; i < af.ae_ids.size(); i++)
	{
		p_b_num += abes_p[af.ae_ids[i].second].he_line.size();
	}
	Eigen::MatrixXd uv_init;
	Eigen::MatrixXd bnd_uv;
	Eigen::VectorXi bnd;
	bnd.resize(p_b_num);
	bnd_uv.resize(p_b_num, 2);
	set<int> p_b_vs_set;
	c_ = 0;
	for (size_t i = 0; i < af.ae_ids.size(); i++)
	{
		if (af.ae_ids[i].first)
		{
			auto &var_ = abes_p[af.ae_ids[i].second].he_line;
			for (size_t j = 0; j < var_.size(); j++)
			{
				int v_src_id_ = mesh.from_vertex_handle(mesh.halfedge_handle(var_[j])).idx();
				bnd[c_] = v_src2s[v_src_id_];
				p_b_vs_set.insert(bnd[c_]);
				bnd_uv.row(c_) = sD.w_uv.row(v_src_id_);
				c_++;
			}
		}
		else
		{
			auto &var_ = abes_p[af.ae_ids[i].second].he_line;
			for (size_t j = 0; j < var_.size(); j++)
			{
				int v_src_id_ = mesh.to_vertex_handle(mesh.halfedge_handle(var_[j])).idx();
				bnd[c_] = v_src2s[v_src_id_];
				p_b_vs_set.insert(bnd[c_]);
				bnd_uv.row(c_) = sD.w_uv.row(v_src_id_);
				c_++;
			}
		}
	}


	Tutte(v_s2src.size(), p_F, bnd, bnd_uv, uv_init);

	for (size_t i = 0; i < v_s2src.size(); i++)
	{
		if (p_b_vs_set.find(i) == p_b_vs_set.end())
			sD.w_uv.row(v_s2src[i]) = uv_init.row(i);
	}


}

void SimplifyMesh::curve_to_line(vector<AFace>& Afs, vector<AbstractEdge>& abes_p, vector<int>& f_cs, vector<int>& f_color_num)
{
	for (int ae_idx = 0; ae_idx < abes_p.size(); ae_idx++)
	{
		auto&ae_ = abes_p[ae_idx];
		int c1 = ae_.color_left;
		int c2 = ae_.color_right;

		if (c1 == -1 || c2 == -1)
			continue;
		if (f_color_num[c1] == 0 || f_color_num[c2] == 0)
			continue;

		int start_id = ae_.start_id;
		int end_id = ae_.end_id;
		//dijkstra find s_path_

		vector<int> s_path_;

		{
			set<int> hes_; 
			int loops_num = 5;
			queue<pair<int, int>>Q;
			set<int> vs_;
			//add half_edges
			for (const auto&var_ : ae_.he_line)
			{
				auto he_h_ = mesh.halfedge_handle(var_);
				hes_.insert(var_ / 2);
				int v_id_ = mesh.from_vertex_handle(he_h_).idx();
				vs_.insert(v_id_);
				Q.emplace(v_id_, 0);
			}
			Q.emplace(ae_.end_id, 0);
			vs_.insert(ae_.end_id);

			while (!Q.empty())
			{
				auto node_ = Q.front();
				Q.pop();
				if (node_.second >= loops_num)
					continue;
				auto v_h_ = mesh.vertex_handle(node_.first);
				for (auto itvoh = mesh.voh_begin(v_h_); itvoh != mesh.voh_end(v_h_); itvoh++)
				{
					hes_.insert(itvoh->idx() / 2);
					int v_t_id_ = mesh.to_vertex_handle(*itvoh).idx();
					if (vs_.find(v_t_id_) == vs_.end())
					{
						Q.emplace(v_t_id_, node_.second + 1);
						vs_.insert(v_t_id_);
					}
				}
			}

			//delete half_edges
			for (int i = 0; i < Afs.size(); i++)
			{
				if (Afs[i].f_color == c1|| Afs[i].f_color == c2)
				{
					for (const auto&var : Afs[i].ae_ids)
					{
						if (var.second == ae_idx)
							continue;
						if (var.first)
						{
							for (const auto&heid_ : abes_p[var.second].he_line)
							{
								auto he_h = mesh.halfedge_handle(heid_);
								hes_.erase(heid_ / 2);
								auto v_h_ = mesh.from_vertex_handle(he_h);
								if (v_h_.idx() != start_id && v_h_.idx() != end_id)
								{
									for (auto itvoh = mesh.voh_begin(v_h_); itvoh != mesh.voh_end(v_h_); itvoh++)
									{
										hes_.erase(itvoh->idx() / 2);
									}
								}
							}
						}
						else
						{
							for (const auto&heid_ : abes_p[var.second].he_line)
							{
								auto he_h = mesh.halfedge_handle(heid_);
								hes_.erase(heid_ / 2);
								auto v_h_ = mesh.to_vertex_handle(he_h);
								if (v_h_.idx() != start_id && v_h_.idx() != end_id)
								{
									for (auto itvoh = mesh.voh_begin(v_h_); itvoh != mesh.voh_end(v_h_); itvoh++)
									{
										hes_.erase(itvoh->idx() / 2);
									}
								}

							}
						}

					}

				}
			}

			shortest_path(end_id, set<int>{start_id}, f_cs, vector<int>{c1, c2}, s_path_,hes_);
		}

		ae_.he_line = s_path_;

		set<int> block_fids_1, block_fids_2;

		for (int i = 0; i < Afs.size(); i++)
		{
			if (Afs[i].f_color == c1)
			{
				for (const auto&var : Afs[i].ae_ids)
				{
					if (var.second == ae_idx)
						continue;
					if (var.first)
					{
						for (const auto&heid_ : abes_p[var.second].he_line)
						{
							auto he_h = mesh.halfedge_handle(heid_);
							int f_id_ = mesh.opposite_face_handle(he_h).idx();
							block_fids_1.insert(f_id_);
						}
					}
					else
					{
						for (const auto&heid_ : abes_p[var.second].he_line)
						{
							auto he_h = mesh.halfedge_handle(heid_);
							int f_id_ = mesh.face_handle(he_h).idx();
							block_fids_1.insert(f_id_);
						}
					}
				}
				break;
			}
		}

		for (int i = 0; i < Afs.size(); i++)
		{
			if (Afs[i].f_color == c2)
			{
				for (const auto&var : Afs[i].ae_ids)
				{
					if (var.second == ae_idx)
						continue;
					if (var.first)
					{
						for (const auto&heid_ : abes_p[var.second].he_line)
						{
							auto he_h = mesh.halfedge_handle(heid_);
							int f_id_ = mesh.opposite_face_handle(he_h).idx();
							block_fids_2.insert(f_id_);
						}
					}
					else
					{
						for (const auto&heid_ : abes_p[var.second].he_line)
						{
							auto he_h = mesh.halfedge_handle(heid_);
							int f_id_ = mesh.face_handle(he_h).idx();
							block_fids_2.insert(f_id_);
						}
					}
				}
				break;
			}
		}

		//change f_cs
		queue<int> Q1, Q2;
		int c1_num = 0;
		int c2_num = 0;
		for (int j = 0; j < s_path_.size(); j++)
		{
			auto he_h = mesh.halfedge_handle(s_path_[j]);
			int f1 = mesh.face_handle(he_h).idx();
			int f2 = mesh.opposite_face_handle(he_h).idx();
			block_fids_2.insert(f1);
			block_fids_1.insert(f2);

			if (-1 != f1 && f_cs[f1] == c2)//== delete_color
			{
				if (f_cs[f1] == c2)
					c2_num--;
				f_cs[f1] = c1;
				Q1.push(f1);
				c1_num++;
			}
			if (-1 != f2 && f_cs[f2] == c1)//== delete_color
			{
				if (f_cs[f2] == c1)
					c1_num--;
				f_cs[f2] = c2;
				Q2.push(f2);
				c2_num++;
			}
		}

		while (!Q1.empty())
		{
			auto node_ = Q1.front();
			Q1.pop();
			auto seed_f_h = mesh.face_handle(node_);
			for (auto itvf = mesh.ff_begin(seed_f_h); itvf != mesh.ff_end(seed_f_h); itvf++)
			{
				if (block_fids_1.find(itvf->idx()) == block_fids_1.end() && f_cs[itvf->idx()] == c2)//== delete_color
				{
					if (f_cs[itvf->idx()] == c2)
						c2_num--;
					f_cs[itvf->idx()] = c1;
					c1_num++;
					Q1.push(itvf->idx());
				}
			}
		}

		while (!Q2.empty())
		{
			auto node_ = Q2.front();
			Q2.pop();
			auto seed_f_h = mesh.face_handle(node_);
			for (auto itvf = mesh.ff_begin(seed_f_h); itvf != mesh.ff_end(seed_f_h); itvf++)
			{
				if (block_fids_2.find(itvf->idx()) == block_fids_2.end() && f_cs[itvf->idx()] == c1)//== delete_color
				{
					if (f_cs[itvf->idx()] == c1)
						c1_num--;
					f_cs[itvf->idx()] = c2;
					c2_num++;
					Q2.push(itvf->idx());
				}
			}
		}

		if ((c1_num + c2_num) != 0)
		{
			cout << "s_path_-size: " << s_path_.size() << " start_id: " << start_id << " end_id: " << end_id << endl;
			cout << "Handle_encircle_case Wrong!!! c1: " << c1_num << " c2: " << c2_num << endl;
		}
		if (-1 != c1)
			f_color_num[c1] += c1_num;
		if (-1 != c2)
			f_color_num[c2] += c2_num;


	}

}

