#pragma once
#include"Common\Common.h"
#include"SharedData.h"

struct AbstractEdge
{
	AbstractEdge(int color_left_,
		int color_right_,
		const vector<int>& he_line_)
	{
		color_left = color_left_;
		color_right = color_right_;
		he_line = he_line_;
	}

	AbstractEdge(int start_id_,
		int end_id_,
		int color_left_,
		int color_right_,
		const vector<int>& he_line_)
	{
		start_id = start_id_;
		end_id = end_id_;
		color_left = color_left_;
		color_right = color_right_;
		he_line = he_line_;
	}
	int start_id;
	int end_id;
	int color_left;
	int color_right;
	vector<int> he_line;

	void set_start_end(int s_id_, int e_id_)
	{
		start_id = s_id_;
		end_id = e_id_;
	}
};

struct AFace
{
	AFace(int f_c_, const vector<pair<bool, int>>& ae_ids_)
	{
		f_color = f_c_;
		ae_ids = ae_ids_;
	}
	vector<pair<bool, int>> ae_ids;
	int f_color;
};

class SimplifyMesh
{
public:
	SimplifyMesh(SharedData& sharedd_);
	~SimplifyMesh();

public:
	void prepare_for_dijkstra();
	void color_src_mesh();


	//////////////////////////

	bool handle_encircle_case(vector<AFace>& Afs, vector<AbstractEdge>& abes_p, vector<int>& f_cs, int color_start_, vector<int>& f_color_num, int f_part_hope);

	void find_inter_line(vector<AbstractEdge>& abes_p, set<int>& b_nodes, const vector<int>& f_cs, const vector<int>& b_hes);

	void find_first_ring_with_b(const vector<int>&b_line, const int color_split, vector<int>& f_cs, set<int>&ring1_vids);
	void shortest_path(int start_id, const set<int>& end_set, const vector<int>& f_cs, const vector<int>& c_vec, vector<int>& path, const set<int>& es_ok_ = set<int>{});
	void split_for_tutte(vector<AbstractEdge>& abes_p, vector<int>& f_cs, vector<int>& f_color_num);
	void color_disk(const vector<int>& b_he_lines, vector<int>& f_cs, int &color_start, int color_split, int f_total_num, int f_part_hope, pair<int, int> b_s_e_idx = pair<int, int>(0, 0));
	void curve_to_line(vector<AFace>& Afs, vector<AbstractEdge>& abes_p, vector<int>& f_cs, vector<int>& f_color_num);


	void tutte_the_parts(vector<AFace>& Afs, vector<AbstractEdge>& abes_p, set<int>& b_nodes_, vector<int>& f_cs, vector<int>& f_color_num);

	void tutte_for_one_part(AFace& af, vector<AbstractEdge>& abes_p, set<int>& fids_, vector<int>&f_cs);


public:
	SharedData& sD;
	Mesh& mesh;
	int V_N;
	int F_N;


	vector<double> e_len;



};

