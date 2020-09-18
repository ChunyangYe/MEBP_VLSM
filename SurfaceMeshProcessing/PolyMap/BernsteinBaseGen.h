#pragma once
#include<tuple>
#include<vector>
struct FunInfo
{
	static enum BASETYPE
	{
		V,E,F
	};
	BASETYPE basetype;//0:v,1:e,2:f
	int n;
	int i;
	int j;
	int vl_od;
	int vr_od;
	int e_id;//012
	int v_id;//012
	int e_op_od;
	static std::vector<std::vector<std::vector<double>>> coef_v;
	static std::vector<std::vector<std::vector<double>>> coef_e;

};
std::vector<std::vector<std::vector<double>>> FunInfo::coef_v = {{},
																{{0.5}},
																{{2.0/3.0,1.0/3.0},{1.0/3.0}},
																{{},{},{}} };

std::vector<std::vector<std::vector<double>>> FunInfo::coef_e = {{},
																 {},
																{{},{2.0 / 3.0,2.0/3.0}},
																{{},{},{}} };

typedef std::tuple<double, double, double> Triplet3d;

Triplet3d operator+(const Triplet3d& a, const Triplet3d& b);
Triplet3d operator-(const Triplet3d& a, const Triplet3d& b);
Triplet3d operator*(const double& a, const Triplet3d& b);
Triplet3d operator/(const Triplet3d& a, const double& b);



double factorial(int n);
Triplet3d bezier_base(const int& n, const int& i, const int& j, const double& kesi0, const double& kesi1);

Triplet3d bezier_base_gen(const FunInfo& funinfo, const double & kesi0, const double & kesi1);

Triplet3d bezier_300(const double& kesi0, const double& kesi1);

Triplet3d bezier_210(const double& kesi0, const double& kesi1);

Triplet3d bezier_201(const double& kesi0, const double& kesi1);

Triplet3d bezier_120(const double& kesi0, const double& kesi1);

Triplet3d bezier_111(const double& kesi0, const double& kesi1);

Triplet3d bezier_102(const double& kesi0, const double& kesi1);

Triplet3d bezier_030(const double& kesi0, const double& kesi1);

Triplet3d bezier_021(const double& kesi0, const double& kesi1);

Triplet3d bezier_012(const double& kesi0, const double& kesi1);

Triplet3d bezier_003(const double& kesi0, const double& kesi1);

Triplet3d bezier_200(const double& kesi0, const double& kesi1);

Triplet3d bezier_110(const double& kesi0, const double& kesi1);

Triplet3d bezier_101(const double& kesi0, const double& kesi1);

Triplet3d bezier_020(const double& kesi0, const double& kesi1);

Triplet3d bezier_011(const double& kesi0, const double& kesi1);

Triplet3d bezier_002(const double& kesi0, const double& kesi1);

Triplet3d bezier_100(const double& kesi0, const double& kesi1);
Triplet3d bezier_010(const double& kesi0, const double& kesi1);
Triplet3d bezier_001(const double& kesi0, const double& kesi1);


//vertex 0, left order-3, right order-2
Triplet3d bezier_3_302(const double& kesi0, const double& kesi1);
Triplet3d bezier_3_301(const double& kesi0, const double& kesi1);
Triplet3d bezier_3_203(const double& kesi0, const double& kesi1);
Triplet3d bezier_3_202(const double& kesi0, const double& kesi1);
Triplet3d bezier_3_201(const double& kesi0, const double& kesi1);
Triplet3d bezier_3_103(const double& kesi0, const double& kesi1);
Triplet3d bezier_3_102(const double& kesi0, const double& kesi1);
Triplet3d bezier_3_101(const double& kesi0, const double& kesi1);

Triplet3d bezier_2_201(const double& kesi0, const double& kesi1);
Triplet3d bezier_2_102(const double& kesi0, const double& kesi1);
Triplet3d bezier_2_101(const double& kesi0, const double& kesi1);

//
Triplet3d bezier_3_312(const double& kesi0, const double& kesi1);
Triplet3d bezier_3_311(const double& kesi0, const double& kesi1);
Triplet3d bezier_3_213(const double& kesi0, const double& kesi1);
Triplet3d bezier_3_212(const double& kesi0, const double& kesi1);
Triplet3d bezier_3_211(const double& kesi0, const double& kesi1);
Triplet3d bezier_3_113(const double& kesi0, const double& kesi1);
Triplet3d bezier_3_112(const double& kesi0, const double& kesi1);
Triplet3d bezier_3_111(const double& kesi0, const double& kesi1);

Triplet3d bezier_2_211(const double& kesi0, const double& kesi1);
Triplet3d bezier_2_112(const double& kesi0, const double& kesi1);
Triplet3d bezier_2_111(const double& kesi0, const double& kesi1);

//
Triplet3d bezier_3_322(const double& kesi0, const double& kesi1);
Triplet3d bezier_3_321(const double& kesi0, const double& kesi1);
Triplet3d bezier_3_223(const double& kesi0, const double& kesi1);
Triplet3d bezier_3_222(const double& kesi0, const double& kesi1);
Triplet3d bezier_3_221(const double& kesi0, const double& kesi1);
Triplet3d bezier_3_123(const double& kesi0, const double& kesi1);
Triplet3d bezier_3_122(const double& kesi0, const double& kesi1);
Triplet3d bezier_3_121(const double& kesi0, const double& kesi1);

Triplet3d bezier_2_221(const double& kesi0, const double& kesi1);
Triplet3d bezier_2_122(const double& kesi0, const double& kesi1);
Triplet3d bezier_2_121(const double& kesi0, const double& kesi1);

//
Triplet3d bezier_3_2_e01(const double& kesi0, const double& kesi1);
Triplet3d bezier_3_2_e12(const double& kesi0, const double& kesi1);
Triplet3d bezier_3_2_e20(const double& kesi0, const double& kesi1);


typedef Triplet3d (*BezierBase)(const double&, const double&);


