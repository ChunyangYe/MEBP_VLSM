#include "BernsteinBase.h"
Triplet3d operator+(const Triplet3d & a, const Triplet3d & b)
{
	return std::make_tuple(std::get<0>(a)+std::get<0>(b), std::get<1>(a) + std::get<1>(b), std::get<2>(a) + std::get<2>(b));
}
Triplet3d operator-(const Triplet3d & a, const Triplet3d & b)
{
	return std::make_tuple(std::get<0>(a) - std::get<0>(b), std::get<1>(a) - std::get<1>(b), std::get<2>(a) - std::get<2>(b));
}
Triplet3d operator*(const double & a, const Triplet3d & b)
{
	return std::make_tuple(a*std::get<0>(b), a*std::get<1>(b), a*std::get<2>(b));
}
Triplet3d operator/(const Triplet3d & a, const double & b)
{
	return std::make_tuple(std::get<0>(a) / b, std::get<1>(a) / b, std::get<2>(a) / b);
}
Triplet3d bezier_300(const double & kesi0, const double & kesi1)
{
	return	std::make_tuple(kesi0*kesi0*kesi0, 3.0*kesi0*kesi0, 0.0);
	//	return kesi0*kesi0*kesi0;
}

Triplet3d bezier_210(const double & kesi0, const double & kesi1)
{
	return	std::make_tuple(3.0*kesi0*kesi0*kesi1, 6.0*kesi0*kesi1, 3.0*kesi0*kesi0);
	//return 3.0*kesi0*kesi0*kesi1;
}

Triplet3d bezier_201(const double & kesi0, const double & kesi1)
{
	return std::make_tuple(3.0*kesi0*kesi0*(1.0 - kesi0 - kesi1), 3.0*kesi0*(2. - 3.*kesi0 - 2.*kesi1), -3.*kesi0*kesi0);
}

Triplet3d bezier_120(const double & kesi0, const double & kesi1)
{
	return std::make_tuple(3.0*kesi0*kesi1*kesi1, 3.*kesi1*kesi1, 6.*kesi0*kesi1);
}

Triplet3d bezier_111(const double & kesi0, const double & kesi1)
{
	return std::make_tuple(6.0*kesi0*kesi1*(1.0-kesi0-kesi1),6.*kesi1*(1.-2.*kesi0-kesi1),6.*kesi0*(1.-kesi0-2.*kesi1));
}

Triplet3d bezier_102(const double & kesi0, const double & kesi1)
{
	double kesi2 = 1.0 - kesi0 - kesi1;
	return std::make_tuple(3.0*kesi0*kesi2*kesi2,3.*kesi2*kesi2-6.*kesi0*kesi2,-6.*kesi0*kesi2);
}

Triplet3d bezier_030(const double & kesi0, const double & kesi1)
{
	return std::make_tuple(kesi1*kesi1*kesi1,0.,3.*kesi1*kesi1);
}

Triplet3d bezier_021(const double & kesi0, const double & kesi1)
{
	return std::make_tuple(3.0*kesi1*kesi1*(1.0-kesi0-kesi1),-3.*kesi1*kesi1,3.*kesi1*(2.-2.*kesi0-3.*kesi1));
}

Triplet3d bezier_012(const double & kesi0, const double & kesi1)
{
	double kesi2 = 1.0 - kesi0 - kesi1;
	return std::make_tuple(3.0*kesi1*kesi2*kesi2,-6.*kesi1*kesi2,3.*kesi2*(kesi2-2.*kesi1));
}

Triplet3d bezier_003(const double & kesi0, const double & kesi1)
{
	double kesi2 = 1.0 - kesi0 - kesi1;
	return std::make_tuple(kesi2*kesi2*kesi2,-3.*kesi2*kesi2,-3.*kesi2*kesi2);
}

Triplet3d bezier_200(const double & kesi0, const double & kesi1)
{
	return std::make_tuple(kesi0*kesi0,2.*kesi0,0.);
}

Triplet3d bezier_110(const double & kesi0, const double & kesi1)
{
	return std::make_tuple(2.0*kesi0*kesi1,2.*kesi1,2.*kesi0);
}

Triplet3d bezier_101(const double & kesi0, const double & kesi1)
{
	return std::make_tuple(2.0*kesi0*(1.0-kesi0-kesi1),2.*(1.-2.*kesi0-kesi1),-2.*kesi0);
}

Triplet3d bezier_020(const double & kesi0, const double & kesi1)
{
	return std::make_tuple(kesi1*kesi1, 0., 2.*kesi1);
}

Triplet3d bezier_011(const double & kesi0, const double & kesi1)
{
	return std::make_tuple(2.0*kesi1*(1.0 - kesi0 - kesi1), -2.*kesi1, 2.*(1. - kesi0 - 2.*kesi1));
}

Triplet3d bezier_002(const double & kesi0, const double & kesi1)
{
	double kesi2 = 1.0 - kesi0 - kesi1;
	return std::make_tuple(kesi2*kesi2,-2.*kesi2,-2.*kesi2);
}

Triplet3d bezier_100(const double & kesi0, const double & kesi1)
{
	return std::make_tuple(kesi0, 1., 0.);
}

Triplet3d bezier_010(const double & kesi0, const double & kesi1)
{
	return std::make_tuple(kesi1, 0., 1.0);
}

Triplet3d bezier_001(const double & kesi0, const double & kesi1)
{
	return std::make_tuple(1.0-kesi0-kesi1,-1.,-1.);
}

Triplet3d bezier_3_302(const double & kesi0, const double & kesi1)
{
	return bezier_300(kesi0, kesi1) + bezier_210(kesi0, kesi1)/3.0;
}

Triplet3d bezier_3_301(const double & kesi0, const double & kesi1)
{
	return bezier_300(kesi0, kesi1) + (2.0 / 3.0)*bezier_210(kesi0, kesi1) + bezier_120(kesi0, kesi1) / 3.0;
}

Triplet3d bezier_3_203(const double & kesi0, const double & kesi1)
{
	return bezier_300(kesi0, kesi1) + bezier_201(kesi0, kesi1)/3.0;
}

Triplet3d bezier_3_202(const double & kesi0, const double & kesi1)
{
	return bezier_300(kesi0, kesi1) + (bezier_210(kesi0, kesi1)+ bezier_201(kesi0, kesi1))/3.0;
}

Triplet3d bezier_3_201(const double & kesi0, const double & kesi1)
{
	return bezier_300(kesi0, kesi1) + bezier_201(kesi0, kesi1)/3.0 + (2.0 / 3.0)*bezier_210(kesi0, kesi1) + bezier_120(kesi0, kesi1) / 3.0;
}

Triplet3d bezier_3_103(const double & kesi0, const double & kesi1)
{
	return bezier_300(kesi0, kesi1) + (2.0 / 3.0)*bezier_201(kesi0, kesi1) + bezier_102(kesi0, kesi1) / 3.0;
}

Triplet3d bezier_3_102(const double & kesi0, const double & kesi1)
{
	return bezier_300(kesi0, kesi1) + (2.0 / 3.0)*bezier_201(kesi0, kesi1) + bezier_102(kesi0, kesi1) / 3.0 + bezier_210(kesi0, kesi1)/ 3.0;
}

Triplet3d bezier_3_101(const double & kesi0, const double & kesi1)
{
	return bezier_300(kesi0, kesi1) + (2.0 / 3.0)*(bezier_210(kesi0, kesi1) + bezier_201(kesi0, kesi1)) + (bezier_102(kesi0, kesi1) + bezier_120(kesi0, kesi1)) / 3.0;
}



Triplet3d bezier_2_201(const double & kesi0, const double & kesi1)
{
	return bezier_200(kesi0, kesi1) + 0.5*bezier_110(kesi0, kesi1);
}

Triplet3d bezier_2_102(const double & kesi0, const double & kesi1)
{
	return bezier_200(kesi0, kesi1) + 0.5*bezier_101(kesi0, kesi1);
}

Triplet3d bezier_2_101(const double & kesi0, const double & kesi1)
{
	return bezier_200(kesi0, kesi1) + 0.5*(bezier_110(kesi0, kesi1) + bezier_101(kesi0, kesi1));
}

// for 1
Triplet3d bezier_3_312(const double & kesi0, const double & kesi1)
{
	return bezier_030(kesi0, kesi1) + bezier_021(kesi0, kesi1)/3.0;
}

Triplet3d bezier_3_311(const double & kesi0, const double & kesi1)
{
	return bezier_030(kesi0, kesi1) + 2.0 / 3.0*bezier_021(kesi0, kesi1) + bezier_012(kesi0, kesi1) / 3.0;
}

Triplet3d bezier_3_213(const double & kesi0, const double & kesi1)
{
	return bezier_030(kesi0, kesi1) + bezier_120(kesi0, kesi1)/3.0;
}

Triplet3d bezier_3_212(const double & kesi0, const double & kesi1)
{
	return bezier_030(kesi0, kesi1) + (bezier_021(kesi0, kesi1) + bezier_120(kesi0, kesi1))/3.0;
}

Triplet3d bezier_3_211(const double & kesi0, const double & kesi1)
{
	return bezier_030(kesi0, kesi1) + bezier_120(kesi0, kesi1)/3.0 + 2.0 / 3.0*bezier_021(kesi0, kesi1) + bezier_012(kesi0, kesi1) / 3.0;
}

Triplet3d bezier_3_113(const double & kesi0, const double & kesi1)
{
	return bezier_030(kesi0, kesi1) + 2.0 / 3.0*bezier_120(kesi0, kesi1) + bezier_210(kesi0, kesi1) / 3.0;
}

Triplet3d bezier_3_112(const double & kesi0, const double & kesi1)
{
	return bezier_030(kesi0, kesi1) + 2.0 / 3.0*bezier_120(kesi0, kesi1) + bezier_210(kesi0, kesi1) / 3.0 + bezier_021(kesi0, kesi1)/3.0;
}

Triplet3d bezier_3_111(const double & kesi0, const double & kesi1)
{
	return bezier_030(kesi0, kesi1) + 2.0 / 3.0*(bezier_120(kesi0, kesi1) + bezier_021(kesi0, kesi1)) + (bezier_210(kesi0, kesi1) + bezier_012(kesi0, kesi1)) / 3.0;
}



Triplet3d bezier_2_211(const double & kesi0, const double & kesi1)
{
	return bezier_020(kesi0, kesi1) + 0.5*bezier_011(kesi0, kesi1);
}

Triplet3d bezier_2_112(const double & kesi0, const double & kesi1)
{
	return bezier_020(kesi0, kesi1) + 0.5*bezier_110(kesi0, kesi1);
}

Triplet3d bezier_2_111(const double & kesi0, const double & kesi1)
{
	return bezier_020(kesi0, kesi1) + 0.5*(bezier_011(kesi0, kesi1) + bezier_110(kesi0, kesi1));
}

//for 2
Triplet3d bezier_3_322(const double & kesi0, const double & kesi1)
{
	return bezier_003(kesi0, kesi1) + bezier_102(kesi0, kesi1)/3.0;
}

Triplet3d bezier_3_321(const double & kesi0, const double & kesi1)
{
	return bezier_003(kesi0, kesi1) + 2.0 / 3.0*bezier_102(kesi0, kesi1) + bezier_201(kesi0, kesi1) / 3.0;
}

Triplet3d bezier_3_223(const double & kesi0, const double & kesi1)
{
	return bezier_003(kesi0, kesi1) + bezier_012(kesi0, kesi1)/3.0;
}

Triplet3d bezier_3_222(const double & kesi0, const double & kesi1)
{
	return bezier_003(kesi0, kesi1) + (bezier_102(kesi0, kesi1) + bezier_012(kesi0, kesi1))/3.0;
}

Triplet3d bezier_3_221(const double & kesi0, const double & kesi1)
{
	return bezier_003(kesi0, kesi1) + bezier_012(kesi0, kesi1)/3.0 + 2.0 / 3.0*bezier_102(kesi0, kesi1) + bezier_201(kesi0, kesi1) / 3.0;
}

Triplet3d bezier_3_123(const double & kesi0, const double & kesi1)
{
	return bezier_003(kesi0, kesi1) + 2.0 / 3.0*bezier_012(kesi0, kesi1) + bezier_021(kesi0, kesi1) / 3.0;
}

Triplet3d bezier_3_122(const double & kesi0, const double & kesi1)
{
	return bezier_003(kesi0, kesi1) + 2.0 / 3.0*bezier_012(kesi0, kesi1) + bezier_021(kesi0, kesi1) / 3.0 + bezier_102(kesi0, kesi1)/3.0;
}

Triplet3d bezier_3_121(const double & kesi0, const double & kesi1)
{
	return bezier_003(kesi0, kesi1) + 2.0 / 3.0*(bezier_012(kesi0, kesi1) + bezier_102(kesi0, kesi1)) + (bezier_021(kesi0, kesi1) + bezier_201(kesi0, kesi1)) / 3.0;
}



Triplet3d bezier_2_221(const double & kesi0, const double & kesi1)
{
	return bezier_002(kesi0, kesi1) + 0.5*bezier_101(kesi0, kesi1);
}

Triplet3d bezier_2_122(const double & kesi0, const double & kesi1)
{
	return bezier_002(kesi0, kesi1) + 0.5*bezier_011(kesi0, kesi1);
}

Triplet3d bezier_2_121(const double & kesi0, const double & kesi1)
{
	return bezier_002(kesi0, kesi1) + 0.5*(bezier_101(kesi0, kesi1) + bezier_011(kesi0, kesi1));
}


//
Triplet3d bezier_3_2_e01(const double & kesi0, const double & kesi1)
{
	return (2.0 / 3.0)*(bezier_210(kesi0, kesi1) + bezier_120(kesi0, kesi1));
}

Triplet3d bezier_3_2_e12(const double & kesi0, const double & kesi1)
{
	return (2.0 / 3.0)*(bezier_021(kesi0, kesi1) + bezier_012(kesi0, kesi1));
}

Triplet3d bezier_3_2_e20(const double & kesi0, const double & kesi1)
{
	return (2.0 / 3.0)*(bezier_102(kesi0, kesi1) + bezier_201(kesi0, kesi1));
}



