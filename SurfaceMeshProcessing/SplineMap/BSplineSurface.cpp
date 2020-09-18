#include "BSplineSurface.h"
#include <ctime>


BSplineSurface::BSplineSurface()
{
	//TestSet();
}


BSplineSurface::~BSplineSurface()
{
}

Eigen::Vector3d BSplineSurface::operator()(const double & u, const double & v) const
{
	return DeBoor(ctrlpoints, u, v);
}

Eigen::Vector3d BSplineSurface::getPoint(const double &u, const double &v)
{
//	return DeBoor(ctrlpoints, u, v);

	size_t m = ctrlpoints.size();
	int rx = 0;
	for (size_t i = degreex; i + 1 <= ctrlpoints.size(); i++)
	{
		if (u >= knotsx[i])
		{
			rx = static_cast<int>(i);
		}
		if (u < knotsx[i + 1])
		{
			break;
		}
	}
	double t = u - knotsx[rx];
	std::vector<double> Baseu(4);
	if (rx < 4)
	{
		Baseu = FirstInterval(t);
	}
	else if (rx < 5)
	{
		Baseu = SecondInterval(t);
	}
	else if (rx < (m - 2))
	{
		Baseu = CommonInterval(t);
	}
	else if (rx < (m - 1))
	{
		Baseu = RSecondInterval(t);
	}
	else
	{
		Baseu = RFirstInterval(t);
	}

	int ry = 0;
	for (size_t i = degreey; i + 1 <= ctrlpoints[0].size(); i++)
	{
		if (v >= knotsy[i])
		{
			ry = static_cast<int>(i);
		}
		if (v < knotsy[i + 1])
		{
			break;
		}
	}
	t = v - knotsy[ry];
	std::vector<double> Basev(4);
	if (ry < 4)
	{
		Basev = FirstInterval(t);
	}
	else if (ry < 5)
	{
		Basev = SecondInterval(t);
	}
	else if (ry < (m - 2))
	{
		Basev = CommonInterval(t);
	}
	else if (ry < (m - 1))
	{
		Basev = RSecondInterval(t);
	}
	else
	{
		Basev = RFirstInterval(t);
	}

	Eigen::Vector3d val(0,0,0);
	for (size_t i = 0; i < degreex + 1; i++)
	{
		for (size_t j = 0; j < degreey + 1; j++)
		{
			val += ctrlpoints[i + rx - 3][j + ry - 3] * Baseu[i] * Basev[j];
		}
	}
	return val;
}

void BSplineSurface::Derivative(const double & u, const double & v, int & rx, int & ry, std::vector<double>& du, std::vector<double>& dv) const
{
	size_t m = ctrlpoints.size();
	rx = 0;

	//rx = (int)(u / len) + degreex;

	for (size_t i = degreex; i + 1 <= ctrlpoints.size(); i++)
	{
		if (u >= knotsx[i])
		{
			rx = static_cast<int>(i);
		}
		if (u < knotsx[i + 1])
		{
			break;
		}
	}
	double t = u - knotsx[rx];
	std::vector<double> Baseu(4);
	std::vector<double> Baseu_D(4);
	if (rx < 4)
	{
		Baseu = FirstInterval(t);
		Baseu_D = FirstInterval_D(t);
	}
	else if (rx < 5)
	{
		Baseu = SecondInterval(t);
		Baseu_D = SecondInterval_D(t);
	}
	else if (rx < (m - 2))
	{
		Baseu = CommonInterval(t);
		Baseu_D = CommonInterval_D(t);
	}
	else if (rx < (m - 1))
	{
		Baseu = RSecondInterval(t);
		Baseu_D = RSecondInterval_D(t);
	}
	else
	{
		Baseu = RFirstInterval(t);
		Baseu_D = RFirstInterval_D(t);
	}

	ry = 0;
	for (size_t i = degreey; i + 1 <= ctrlpoints[0].size(); i++)
	{
		if (v >= knotsy[i])
		{
			ry = static_cast<int>(i);
		}
		if (v < knotsy[i + 1])
		{
			break;
		}
	}
	
	t = v - knotsy[ry];
	std::vector<double> Basev(4);
	std::vector<double> Basev_D(4);
	if (ry < 4)
	{
		Basev = FirstInterval(t);
		Basev_D = FirstInterval_D(t);
	}
	else if (ry < 5)
	{
		Basev = SecondInterval(t);
		Basev_D = SecondInterval_D(t);
	}
	else if (ry < (m - 2))
	{
		Basev = CommonInterval(t);
		Basev_D = CommonInterval_D(t);
	}
	else if (ry < (m - 1))
	{
		Basev = RSecondInterval(t);
		Basev_D = RSecondInterval_D(t);
	}
	else
	{
		Basev = RFirstInterval(t);
		Basev_D = RFirstInterval_D(t);
	}

	du.resize((degreex + 1) * (degreey + 1));
	dv.resize((degreex + 1) * (degreey + 1));
	for (size_t i = 0; i < degreex + 1; i++)
	{
		for (size_t j = 0; j < degreey + 1; j++)
		{
			du[i * (degreey + 1) + j] = Baseu_D[i] * Basev[j];
			dv[i * (degreey + 1) + j] = Baseu[i] * Basev_D[j];
		}
	}
}

Eigen::Vector3d BSplineSurface::DerivativeU(const double & u, const double & v) const
{
	if (ctrlpoints.empty()) return Point();
	BSplineSurface tempsurface;
	std::vector<std::vector<Point>> &derivativepoints = tempsurface.ControlPoints();
	derivativepoints.resize(ctrlpoints.size() - 1, std::vector<Point>(ctrlpoints[0].size()));
	for (size_t j = 0; j < ctrlpoints[0].size(); j++)
	{
		for (size_t i = 1; i < ctrlpoints.size(); i++)
		{
			derivativepoints[i - 1][j] = (ctrlpoints[i][j] - ctrlpoints[i - 1][j]) / (knotsx[i + degreex] - knotsx[i]);
		}
	}
	tempsurface.SetDegree(degreex - 1, degreey);
	return Point();
}

void BSplineSurface::SetDegree(const size_t & degx, const size_t & degy)
{
	degreex = degx;
	degreey = degy;
	UpdateKnots();
}

void BSplineSurface::SetNumberControlPoints(const size_t & mp, const size_t & np)
{
	ctrlpoints.clear();
	ctrlpoints.resize(mp, std::vector<Point>(np, Point(0.0, 0.0, 0.0)));
	UpdateKnots();
}

void BSplineSurface::SetDataPoints(const std::vector<Point>& points)
{
	datapoints = points;
}

void BSplineSurface::SetFlipPoints(const std::vector<Point>& points)
{
	flippoints = points;
}

void BSplineSurface::PreprocessData(void)
{
	std::vector<Point> datanew;
	std::vector<double> lens;
	double totallen = 0.0;
	for (size_t i = 0; i + 1 < datapoints.size(); i++)
	{
		double len = (datapoints[i] - datapoints[i + 1]).norm();
		if (len > 1e-12)
		{
			datanew.push_back(datapoints[i]);
			lens.push_back(len);
			totallen += len;
		}
	}
	double len = (datapoints[datapoints.size() - 1] - datapoints[0]).norm();
	if (len > 1e-12)
	{
		datanew.push_back(datapoints[datapoints.size() - 1]);
		lens.push_back(len);
		totallen += len;
	}
	double curvelen = totallen / 4.0;
	size_t i1 = 0, i2 = 0, i3 = 0;
	double len1 = 0.0;
	for (size_t i = 0; i < datanew.size(); i++)
	{
		len1 += lens[i];
		if (len1 > curvelen)
		{

		}
	}
}

void BSplineSurface::Parameterization(std::vector<double>& t)
{
}

void BSplineSurface::Regression(void)
{
}

void BSplineSurface::UpdateKnots(double interval)
{
	size_t m = ctrlpoints.size();
	if (m == 0) return;
	size_t n = ctrlpoints[0].size();
	assert(m > degreex && n > degreey);
	knotsx.resize(m + degreex + 1);
	knotsy.resize(n + degreey + 1);
	size_t dtx = m - degreex;
	size_t dty = n - degreey;
	len = interval / double(dtx);
	len2 = len*len;
	len3 = len*len*len;
	for (size_t i = 0; i < degreex; i++)
	{
		knotsx[i] = 0.0;
	}
	for (size_t i = degreex; i < m; i++)
	{
		knotsx[i] = interval * static_cast<double>(i - degreex) / dtx;
	}
	for (size_t i = m; i <= m + degreex; i++)
	{
		knotsx[i] = interval;
	}
	for (size_t i = 0; i < degreey; i++)
	{
		knotsy[i] = 0.0;
	}
	for (size_t i = degreey; i < n; i++)
	{
		knotsy[i] = interval * static_cast<double>(i - degreey) / dty;
	}
	for (size_t i = n; i <= n + degreey; i++)
	{
		knotsy[i] = interval;
	}
}

void BSplineSurface::TestSet(void)
{
	degreex = degreey = 3;
	ctrlpoints.clear();
	size_t m = 5;
	size_t n = 5;
	ctrlpoints.resize(m, std::vector<Point>(n));
	srand((unsigned)time(NULL));
	for (size_t i = 0; i < m; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			ctrlpoints[i][j] << (i + rand() / double(RAND_MAX)) * 2.0 / m - 1.0, (j + rand() / double(RAND_MAX)) * 2.0 / n - 1.0, 0.0;
			//ctrlpoints[i][j] << i * 2.0 / m - 1.0, j * 2.0 / n - 1.0, 0.0;
		}
	}
	int i = 0;
	for (const auto & ti : { 0.0, 1.0 / 6.0, 0.5, 5.0 / 6.0, 1.0 })
	{
		int j = 0;
		for (const auto & tj : { 0.0, 1.0 / 6.0, 0.5, 5.0 / 6.0, 1.0 })
		{
			ctrlpoints[i][j] << ti, tj, 0.0;
			j++;
		}
		i++;
	}
	UpdateKnots();
}

void BSplineSurface::SetControlPoints(std::vector<std::vector<Point>> & points)
{
	ctrlpoints = points;
}

void BSplineSurface::SetIdentity()
{
	size_t m = ctrlpoints.size();
	if (m == 0) return;
	size_t n = ctrlpoints[0].size();

	for (int i = 1; i < m+1; ++i)
	{
		for (int j = 1; j < n+1; ++j)
		{
			ctrlpoints[i - 1][j - 1] = Point((knotsx[i] + knotsx[i + 1] + knotsx[i + 2])/3, (knotsy[j] + knotsy[j + 1] + knotsy[j + 2])/3, 0.);
		}
	}
}

std::vector<double> BSplineSurface::CommonInterval(double t) const
{
	double t2 = t*t;
	double t3 = t2*t;
	std::vector<double> value(4);
	value[0] = (1 - 3 * t / len + 3 * t2 / len2 -     t3 / len3) / 6.0;
	value[1] = (4               - 6 * t2 / len2 + 3 * t3 / len3) / 6.0;
	value[2] = (1 + 3 * t / len + 3 * t2 / len2 - 3 * t3 / len3) / 6.0;
	value[3] = (                                      t3 / len3) / 6.0;
	return value;
}
std::vector<double> BSplineSurface::CommonInterval_D(double t) const
{
	double t2 = t*t;
	std::vector<double> value(4);
	value[0] = (- 1.0 / len + 2 * t / len2 -     t2 / len3) / 2.0;
	value[1] = (            - 4 * t / len2 + 3 * t2 / len3) / 2.0;
	value[2] = (  1.0 / len + 2 * t / len2 - 3 * t2 / len3) / 2.0;
	value[3] = (                                 t2 / len3) / 2.0;
	return value;
}
std::vector<double> BSplineSurface::FirstInterval(double t) const
{
	double t2 = t*t;
	double t3 = t2*t;
	std::vector<double> value(4);
	value[0] = (1 - 3 * t / len + 3 * t2 / len2 -     t3 / len3);
	value[1] = (   12 * t / len -18 * t2 / len2 + 7 * t3 / len3) / 4.0;
	value[2] = (                 18 * t2 / len2 -11 * t3 / len3) / 12.0;
	value[3] = (                                      t3 / len3) / 6.0;
	return value;
}
std::vector<double> BSplineSurface::FirstInterval_D(double t) const
{
	double t2 = t*t;
	std::vector<double> value(4);
	value[0] = (-3 / len + 6 * t / len2 - 3 * t2 / len3);
	value[1] = (12 / len -36 * t / len2 +21 * t2 / len3) / 4.0;
	value[2] = (          36 * t / len2 -33 * t2 / len3) / 12.0;
	value[3] = (                              t2 / len3) / 2.0;
	return value;
}
std::vector<double> BSplineSurface::SecondInterval(double t) const
{
	double t2 = t*t;
	double t3 = t2*t;
	std::vector<double> value(4);
	value[0] = (1 - 3 * t / len + 3 * t2 / len2 -     t3 / len3) / 4.0;
	value[1] = (7 + 3 * t / len -15 * t2 / len2 + 7 * t3 / len3) / 12.0;
	value[2] = (1 + 3 * t / len + 3 * t2 / len2 - 3 * t3 / len3) / 6.0;
	value[3] = (                                      t3 / len3) / 6.0;
	return value;
}
std::vector<double> BSplineSurface::SecondInterval_D(double t) const
{
	double t2 = t*t;
	std::vector<double> value(4);
	value[0] = (-3/ len + 6 * t / len2 - 3 * t2 / len3) / 4.0;
	value[1] = (1 / len - 10* t / len2 + 7 * t2 / len3) / 4.0;
	value[2] = (1 / len + 2 * t / len2 - 3 * t2 / len3) / 2.0;
	value[3] = (                             t2 / len3) / 2.0;
	return value;
}
std::vector<double> BSplineSurface::RFirstInterval(double t) const
{
	double t2 = t*t;
	double t3 = t2*t;
	std::vector<double> value(4);
	value[0] = (1 - 3 * t / len + 3 * t2 / len2 -     t3 / len3) / 6.0;
	value[1] = (7 - 3 * t / len - 15* t2 / len2 + 11* t3 / len3) /12.0;
	value[2] = (1 + 3 * t / len + 3 * t2 / len2 - 7 * t3 / len3) / 4.0;
	value[3] = (                                      t3 / len3);
	return value;
}
std::vector<double> BSplineSurface::RFirstInterval_D(double t) const
{
	double t2 = t*t;
	double t3 = t2*t;
	std::vector<double> value(4);
	value[0] = (-1/ len + 2 * t / len2 -      t2 / len3) / 2.0;
	value[1] = (-1/ len - 10* t / len2 + 11 * t2 / len3) / 4.0;
	value[2] = (3 / len + 6 * t / len2 - 21 * t2 / len3) / 4.0;
	value[3] = (                            3*t2 / len3);
	return value;
}
std::vector<double> BSplineSurface::RSecondInterval(double t) const
{
	double t2 = t*t;
	double t3 = t2*t;
	std::vector<double> value(4);
	value[0] = (1 - 3 * t / len + 3 * t2 / len2 -     t3 / len3) / 6.0;
	value[1] = (4               - 6 * t2 / len2 + 3 * t3 / len3) / 6.0;
	value[2] = (2 + 6 * t / len + 6 * t2 / len2 - 7 * t3 / len3) / 12.0;
	value[3] = (                                      t3 / len3) / 4.0;
	return value;
}
std::vector<double> BSplineSurface::RSecondInterval_D(double t) const
{
	double t2 = t*t;
	std::vector<double> value(4);
	value[0] = (-1/ len + 2 * t / len2 -     t2 / len3) / 2.0;
	value[1] = (         -4 * t / len2 + 3 * t2 / len3) / 2.0;
	value[2] = (2 / len + 4 * t / len2 - 7 * t2 / len3) / 4.0;
	value[3] = (                         3 * t2 / len3) / 4.0;
	return value;
}