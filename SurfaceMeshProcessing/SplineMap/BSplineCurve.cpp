#include "BSplineCurve.h"
#include <ctime>


BSplineCurve::BSplineCurve(void)
	:degree(3)
{
	//TestSet();
}

BSplineCurve::BSplineCurve(const size_t & deg, const size_t &nctrl)
	:degree(deg),
	ctrlpoints(nctrl)
{
	UpdateKnots();
}

BSplineCurve::~BSplineCurve(void)
{
}

Eigen::Vector3d BSplineCurve::operator()(const double & t) const
{
	size_t r = 0;
	for (size_t i = degree; i + 1 <= ctrlpoints.size(); i++)
	{
		if (t >= knots[i])
		{
			r = i;
		}
		if (t < knots[i + 1])
		{
			break;
		}
	}
	std::vector<Eigen::Vector3d> d(degree + 1);
	for (size_t i = 0; i < degree + 1; i++)
	{
		d[i] = ctrlpoints[i + r - degree];
	}
	for (size_t j = 1; j <= degree; j++)
	{
		for (size_t i = degree; i >= j; i--)
		{
			double alpha = (t - knots[i + r - degree]) / (knots[i + 1 + r - j] - knots[i + r - degree]);
			d[i] = (1 - alpha) * d[i - 1] + alpha * d[i];
		}
	}
	return d[degree];
}

void BSplineCurve::SetDegree(const size_t & deg)
{
	degree = deg;
	UpdateKnots();
}

void BSplineCurve::SetNumberControlPoints(const size_t & np)
{
	ctrlpoints.resize(np);
	UpdateKnots();
}

void BSplineCurve::SetDataPoints(const std::vector<Point>& points)
{
	datapoints = points;
}

void BSplineCurve::Parameterization(std::vector<double>& t)
{
	t.resize(datapoints.size());
	double totallen = 0.0;
	for (size_t i = 0; i + 1 < datapoints.size(); i++)
	{
		totallen += (datapoints[i] - datapoints[i + 1]).norm();
	}
	double currentlen = 0.0;
	t[0] = 0.0;
	for (size_t i = 0; i + 1 < datapoints.size(); i++)
	{
		currentlen += (datapoints[i] - datapoints[i + 1]).norm();
		t[i + 1] = currentlen / totallen;
	}
}

void BSplineCurve::Regression(void)
{
	if (datapoints.empty())
	{
		return;
	}
	std::vector<double> t;
	Parameterization(t);
	Eigen::MatrixXd A(datapoints.size(), ctrlpoints.size());
	
	for (size_t i = 0; i < datapoints.size(); i++)
	{
		for (size_t j = 0; j < ctrlpoints.size(); j++)
		{
			size_t r = 0;
			for (size_t k = degree; k + 1 <= ctrlpoints.size(); k++)
			{
				if (t[i] >= knots[k])
				{
					r = k;
				}
				if (t[i] < knots[k + 1])
				{
					break;
				}
			}
			std::vector<double> d(degree + 1, 0.0);
			if (j + degree >= r && j <= r)
			{
				d[j + degree - r] = 1.0;
				for (size_t k = 1; k <= degree; k++)
				{
					for (size_t l = degree; l >= k; l--)
					{
						double alpha = (t[i] - knots[l + r - degree]) / (knots[l + 1 + r - k] - knots[l + r - degree]);
						d[l] = (1 - alpha) * d[l - 1] + alpha * d[l];
					}
				}
			}
			A(i, j) = d[degree];
		}
	}
	Eigen::VectorXd bx(datapoints.size()), by(datapoints.size());
	for (size_t i = 0; i < datapoints.size(); i++)
	{
		bx(i) = datapoints[i](0);
		by(i) = datapoints[i](1);
	}
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
	Eigen::VectorXd x = svd.solve(bx);
	Eigen::VectorXd y = svd.solve(by);
	for (size_t i = 0; i < ctrlpoints.size(); i++)
	{
		ctrlpoints[i][0] = x(i);
		ctrlpoints[i][1] = y(i);
		ctrlpoints[i][2] = 0.0;
	}
}

void BSplineCurve::UpdateKnots(void)
{
	size_t n = ctrlpoints.size();
	assert(n > degree);
	knots.resize(n + degree + 1);
	size_t dt = n - degree;
	for (size_t i = 0; i < degree; i++)
	{
		knots[i] = 0.0;
	}
	for (size_t i = degree; i < n; i++)
	{
		knots[i] = static_cast<double>(i - degree) / dt;
	}
	for (size_t i = n; i <= n + degree; i++)
	{
		knots[i] = 1.0;
	}
}

void BSplineCurve::TestSet(void)
{
	degree = 3;
	ctrlpoints.clear();
	size_t n = 7;
	srand((unsigned)time(NULL));
	for (size_t i = 0; i < n; i++)
	{
		Eigen::Vector3d p;
		//p << std::cos(i * 2.0 * M_PI / n), std::sin(i * 2.0 * M_PI / n), 0.0;
		p << rand() / double(RAND_MAX) * 2.0 - 1.0, rand() / double(RAND_MAX) * 2.0 - 1.0, 0.0;
		ctrlpoints.push_back(p);
	}
	UpdateKnots();
}
