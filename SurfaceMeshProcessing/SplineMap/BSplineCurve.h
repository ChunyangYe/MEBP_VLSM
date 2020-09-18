#pragma once
#include <Eigen/Eigen>
class BSplineCurve
{
private:
	typedef Eigen::Vector3d Point;
public:
	BSplineCurve();
	BSplineCurve(const size_t &deg, const size_t &nctrl);
	~BSplineCurve();
	
	Point operator()(const double &t) const;
	const std::vector<Point> & ControlPoints(void) const { return ctrlpoints; }
	const std::vector<Point> & DataPoints(void) const { return datapoints; }
	void SetDegree(const size_t &deg);
	void SetNumberControlPoints(const size_t &np);
	void SetDataPoints(const std::vector<Point> &points);
	void Parameterization(std::vector<double> &t);
	void Regression(void);
private:
	void UpdateKnots(void);
	void TestSet(void);
private:
	size_t degree;
	std::vector<Point> ctrlpoints;
	std::vector<double> knots;
	std::vector<Point> datapoints;
};

