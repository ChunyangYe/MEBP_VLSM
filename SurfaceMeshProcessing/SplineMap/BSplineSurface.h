#pragma once
#include <Eigen\Eigen>
class BSplineSurface
{
protected:
	typedef Eigen::Vector3d Point;
public:
	BSplineSurface();
	~BSplineSurface();
	Point operator()(const double &u, const double &v) const;
	Point getPoint(const double &u, const double &v);
	double get_inter_len() { return len; }
	void Derivative(const double &u, const double &v, int &rx, int &ry, std::vector<double> &du, std::vector<double> &dv) const;
	Point DerivativeU(const double &u, const double &v) const;
	Point DerivativeV(const double &u, const double &v) const;
	double IsometricDistortion(const double &u, const double &v) const;
	double ConformalDistortion(const double &u, const double &v) const;
	
	const std::vector<std::vector<Point>> & ControlPoints(void) const { return ctrlpoints; }
	std::vector<std::vector<Point>> & ControlPoints(void) { return ctrlpoints; }
	const std::vector<Point> & DataPoints(void) const { return datapoints; }
	const std::vector<Point> & FlipPoints(void) const { return flippoints; }
	size_t DegreeX(void) const { return degreex; }
	size_t DegreeY(void) const { return degreey; }
	void SetDegree(const size_t &degx, const size_t &degy);
	void SetNumberControlPoints(const size_t &mp, const size_t &np);
	void SetDataPoints(const std::vector<Point> &points);
	void SetFlipPoints(const std::vector<Point> &points);
	void PreprocessData(void);
	void Parameterization(std::vector<double> &t);
	void Regression(void);
	void SetControlPoints(std::vector<std::vector<Point>> & points);
	template<typename T>
	T DeBoor(const std::vector<std::vector<T>> &points, const double &u, const double &v) const
	{
		if (points.empty()) return T();
		size_t rx = 0;
		for (size_t i = degreex; i + 1 <= points.size(); i++)
		{
			if (u >= knotsx[i])
			{
				rx = i;
			}
			if (u < knotsx[i + 1])
			{
				break;
			}
		}
		size_t ry = 0;
		for (size_t i = degreey; i + 1 <= points[0].size(); i++)
		{
			if (v >= knotsy[i])
			{
				ry = i;
			}
			if (v < knotsy[i + 1])
			{
				break;
			}
		}
		std::vector<std::vector<T>> d(degreex + 1, std::vector<T>(degreey + 1));
		for (size_t i = 0; i < degreex + 1; i++)
		{
			for (size_t j = 0; j < degreey + 1; j++)
			{
				d[i][j] = points[i + rx - degreex][j + ry - degreey];
			}
		}
		for (size_t jx = 1; jx <= degreex; jx++)
		{
			for (size_t iy = 0; iy <= degreey; iy++)
			{
				for (size_t ix = degreex; ix >= jx; ix--)
				{
					double alphax = (u - knotsx[ix + rx - degreex]) / (knotsx[ix + 1 + rx - jx] - knotsx[ix + rx - degreex]);
					d[ix][iy] = (1 - alphax) * d[ix - 1][iy] + alphax * d[ix][iy];
				}
			}
		}
		for (size_t jy = 1; jy <= degreey; jy++)
		{
			for (size_t iy = degreey; iy >= jy; iy--)
			{
				double alphay = (v - knotsy[iy + ry - degreey]) / (knotsy[iy + 1 + ry - jy] - knotsy[iy + ry - degreey]);
				d[degreex][iy] = (1 - alphay) * d[degreex][iy - 1] + alphay * d[degreex][iy];
			}
		}
		return d[degreex][degreey];
	}

	void SetIdentity();
	void UpdateKnots(double interval = 1.0);

protected:
	
	void TestSet(void);
	std::vector<double> CommonInterval(double t) const;
	std::vector<double> CommonInterval_D(double t) const;
	std::vector<double> FirstInterval(double t) const;
	std::vector<double> FirstInterval_D(double t) const;
	std::vector<double> SecondInterval(double t) const;
	std::vector<double> SecondInterval_D(double t) const;
	std::vector<double> RFirstInterval(double t) const;
	std::vector<double> RFirstInterval_D(double t) const;
	std::vector<double> RSecondInterval(double t) const;
	std::vector<double> RSecondInterval_D(double t) const;
	double len;			
	double len2;
	double len3;
protected:
	size_t degreex, degreey;
	std::vector<double> knotsx, knotsy;
	std::vector<std::vector<Point>> ctrlpoints;
	std::vector<Point> datapoints;
	std::vector<Point> flippoints;
};

