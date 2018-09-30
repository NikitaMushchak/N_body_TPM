#pragma once
#include <vector>
#include <cstdio>
#define BOOL int

#define Vect const std::vector<double>&

class vect
{
public:
	vect()								 {}
	virtual ~vect();
	//vect(double x, double y, double z) { Set(x, y, z); }
	std::vector<double> r(double x, double y, double z) { Set(x,y,z); };
	vect(Vect r) { x = r[0]; y = r[1]; z = r[2]; }
	
	void Set(double x, double y, double z) { this->x = x; this->y = y; this->z = z; }
	
	vect& operator +=(Vect r)   { x += r[0]; y += r[1]; z += r[2];	}
	vect& operator -=(Vect r)   { x -= r[0]; y -= r[1]; z -= r[2];	}
	vect& operator *=(Vect r)   { x *= r[0]; y *= r[1]; z *= r[2];	}
	vect& operator /=(Vect r)   { x /= r[0]; y /= r[1]; z /= r[2];	}

	//vect operator -() { vect r(-x, -y, -z); return r; }
	BOOL operator ==(Vect r)			const { return x == r[0] && y == r[1] && z == r[2]; }
	BOOL operator !=(Vect r)			const { return !operator==(r); }


	double Sqr()						const { return x * x + y * y + z * z; }
	double Abs()						const;
	double GetNorm1()					const;


	vect& Sum(Vect r1, Vect r2)
	{
		x = r1[0] + r2[0]; y = r1[1] + r2[1]; z = r1[2] + r2[2]; return *this;
	}

	vect& AddSum(Vect r1, Vect r2)
	{
		x += r1[0] + r2[0]; y += r1[1] + r2[1]; z += r1[2] + r2[2]; return *this;
	}

	vect& SubtractSum(Vect r1, Vect r2)
	{
		x -= r1[0] + r2[0]; y -= r1[1] + r2[1]; z -= r1[2] + r2[2]; return *this;
	}

	vect& Difference(Vect r1, Vect r2)
	{
		x = r1[0] - r2[0]; y = r1[1] - r2[1]; z = r1[2] - r2[2]; return *this;
	}

	vect& AddDifference(Vect r1, Vect r2)
	{
		x += r1[0] - r2[0]; y += r1[1] - r2[1]; z += r1[2] - r2[2]; return *this;
	}

	vect& Product(Vect r0, double k)
	{
		x = r0[0] * k; y = r0[1] * k; z = r0[2] * k; return *this;
	}
	// method multiplication of the Vector and scalar
	vect& Product(double k, Vect r0)
	{
		x = r0[0] * k; y = r0[1] * k; z = r0[2] * k; return *this;
	}
	//
	vect& AddProduct(Vect r0, double k)
	{
		x += r0[0] * k; y += r0[1] * k; z += r0[2] * k; return *this;
	}

	vect& AddProduct(double k, Vect r0)
	{
		x += r0[0] * k; y += r0[1] * k; z += r0[2] * k; return *this;
	}

	vect& SubtractProduct(Vect r0, double k)
	{
		x -= r0[0] * k; y -= r0[1] * k; z -= r0[2] * k; return *this;
	}

	vect& SubtractProduct(double k, Vect r0)
	{
		x -= r0[0] * k; y -= r0[1] * k; z -= r0[2] * k; return *this;
	}

	vect& Quotient(Vect r0, double k)
	{
		x = r0[0] / k; y = r0[1] / k; z = r0[2] / k; return *this;
	}

	vect& VectorProduct(vect r1, vect r2)
	{
		x = r1.y * r2.z - r1.z * r2.y; y = r1.z * r2.x - r1.x * r2.z;
		z = r1.x * r2.y - r1.y * r2.x; return *this;
	}

	vect& AddVectorProduct(vect r1, vect r2)
	{
		x += r1.y * r2.z - r1.z * r2.y; y += r1.z * r2.x - r1.x * r2.z;
		z += r1.x * r2.y - r1.y * r2.x;	return *this;
	}

	double x, y, z;

//private:
//	double rx, ry, rz;
//	double vx, vy, vz;
//	double ax, ay, az;
//	std::vector<double> r{ rx, ry, rz };
//	std::vector<double> v{ vx, vy, vz };
//	std::vector<double> a{ ax, ay, az };

};

static const vect VECT_NULL(0., 0., 0.); //NULL vector

vect::vect()
{
}


vect::~vect()
{
}
