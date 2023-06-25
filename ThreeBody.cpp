#include <iostream>
#include <cmath>

using namespace std;
const double G = 6.67408e-11;
const double MIN_DIS = 1e-6;

struct R3 {
	double x, y, z;
	
	R3(double sx, double sy, double sz): x(sx), y(sy), z(sz) { }
	
	R3():x(0), y(0), z(0) { }
	
	inline R3 operator=(R3 r1) {
		x = r1.x; y = r1.y; z = r1.z;
		return *this;
	}
	
	inline R3 operator/(double n1) {
		return R3(x / n1, y / n1, z / n1);
	}
	
	inline R3 operator*(double n1) {
		return R3(x * n1, y * n1, z * n1);
	}
	
	inline double operator*(R3 r1) {
		return x * r1.x + y * r1.y + z * r1.z;
	}
	
	inline R3 operator+(R3 r1) {
		return R3(x + r1.x, y + r1.y, z + r1.z);
	}
	
	inline R3 operator-(R3 r1) {
		return R3(x - r1.x, y - r1.y, z - r1.z);
	}
	
	inline R3 operator-() {
		return R3(-x, -y, -z);
	}
	
	inline R3 operator+=(R3 r1) {
		*this = *this + r1;
		return *this;
	}
	
	inline R3 operator*=(double n1) {
		*this = *this * n1;
		return *this;
	}
};

struct Body {
	R3 s, v, a;
	double m;

	inline bool operator==(Body s1) {
		R3 ds = this->s - s1.s;
		double dis = sqrt(ds * ds);
		return dis <= MIN_DIS;
	}
};
Body *star;
int n;
double per, t, display, kper;

void update_a()
{
	R3 init_a;
	for (int i = 0; i < n; i++) star[i].a = init_a;
	
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			R3 ds = star[j].s - star[i].s;
			double ds2 = ds * ds;
			double dis = sqrt(ds2);
			ds *= G / (ds2 * dis);
			star[i].a += ds * star[j].m;
			star[j].a += -ds * star[i].m;
		}
	}
}

void update_v()
{
	for (int i = 0; i < n; i++) {
		star[i].v += star[i].a * per;
	}
}

void update_s()
{
	for (int i = 0; i < n; i++) {
		star[i].s += star[i].v * per + star[i].a * kper;
	}
}

void displays(double time)
{
	cout << "At " << time << " sec" << endl;
	for (int i = 0; i < n; i++) {
		cout << "Star " << i << ": "
				 << '(' << star[i].s.x << ", "
				 << star[i].s.y << ", "
				 << star[i].s.z << ')' << endl;
	}
	cout << endl;
}

bool is_close(int last)
{
	for (int i = 0; i < last - 1; i++)
		if (star[i] == star[last]) return true;
	return false;
}

bool is_crashed()
{
	for (int i = 1; i < n; i++)
		if (is_close(i)) return true;
	return false;
}

int main()
{
	cout << "Number of stars: ";
	cin >> n;
	cout << "Calculation interval(sec): ";
	cin >> per;
	cout << "Finish time(sec): ";
	cin >> t;
	cout << "Display interval(sec): ";
	cin >> display;
	
	star = new Body[n];
	if (!star) {
		cout << "Initializing failed!" << endl;
		return 0;
	}
	
	cout << endl;
	for (int i = 0; i < n; i++) {
		cout << "Star " << i << ':' << endl;
		cout << "Mass: ";
		cin >> star[i].m;
		cout << "Position: ";
		cin >> star[i].s.x >> star[i].s.y >> star[i].s.z;
		cout << "Velocity: ";
		cin >> star[i].v.x >> star[i].v.y >> star[i].v.z;
		cout << endl;

		if (is_close(n)) {
			cout << "Position cannot be the same or too close!" << endl;
			return 0;
		}
	}
	
	kper = 0.5 * per * per;		// Reduce calc complexation
	double interval = per;
	 for (double now = per; now <= t || interval >= display; now += per)  {
	 	update_a();
	 	update_s();
	 	update_v();

		// Annotate this if you don't want crash detection
		if (is_crashed()) {
			cout << "Crashed!" << endl;
			return 0;
		}
	 	
	 	if (interval >= display) {
			interval = 0;
			displays(now);
		}
		interval += per;
	}
	
	delete[] star;
	return 0;
}
