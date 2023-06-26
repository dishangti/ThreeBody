#include <iostream>
#include <cmath>
#include <utility>

using namespace std;
const double G = 6.67408e-11;
const double MIN_DIS = 1e-6;
const double MIN_TIME = 1e-6;

struct R3 {
	double x, y, z;
	
	R3(double sx, double sy, double sz): x(sx), y(sy), z(sz) { }
	
	R3():x(0), y(0), z(0) { }
	
	inline R3 &operator=(const R3 &r1) {
		x = r1.x; y = r1.y; z = r1.z;
		return *this;
	}
	
	inline R3 operator/(double n1) {
		return R3(x / n1, y / n1, z / n1);
	}
	
	inline R3 operator*(double n1) {
		return R3(x * n1, y * n1, z * n1);
	}
	
	inline double operator*(const R3 &r1) {
		return x * r1.x + y * r1.y + z * r1.z;
	}
	
	inline R3 operator+(const R3 &r1) {
		return R3(x + r1.x, y + r1.y, z + r1.z);
	}
	
	inline R3 operator-(const R3 &r1) {
		return R3(x - r1.x, y - r1.y, z - r1.z);
	}
	
	inline R3 operator-() {
		return R3(-x, -y, -z);
	}
	
	inline R3 &operator+=(const R3 &r1) {
		*this = *this + r1;
		return *this;
	}
	
	inline R3 &operator*=(double n1) {
		*this = *this * n1;
		return *this;
	}

	inline R3 &operator/=(double n1) {
		*this = *this / n1;
		return *this;
	}
};

struct Body {
	R3 s, v;
	double m;

	inline bool operator==(const Body &s1) {
		R3 &&ds = this->s - s1.s;
		double dis = sqrt(ds * ds);
		return dis <= MIN_DIS;
	}
};
Body *star;
pair<R3, R3> *sv_tmp;
int n;
double per, per2;

R3 a_func(R3 p, int except)
// Function of acceleration
// "except" is the number of the star being calculated
{
	R3 a;

	for (int i = 0; i < n; i++) {
		if (i == except) continue;

		R3 &&ds = star[i].s - p;
		double ds2 = ds * ds;
		double dis = sqrt(ds2);
		ds /= ds2 * dis;
		a += ds * star[i].m;
	}

	return a * G;
}

pair<R3, R3> RK4(R3 &s0, R3 &v0, int except)
// Order-4 Runge-Kutta method for velocity and position
{
	R3 &&k1v = a_func(s0, except);
	R3 &&k1s = move(v0);

	R3 &&k2v = a_func(s0 + k1s * per2, except);
	R3 &&k2s = v0 + k1v * per2;

	R3 &&k3v = a_func(s0 + k2s * per2, except);
	R3 &&k3s = v0 + k2v * per2;

	R3 &&k4v = a_func(s0 + k3s * per, except);
	R3 &&k4s = v0 + k3v * per;

	return pair<R3, R3>(s0 + (k1s + (k2s + k3s) * 2 + k4s) * (per / 6),
						v0 + (k1v + (k2v + k3v) * 2 + k4v) * (per / 6));
}

void update()
{
	for (int i = 0; i < n; i++) {
		Body &st = star[i];
		sv_tmp[i] = RK4(st.s, st.v, i);
	}

	for (int i = 0; i < n; i++) {
		Body &st = star[i];
		st.s = sv_tmp[i].first;
		st.v = sv_tmp[i].second;
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
	for (int i = 0; i < last; i++)
		if (star[i] == star[last]) return true;
	return false;
}

bool is_crashed()
{
	for (int i = 0; i < n; i++)
		if (is_close(i)) return true;
	return false;
}

int main()
{
	double t, display;

	cout << "Number of stars: ";
	cin >> n;
	cout << "Calculation interval(sec): ";
	cin >> per;
	cout << "Finish time(sec): ";
	cin >> t;
	cout << "Display interval(sec): ";
	cin >> display;
	
	star = new Body[n];
	sv_tmp = new pair<R3, R3>[n];
	if (!star || !sv_tmp) {
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

		if (is_close(i)) {
			cout << "Position cannot be the same or too close!" << endl;
			return 0;
		}
	}
	
	per2 = per / 2;		// Reduce calc complexation
	double now = per, interval = per;
	for (; now <= t - MIN_TIME; now += per, interval += per) {
		update();

		// Annotate this if you don't want crash detection
		if (is_crashed()) {
			cout << "Crashed!" << endl;
			displays(now);
			return 0;
		}
	 	
	 	if (interval >= display - MIN_TIME) {
			interval = 0;
			displays(now);
		}
		
	}
	displays(now);

	delete[] star;
	delete[] sv_tmp;
	return 0;
}
