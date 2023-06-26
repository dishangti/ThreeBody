#include <iostream>
#include <cmath>
#include <utility>
#include <stdlib>

using namespace std;
const double G = 6.67408e-11;
const double MIN_DIS = 1e-6;
const double MIN_TIME = 1e-6;
const double RK45_EPS = 1e-6;
const double RK45_DT_UPPER_BD = 1e-2;
const double RK45_DT_LOWER_BD = 1e-8;

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
double dt = 0.01;

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

pair<R3, R3> RKF45(R3 &s0, R3 &v0, int except)
// Order-4 Runge-Kutta method for velocity and position
{
	R3 v1z, s1z;

	while (true) {
		R3 &&k1v = a_func(s0, except) * dt;
		R3 &&k1s = v0 * dt;

		R3 &&k2v = a_func(s0 + k1s * (1.0/4.0), except) * dt;
		R3 &&k2s = (v0 + k1v * (1.0/4.0)) * dt;

		R3 &&k3v = a_func(s0 + k1s * (3.0/32.0) + k2s * (9.0/32.0), except) * dt;
		R3 &&k3s = (v0 + k1v * (3.0/32.0) + k2v * (9.0/32.0)) * dt;

		R3 &&k4v = a_func(s0 + k1s * (1932.0/2197.0) - k2s * (7200.0/2197.0) + k3s * (7296.0/2197.0), except) * dt;
		R3 &&k4s = (v0 + k1v * (1932.0/2197.0) - k2v * (7200.0/2197.0) + k3v * (7296.0/2197.0)) * dt;

		R3 &&k5v = a_func(s0 + k1s * (439.0/216.0) - k2s * 8 + k3s * (3680.0/513.0) - k4s * (845.0/4104.0), except) * dt;
		R3 &&k5s = (v0 + k1v * (439.0/216.0) - k2v * 8 + k3v * (3680.0/513.0) - k4v * (845.0/4104.0)) * dt;

		R3 &&k6v = a_func(s0 - k1s * (8.0/27.0) + k2s * 2 - k3s * (3544.0/2565.0) + k4s * (1859.0/4104.0) - k5s * (11.0/40.0), except) * dt;
		R3 &&k6s = (v0 - k1v * (8.0/27.0) + k2v * 2 - k3v * (3544.0/2565.0) + k4v * (1859.0/4104.0) - k5v * (11.0/40.0)) * dt;

		R3 &&v1y = v0 + k1v * (25.0/216.0) + k3v * (1408.0/2565.0) + k4v * (219.0/4104.0) - k5v * (1.0/5.0);
		R3 &&s1y = s0 + k1s * (25.0/216.0) + k3s * (1408.0/2565.0) + k4s * (2197.0/4104.0) - k5s * (1.0/5.0);

		v1z = v0 + k1v * (16.0/135.0) + k3v * (6656.0/12825.0) + k4v * (28561.0/56430.0) - k5v * (9.0/50.0) + k6v * (2.0/55.0);
		s1z = s0 + k1s * (16.0/135.0) + k3s * (6656.0/12825.0) + k4s * (28561.0/56430.0) - k5s * (9.0/50.0) + k6s * (2.0/55.0);

		R3 &&s_deps = s1y - s1z;

		double s_eps = s_deps * s_deps; s_eps = sqrt(s_eps);

		double s_h = sqrt(sqrt((RK45_EPS / 2.0 / s_eps)));

		double n_dt = dt / 2.0 * s_h;
		if (n_dt >= RK45_DT_LOWER_BD && n_dt <= RK45_DT_UPPER_BD) {
			dt = n_dt;
		} else if (n_dt < RK45_DT_LOWER_BD) {
			cout << "Instable solution!"
			exit(0);
		}
			
		
		if (s_eps <= RK45_EPS)
			break;
	}

	return pair<R3, R3>(s1z, v1z);
}

void update()
{
	for (int i = 0; i < n; i++) {
		Body &st = star[i];
		sv_tmp[i] = RKF45(st.s, st.v, i);
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
	
	double now = 0, interval = 0;
	while (now <= t - MIN_TIME) {
		update();

		// Annotate this if you don't want crash detection
		if (is_crashed()) {
			cout << "Crashed!" << endl;
			displays(now);
			return 0;
		}
	 	
		now += dt;
		interval += dt;
	 	if (interval >= display - MIN_TIME) {
			interval = 0;
			displays(now);
		}
		
	}

	delete[] star;
	delete[] sv_tmp;
	return 0;
}
