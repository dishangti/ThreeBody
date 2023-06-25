#include <iostream>
#include <cmath>

using namespace std;
const double G = 6.67408e-11;
const double MIN_DIS = 1e-6;

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
int n;
double per, per2, t, display;

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
		ds *= G / (ds2 * dis);
		a += ds * star[i].m;
	}

	return a;
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

R3 RK4_ds(R3 &v0, R3 &acc)
// Order-4 Runge-Kutta method for position
{
	R3 &&a_per2 = acc * per2;
							
	R3 &&k2 = v0 + v0 * per2 + a_per2;		// k1 = v0
	R3 &&k3 = v0 + k2 * per2 + a_per2;		// v1 = v0 + a * dt
	R3 &&k4 = v0 + k3 * per + acc * per;

	return (v0 + k2 * 2 + k3 * 2 + k4) * (per / 6);
}

void update()
{
	for (int i = 0; i < n; i++) {
		Body &st = star[i];
		R3 &&acc = a_func(st.s, i);	// Assume "acc" is a constant during each interval
		st.s = st.s + RK4_ds(st.v, acc);
		st.v += acc * per;
	}
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

		if (is_close(i)) {
			cout << "Position cannot be the same or too close!" << endl;
			return 0;
		}
	}
	
	per2 = per / 2;		// Reduce calc complexation
	 for (double now = per, interval = per; now <= t || interval >= display; now += per, interval += per) {
		update();

		// Annotate this if you don't want crash detection
		if (is_crashed()) {
			cout << "Crashed!" << endl;
			displays(now);
			return 0;
		}
	 	
	 	if (interval >= display) {
			interval = 0;
			displays(now);
		}
		
	}
	
	delete[] star;
	return 0;
}
