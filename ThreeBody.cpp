#include <iostream>

#include <cmath>

using namespace std;

const double G = 6.67408e-11;

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

};

Body *star;

int n;

double per, t, display, kper;

void updatea()

{

	R3 inita(0, 0, 0);

	for (int i = 0; i < n; i++) star[i].a = inita;

	

	for (int i = 0; i < n; i++) {

		for (int j = i + 1; j < n; j++) {

			R3 tvct = star[j].s - star[i].s;

			double dis = sqrt(tvct * tvct);

			tvct = tvct * G / (dis * dis * dis);

			star[i].a += tvct * star[j].m;

			star[j].a += -tvct * star[i].m;

		}

	}

}

void updatev()

{

	for (int i = 0; i < n; i++) {

		star[i].v += star[i].a * per;

	}

}

void updates()

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

int main()

{

	cout << "Number of stars:";

	cin >> n;

	cout << "Time interval for calculating(sec):";

	cin >> per;

	cout << "Finish time(sec):";

	cin >> t;

	cout << "Time interval for display(sec):";

	cin >> display;

	

	star = new Body[n];

	if (!star) {

		cout << "Initializing failed." << endl;

		return 0;

	}

	

	cout << endl; 

	for (int i = 0; i < n; i++) {

		cout << "Star " << i << ':' << endl;

		cout << "Mass:";

		cin >> star[i].m;

		cout << "Position:";

		cin >> star[i].s.x >> star[i].s.y >> star[i].s.z;

		cout << "Velocity:";

		cin >> star[i].v.x >> star[i].v.y >> star[i].v.z;

		cout << endl;

	}

	

	kper = 0.5 * per * per;

	double interval = per;

	 for (double now = per; now <= t || interval >= display; now += per)  {

	 	updatea();

	 	updates();

	 	updatev();

	 	

	 	if (interval >= display) {

			interval = 0;

			displays(now);

		}

		interval += per;

	}

	

	delete[] star;

	return 0;

}
