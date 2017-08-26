#include <iostream>
#include <utility>

using namespace std;

typedef pair<double, double> dPair;
const double PI = 4 * atan(1);
double toRadians(double degree) {return degree * PI / 180.;}
constexpr double EARTH_EQUATORIAL_RADIUS = 6378137.0; //m
constexpr double EARTH_ECCENTRICITY_SQUARED = 0.00669437999014;

dPair distanceToXY1(dPair gps1, dPair gps2) {
	double earthRadius = 6371000; //meters
	double lon1 = gps1.first;
	double lon2 = gps2.first;
	double lat1 = gps1.second;
	double lat2 = gps2.second;

	double dLat = toRadians(lat2 - lat1);
	double dLng = toRadians(lon2 - lon1);
	double dx = cos(toRadians(lat1)) * cos(toRadians(lat1)) * sin(dLng / 2) * sin(dLng / 2);
	double dy = sin(dLat / 2) * sin(dLat / 2);

	dx = 2 * atan2(sqrt(dx), sqrt(1 - dx));
	dy = -2 * atan2(sqrt(dy), sqrt(1 - dy));
	dx = earthRadius * dx;
	dy = earthRadius * dy;

	return dPair(dx, dy);
}

dPair distanceToLonLat1(double lat1, dPair d) {
	double dx = d.first;
	double dy = d.second;
	double earthRadius = 6371000; //meters

	dx = dx / earthRadius;
	dy = dy / earthRadius;
	dx = tan(dx / 2) * tan(dx / 2) / (1 + tan(dx / 2) * tan(dx / 2));
	dy = tan(dy / 2) * tan(dy / 2) / (1 + tan(dy / 2) * tan(dy / 2));

	double dLat = 2 * asin(sqrt(dy));
	double dLon = 2 * asin(sqrt(dx) / (cos(toRadians(lat1) * cos(toRadians(lat1)))));

	return dPair(dLon, dLat);
}

std::pair<double, double> distanceToXY2(dPair gps1, dPair gps2) {
	double lon1 = gps1.first;
	double lat1 = gps1.second;
	double lon2 = gps2.first;
	double lat2 = gps2.second;
	double lat_radians = lat1 * PI / 180;
	double lat_length = (PI * EARTH_EQUATORIAL_RADIUS * (1.0 - EARTH_ECCENTRICITY_SQUARED)) / (180.0 * pow(1 - (EARTH_ECCENTRICITY_SQUARED * pow(sin(lat_radians), 2)), 1.5));
	double long_length = (PI * EARTH_EQUATORIAL_RADIUS * cos(lat_radians)) / (180.0 * pow(1 - (EARTH_ECCENTRICITY_SQUARED * pow(sin(lat_radians), 2)), 0.5));

	return dPair(-long_length * (lon2 - lon1), -lat_length * (lat2 - lat1));
}

std::pair<double, double> distanceToLonLat2(double lat, dPair d) {
	double dx = d.first;
	double dy = d.second;
	double lat_radians = lat * PI / 180;
	double lat_length = (PI * EARTH_EQUATORIAL_RADIUS * (1.0 - EARTH_ECCENTRICITY_SQUARED)) / (180.0 * pow(1 - (EARTH_ECCENTRICITY_SQUARED * pow(sin(lat_radians), 2)), 1.5));
	double long_length = (PI * EARTH_EQUATORIAL_RADIUS * cos(lat_radians)) / (180.0 * pow(1 - (EARTH_ECCENTRICITY_SQUARED * pow(sin(lat_radians), 2)), 0.5));

	return dPair(-dx / long_length, -dy / lat_length);
}

int main(int argc, char* argv[]) {
	//correct number: 
	//dx = 175.5
	//dy = -2890.2
	//d = 2895.92

	double lon1 = 42.323388;
	double lat1 = -83.433860;
	double lon2 = 42.309647;
	double lat2 = -83.407982;
	
	dPair d = distanceToXY2(dPair(lon1, lat1), dPair(lon2, lat2));
	cout << "dx= " << d.first << " dy=" << d.second << " d=" << sqrt(d.first * d.first + d.second * d.second) << endl;

	d = distanceToLonLat2(lat1, d);
	cout << "Lon2= " << lon1 + d.first << endl;
	cout << "Lat2= " << lat1 + d.second << endl;

	cout << "-------------------------" << endl;
	d = distanceToXY1(dPair(lon1, lat1), dPair(lon2, lat2));
	cout << "dx= " << d.first << " dy=" << d.second << " d=" << sqrt(d.first * d.first + d.second * d.second) << endl;

	d = distanceToLonLat1(lat1, d);
	cout << "Lon2= " << lon1 + d.first << endl;
	cout << "Lat2= " << lat1 + d.second << endl;

	cout << endl << "Standard dx=" << 175.5 << " dy=" << -2890.2 << " d=" << 2895.92 << endl;
	cout << "Standard Lon2 = " << lon2 << " Lat2=" << lat2 << endl;

	system("pause");
}
