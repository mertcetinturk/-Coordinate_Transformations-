import math
import pyproj

lat = 39.86561
lon = 32.73388

print(f'Latitude ---> {lat}')
print(f'Longitude ---> {lon}')

z = int(31 + ((180 * math.radians(lon)) / (6 * math.pi)))
print(f'Zone ---> {z}')

central_meridian = ((6 * z) - 183) * (math.pi/180)
print(f'Central Meridian ---> {central_meridian}')


# ELLIPSOIDAL PARAMETERS

# We can find these values from AHRS (Attitude and Heading Reference Systems) library
a = 6378137.0  # semi-major axis of the ellipsoid
print(f'Semi-Major Axis of the Ellipsoid: {a}')

b = 6356752.314245179  # semi-minor axis of the ellipsoid
print(f'Semi-Minor Axis of the Ellipsoid: {b}')

# We can also calculate b from f
# f = 1/298.257223563
# b = a * (1 - f)

e_square = ((a ** 2) - (b ** 2)) / (a ** 2)  # first eccentricity

e_square_second = ((a ** 2) - (b ** 2)) / (b ** 2)  # second eccentricity

n = (a - b) / (a + b)

lat_rad = math.radians(lat)
lon_rad = math.radians(lon)

p = (a * (1 - e_square)) / (1 - e_square * (math.sin(lat_rad) ** 2)) ** (3 / 2)  # Radius of curvature in the meridian

v = a / (1 - e_square * (math.sin(lat_rad) ** 2)) ** (1 / 2)  # Radius of curvature in the prime vertical


# A' B' C' D' E'

A = a * ((1 - n) + ((5 / 4) * ((n ** 2) - (n ** 3))) + ((81 / 64) * ((n ** 4) - (n ** 5))))
B = ((3 / 2) * a) * ((n - (n ** 2)) + ((7 / 8) * ((n ** 3) - (n ** 4))) + ((55 / 64) * (n ** 5)))
C = ((15 / 16) * a) * (((n ** 2) - (n ** 3)) + ((3 / 4) * ((n ** 4) - (n ** 5))))
D = ((35 / 48) * a) * (((n ** 3) - (n ** 4)) + ((11 / 16) * (n ** 5)))
E = ((315 / 512) * a) * ((n ** 4) - (n ** 5))

# Meridional arc, the true meridional distance on the ellipsoid from the equator
S = ((A * lat_rad) -
     (B * math.sin(2 * lat_rad)) +
     (C * math.sin(4 * lat_rad)) -
     (D * math.sin(6 * lat_rad)) +
     (E * math.sin(8 * lat_rad)))


# UNIVERSAL TRANSVERSE MERCATOR PROJECTION PARAMETERS

delta_lon = lon_rad - central_meridian  # Difference of longitude from the central meridian
k_zero = 0.9996
FN = 0  # False Northing 0 for the Northern Hemisphere; 10,000.000 for the Southern Hemisphere
FE = 500000  # False Easting


# Terms Used to Calculate General Equations

T1 = S * k_zero

T2 = ((v * math.sin(lat_rad) * math.cos(lat_rad) * k_zero) / 2)

T3 = (((v * math.sin(lat_rad) * (math.cos(lat_rad) ** 3) * k_zero) / 24) *
      (5 - (math.tan(lat_rad) ** 2) +
       ((9 * e_square_second) * (math.cos(lat_rad) ** 2)) +
       (4 * (e_square_second ** 2) * (math.cos(lat_rad) ** 4))))

T4 = (((v * math.sin(lat_rad) * math.cos(lat_rad) ** 5 * k_zero) / 720) *
      (61 - 58 * (math.tan(lat_rad) ** 2) +
       (math.tan(lat_rad) ** 4) + 270 * e_square_second * (math.cos(lat_rad) ** 2) -
       330 * (math.tan(lat_rad) ** 2) * e_square_second * (math.cos(lat_rad) ** 2) +
       445 * (e_square_second ** 2) *
       (math.cos(lat_rad) ** 2) +
       324 * (e_square_second ** 3) * (math.cos(lat_rad) ** 6) -
       680 * (math.tan(lat_rad) ** 2) *
       (e_square_second ** 2) * (math.cos(lat_rad) ** 4) +
       88 * (e_square_second ** 4) - 600 * (math.tan(lat_rad) ** 2) *
       (e_square_second ** 3) * (math.cos(lat_rad) ** 6) -
       192 * (math.tan(lat_rad) ** 2) * (e_square_second ** 4) *
       (math.cos(lat_rad) ** 8))
      )

T5 = (((v * math.sin(lat_rad) * (math.cos(lat_rad) ** 7) * k_zero) / 40320) *
      (1385 - (3111 * (math.tan(lat_rad) ** 2)) + (543 * (math.tan(lat_rad) ** 4)) - (math.tan(lat_rad) ** 6)))

T6 = (v * math.cos(lat_rad) * k_zero)

T7 = (((v * (math.cos(lat_rad) ** 3) * k_zero) / 6) *
      (1 - (math.tan(lat_rad) ** 2) + (e_square_second * (math.cos(lat_rad) ** 2))))

T8 = (((v * (math.cos(lat_rad) ** 5) * k_zero) / 120) *
      (5 - (18 * (math.tan(lat_rad) ** 2)) + (math.tan(lat_rad) ** 4) +
       (14 * e_square_second * (math.cos(lat_rad) ** 2)) -
       (58 * (math.tan(lat_rad) ** 2) * e_square_second * (math.cos(lat_rad) ** 2)) +
       (13 * (e_square_second ** 2) * (math.cos(lat_rad) ** 4)) +
       (4 * (e_square_second ** 3) * (math.cos(lat_rad) ** 6)) -
       (64 * (math.tan(lat_rad) ** 2) * (e_square_second ** 2) * (math.cos(lat_rad) ** 4)) -
       (24 * (math.tan(lat_rad) ** 2) * (e_square_second ** 3) * (math.cos(lat_rad) ** 6))))

T9 = (((v * (math.cos(lat_rad) ** 7) * k_zero) / 5040) *
      (61 - (479 * (math.tan(lat_rad) ** 2)) + (179 * (math.tan(lat_rad) ** 4)) - (math.tan(lat_rad) ** 6)))


# CONVERSION OF GEOGRAPHIC COORDINATES TO GRID COORDINATES

# The general formulas for the computation of N and E are:
N = FN + (T1 + ((delta_lon ** 2) * T2) + ((delta_lon ** 4) * T3) + ((delta_lon ** 6) * T4) + ((delta_lon ** 8) * T5))
E = FE + ((delta_lon * T6) + ((delta_lon ** 3) * T7) + ((delta_lon ** 5) * T8) + ((delta_lon ** 7) * T9))

print('\nNorth is ---> ', N)
print('East is ---> ', E)


def latlon_to_utm(lat, lon):
    return pyproj.Transformer.from_crs("EPSG:4326", "EPSG:32636").transform(lat, lon)


utm_coordinates = latlon_to_utm(lat, lon)

print(f'\npyproj North ---> {utm_coordinates[1]}')
print(f'pyproj East ---> {utm_coordinates[0]}\n')

print(f'● Difference between Finding Northing and Proj Northing: {N - utm_coordinates[1]}')

if abs(N - utm_coordinates[1]) < 1 * (10 ** (-5)):
    print(f'\n⮞Is difference smaller then 1.0 x (10**-5)? {"Yes"}')

print("----------------------------------------------------------------------")
print(f'● Difference between Finding Easting and Proj Easting: {E - utm_coordinates[0]}')

if abs(E - utm_coordinates[0]) < 1 * (10 ** (-5)):
    print(f'\n⮞Is difference smaller then 1.0 x (10**-5)? {"Yes"}')
