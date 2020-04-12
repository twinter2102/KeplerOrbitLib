# Judicious Universe Navigation Instrument v1.0

from math import sin, cos, tan, acos, atan, atan2, sqrt, pi
from math import radians as rad
from math import degrees as deg
import datetime

# ------Global Variables--------

cbody_GM = 3.986004418 * 10 ** 14
cbody_radius = 6371000

# ------Functions---------------

def test():
    print("Judicious Universe Navigation Instrument v1.0\nTest successful")

# Sets the standard gravitational parameter for central body.
def set_cbody_GM(NewGM):
    global cbody_GM
    cbody_GM = NewGM

# Sets the radius for central body.
def set_cbody_radius(NewRadius):
    global cbody_radius
    cbody_radius = NewRadius

# Returns mjd from an instance of datetime.datetime().
def utc_to_mjd(someutc):
    Y = someutc.year
    M = someutc.month
    D = someutc.day

    h = someutc.hour
    m = someutc.minute
    s = someutc.second

    jdn = (2 - (Y // 100) + ((Y // 100) // 4)) + D + ((365.25 * (Y + 4716)) - ((365.25 * (Y + 4716)) % 1))\
        + (30.6001 * (M + 1)) - ((30.6001 * (M + 1)) % 1) - 1524
    jd = jdn + ((h - 12) / 24) + (m / 1440) + (s / 86400)
    mjd = jd - 2400000.5
    return mjd

# Returns utc from a given mjd.
def mjd_to_utc(mjd):
    return datetime.datetime(1858, 11, 17) + datetime.timedelta(days=mjd)

# Returns the current mjd.
def current_mjd():
    return utc_to_mjd(datetime.datetime.utcnow())


# ------Classes-----------------


class Orbit:
    def __init__(self, a, e, i, aop, lan, M0, epoch):
        self.a = a
        self.e = e
        self.i = i
        self.aop = aop
        self.lan = lan
        self.M0 = M0
        self.epoch = epoch

        # Add apoapsis, periapsis, and period

        ap = 0
        pe = 0
        P = 0

    def get_elements(self):
        return self.a, self.e, self.i, self.aop, self.lan, self.M0, self.epoch

    # Returns an instance of Orbit based on cartesian state vectors.
    @staticmethod
    def from_csv(r, rdot, epoch):
        # Calculates orbital momentum and eccentricity vector
        h = Vector.cross(r, rdot)
        ev = Vector.subtract(Vector.divide(Vector.cross(rdot, h), cbody_GM), Vector.divide(r, Vector.norm(r)))

        # Calculates true anomaly and arbitrary vector n
        if Vector.dot(r, rdot) >= 0:
            v = acos(Vector.dot(ev, r) / (Vector.norm(ev) * Vector.norm(r)))
        else:
            v = (2 * pi) - acos(Vector.dot(ev, r) / (Vector.norm(ev) * Vector.norm(r)))
        n = Vector(-h.y, h.x, 0)

        # Calculates inclination eccentricity and Eccentric anomaly
        i = acos(h.z / h.norm())
        e = ev.norm()
        E = 2 * atan((tan(v / 2)) / sqrt((1 + e) / (1 - e)))

        # Calculates longitude of ascending node
        if n.y >= 0:
            lan = acos(n.x / n.norm())
        else:
            lan = (2 * pi) - acos(n.x / n.norm())

        # Calculates argument of periapsis
        if ev.z >= 0:
            aop = acos(Vector.dot(n, ev) / (Vector.norm(n) * Vector.norm(ev)))
        else:
            aop = (2 * pi) - acos(Vector.dot(n, ev) / (Vector.norm(n) * Vector.norm(ev)))

        # Calculates semi-major axis and mean anomaly
        a = 1 / ((2 / r.norm()) - ((rdot.norm() ** 2) / cbody_GM))
        if E >= (e * sin(E)):
            M = E - (e * sin(E))
        else:
            M = (2 * pi) + (E - (e * sin(E)))
        return Orbit(a, e, deg(i), deg(aop), deg(lan), deg(M), epoch)

    # Returns position at time t.
    def get_pos(self, t):
        e = self.e
        a = self.a
        i = rad(self.i)
        lan = rad(self.lan)
        aop = rad(self.aop)
        M0 = rad(self.M0)
        epoch = self.epoch

        # Calculates eccentric anomaly from mean anomaly
        dt = 86400 * (t - epoch)
        M = M0 + (dt * sqrt(cbody_GM / (a ** 3)))
        Eguess = M

        while True:
            Enew = Eguess - (Eguess - e * sin(Eguess) - M) / (1 - e * cos(Eguess))
            if round(Eguess, 8) != round(Enew, 8):
                Eguess = Enew
            else:
                E = Enew
                break

        # Calculates true anomaly and length from central body
        v = 2 * atan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2))
        rc = a * (1 - e * cos(E))

        # Calculates position in orbital frame
        o = Vector.multiply(Vector(cos(v), sin(v), 0), rc)

        # Translates positional components to reference frame
        rx = o.x * (cos(aop) * cos(lan) - sin(aop) * cos(i) * sin(lan)) - o.y \
             * (sin(aop) * cos(lan) + cos(aop) * cos(i) * sin(lan))
        ry = o.x * (cos(aop) * sin(lan) + sin(aop) * cos(i) * cos(lan)) + o.y \
             * (cos(aop) * cos(i) * cos(lan) - sin(aop) * sin(lan))
        rz = o.x * (sin(aop) * sin(i)) + o.y * (cos(aop) * sin(i))

        # Creates position vector and returns it
        r = Vector(rx, ry, rz)
        return r

    # Returns latitude and longitude at time t.
    def get_subpoint(self, t):
        pos = self.get_pos(t)

        # Converts from mjd
        jd = t + 2400000.5
        tu = jd - 2451545.0

        # Calculates era and decides current rotation of earth
        era = (0.7790572732640 + 1.00273781191135448 * tu)

        rot = 2 * pi * (era % 1)

        # Calculates latitude and reference longitude
        lat = atan2(pos.z, sqrt(pos.x ** 2 + pos.y ** 2))
        reflon = atan2(pos.y, pos.x)

        # Adjusts for negative reference longitude
        if reflon < 0:
            reflon = (2 * pi) + reflon

        # Adjusts for longitude above 180 degrees
        lon = (reflon - rot) % (2 * pi)
        if lon > pi:
            lon = lon - (2 * pi)

        return deg(lat), deg(lon)

    # Returns position at time t.
    def get_vel(self, t):
        e = self.e
        a = self.a
        i = rad(self.i)
        lan = rad(self.lan)
        aop = rad(self.aop)
        M0 = rad(self.M0)
        epoch = self.epoch

        # Calculates eccentric anomaly from mean anomaly
        dt = 86400 * (t - epoch)
        M = M0 + (dt * sqrt(cbody_GM / (a ** 3)))
        Eguess = M

        while True:
            Enew = Eguess - (Eguess - e * sin(Eguess) - M) / (1 - e * cos(Eguess))
            if round(Eguess, 8) != round(Enew, 8):
                Eguess = Enew
            else:
                E = Enew
                break

        # Calculates length from central body
        rc = a * (1 - e * cos(E))

        # Calculates velocity vector in orbital frame
        odot = Vector.multiply(Vector(- sin(E), sqrt(1 - e ** 2) * cos(E), 0), sqrt(cbody_GM * a) / rc)

        # Translates velocity components to reference frame
        rdotx = odot.x * (cos(aop) * cos(lan) - sin(aop) * cos(i) * sin(lan)) - odot.y\
            * (sin(aop) * cos(lan) + cos(aop) * cos(i) * sin(lan))
        rdoty = odot.x * (cos(aop) * sin(lan) + sin(aop) * cos(i) * cos(lan)) + odot.y\
            * (cos(aop) * cos(i) * cos(lan) - sin(aop) * sin(lan))
        rdotz = odot.x * (sin(aop) * sin(i)) + odot.y * (cos(aop) * sin(i))

        # Creates velocity vector and returns it
        rdot = Vector(rdotx, rdoty, rdotz)
        return rdot

    def to_csv(self, t):
        e = self.e
        a = self.a
        i = rad(self.i)
        lan = rad(self.lan)
        aop = rad(self.aop)
        M0 = rad(self.M0)
        epoch = self.epoch

        # Calculates eccentric anomaly from mean anomaly
        dt = 86400 * (t - epoch)
        M = M0 + (dt * sqrt(cbody_GM / (a ** 3)))
        Eguess = M

        while True:
            Enew = Eguess - (Eguess - e * sin(Eguess) - M) / (1 - e * cos(Eguess))
            if round(Eguess, 8) != round(Enew, 8):
                Eguess = Enew
            else:
                E = Enew
                break

        # Calculates true anomaly and length from central body
        v = 2 * atan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2))
        rc = a * (1 - e * cos(E))

        # Calculates position and velocity vector in orbital frame
        o = Vector.multiply(Vector(cos(v), sin(v), 0), rc)
        odot = Vector.multiply(Vector(- sin(E), sqrt(1 - e ** 2) * cos(E), 0), sqrt(cbody_GM * a) / rc)

        # Translates positional components to reference frame
        rx = o.x * (cos(aop) * cos(lan) - sin(aop) * cos(i) * sin(lan)) - o.y\
            * (sin(aop) * cos(lan) + cos(aop) * cos(i) * sin(lan))
        ry = o.x * (cos(aop) * sin(lan) + sin(aop) * cos(i) * cos(lan)) + o.y\
            * (cos(aop) * cos(i) * cos(lan) - sin(aop) * sin(lan))
        rz = o.x * (sin(aop) * sin(i)) + o.y * (cos(aop) * sin(i))

        # Translates velocity components to reference frame
        rdotx = odot.x * (cos(aop) * cos(lan) - sin(aop) * cos(i) * sin(lan)) - odot.y\
            * (sin(aop) * cos(lan) + cos(aop) * cos(i) * sin(lan))
        rdoty = odot.x * (cos(aop) * sin(lan) + sin(aop) * cos(i) * cos(lan)) + odot.y\
            * (cos(aop) * cos(i) * cos(lan) - sin(aop) * sin(lan))
        rdotz = odot.x * (sin(aop) * sin(i)) + odot.y * (cos(aop) * sin(i))

        # Creates position and velocity vector and returns them
        r = Vector(rx, ry, rz)
        rdot = Vector(rdotx, rdoty, rdotz)
        return r, rdot


class Spacecraft:
    pass


class Burn:
    def __init__(self, provel, norvel, radvel):
        self.provel = provel
        self.norvel = norvel
        self.radvel = radvel

    # Applies burn to orbit, and returns new orbit instance.
    def apply(self, baseorbit, t):
        r, rdot = baseorbit.to_csv(t)
        pro = self.provel
        nor = self.norvel
        rad = self.radvel

        unitp = Vector.multiply(Vector.divide(rdot, rdot.norm()), pro)
        unitn = Vector.multiply(Vector.divide(Vector.cross(r, rdot), Vector.cross(r, rdot).norm()), nor)
        unitr = Vector.multiply(Vector.divide(r, r.norm()), rad)

        burnvector = Vector.add(Vector.add(unitp, unitn), unitr)
        return Orbit.from_csv(r, Vector.add(rdot, burnvector), t)


# Class for vector handling.
class Vector:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def components(self):
        x = self.x
        y = self.y
        z = self.z
        return x, y, z

    def norm(self):
        n = sqrt((self.x ** 2) + (self.y ** 2) + (self.z ** 2))
        return n

    @staticmethod
    def dot(u, v):
        p = (u.x * v.x) + (u.y * v.y) + (u.z * v.z)
        return p

    @staticmethod
    def cross(u, v):
        x = (u.y * v.z) - (u.z * v.y)
        y = (u.z * v.x) - (u.x * v.z)
        z = (u.x * v.y) - (u.y * v.x)
        return Vector(x, y, z)

    @staticmethod
    def add(u, v):
        x = u.x + v.x
        y = u.y + v.y
        z = u.z + v.z
        return Vector(x, y, z)

    @staticmethod
    def subtract(u, v):
        x = u.x - v.x
        y = u.y - v.y
        z = u.z - v.z
        return Vector(x, y, z)

    @staticmethod
    def multiply(u, s):
        x = u.x * s
        y = u.y * s
        z = u.z * s
        return Vector(x, y, z)

    @staticmethod
    def divide(u, s):
        if s == 0:
            print("Division by zero detected")
            return Vector(0, 0, 0)
        else:
            x = u.x * (s ** (-1))
            y = u.y * (s ** (-1))
            z = u.z * (s ** (-1))
            return Vector(x, y, z)


# Has attributes xs, ys, and zs, where all three are lists. Path in space between ti and tf.
class Path:
    def __init__(self, ti):
        self.proptime = ti
        self.xs = []
        self.ys = []
        self.zs = []

        self.lats = []
        self.lons = []

    def prop(self, proporbit, step):
        self.proptime += step * (1 / 86400)

        x, y, z = proporbit.get_pos(self.proptime).components()
        self.xs.append(x)
        self.ys.append(y)
        self.zs.append(z)

        lat, lon = proporbit.get_subpoint(self.proptime)
        self.lats.append(lat)
        self.lons.append(lon)