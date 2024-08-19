// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the BSD 2-Clause License.

/*
 * Simple structs used for representing ECEF and geodetic coordinates.
 */

#ifndef COORDINATE_STRUCTS_HPP
#define COORDINATE_STRUCTS_HPP

struct xyz {
    double x;
    double y;
    double z;
};

struct nva {
    xyz n;
    double alt;
};

struct lla {
    double lat;
    double lon;
    double alt;
};

#endif // COORDINATE_STRUCTS_HPP
