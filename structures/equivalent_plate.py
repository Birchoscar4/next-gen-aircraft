import math

def equivalent_plate_properties(fw, tw, tb, ts, d, E, w, rho, v):
    A_stringer = (d - ts) * tw + 2 * tb * fw
    A_plate = ts * w

    term1 = ((d - ts - 2 * tb) ** 3) * tw / 12
    term2 = 2 * ((tb ** 3 * fw) / 12 + tb * fw * ((d - ts - tb) ** 2) / 4)
    Ix_stringer = term1 + term2

    Ix_plate = w * ts ** 3 / 12

    centroid = ((ts / 2) * A_plate + ((d - ts) / 2 + ts) * A_stringer) / (A_stringer + A_plate)

    dy_stringer = centroid - ((d - ts) / 2 + ts)
    dy_plate = centroid - ts / 2

    I = Ix_plate + A_plate * dy_plate ** 2 + Ix_stringer + A_stringer * dy_stringer ** 2

    te = math.sqrt(12 * I / (A_stringer + A_plate))
    E1e = E * (A_stringer + A_plate) / (w * te)
    rhoe = rho * (A_stringer + A_plate) / (w * te)
    E2e = E * ts / te
    G12 = E / (2 * (1 + v))
    G12e = G12 * ts / te

    return te, E1e, rhoe, E2e, G12e
