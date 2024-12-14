import numpy as np
import backend.structure.residue_constants as rc


def calculate_o_pos(n_pos, ca_pos, c_pos):
    cca = ca_pos - c_pos
    u_cca = cca / np.linalg.norm(cca)

    cn = n_pos - c_pos
    u_cn = cn / np.linalg.norm(cn)

    u_co = adjust_angle(u_cca, u_cn, 92.4)
    return (u_co + c_pos) * rc.c_o


def adjust_angle(u_cca, u_co, alpha):
    alpha = np.deg2rad(alpha)

    cos_theta = np.dot(u_cca, u_co)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    theta = np.arccos(cos_theta)  # current angle
    delta_alpha = alpha - theta  # rotation angle

    r = np.cross(u_cca, u_co)  # rotation axis

    r_unit = r / np.linalg.norm(r)

    return (
            u_co * np.cos(delta_alpha) +
            np.cross(r_unit, u_co) * np.sin(delta_alpha) +
            r_unit * np.dot(r_unit, u_co) * (1 - np.cos(delta_alpha))
    )


def main():
    for key, data in rc.rigid_group_atom_positions.items():
        n = np.array(data[0][2])
        ca = np.array(data[1][2])
        c = np.array(data[2][2])

        o = np.round(calculate_o_pos(n, ca, c), 3)
        print(f"{key} {o[0]:>6.3f}, {o[1]:>6.3f}, {o[2]:>6.3f}")



if __name__ == '__main__':
    main()