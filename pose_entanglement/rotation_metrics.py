import numpy as np
from scipy.spatial.transform import Rotation as R
import jax.numpy as jnp


def geodesic_angle_deg(R1, R2):
    R_rel = jnp.matmul(jnp.transpose(R1, (0, 2, 1)), R2)
    trace = jnp.trace(R_rel, axis1=1, axis2=2)
    angle_rad = jnp.arccos(jnp.clip((trace - 1) / 2, -1.0, 1.0))
    return jnp.degrees(angle_rad)

def rotation_geodesic_angle(R1, R2):
    """
    Geodesic angle(s) on SO(3) between rotation matrices R1 and R2.
    Inputs: arrays with shape (3,3) or (...,3,3).
    Returns: scalar or array of angles in radians (same leading shape).
    """
    R1 = np.asarray(R1, dtype=float)
    R2 = np.asarray(R2, dtype=float)
    if R1.shape[-2:] != (3,3) or R2.shape[-2:] != (3,3):
        raise ValueError("R1 and R2 must have shape (...,3,3)")
    Rrel = np.matmul(np.swapaxes(R1, -1, -2), R2)   # R1^T @ R2
    trace = np.einsum('...ii', Rrel)
    cos_theta = (trace - 1.0) / 2.0
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    return np.arccos(cos_theta)


def geodesic_from_identity(R):
    """
    Geodesic distance (radians, in [0, pi]) from the identity for rotation matrix/matrices R.
    Accepts R shape (3,3) or (...,3,3). Returns scalar or array of angles with leading shape.
    """
    R = np.asarray(R, dtype=float)
    if R.shape[-2:] != (3,3):
        raise ValueError("R must have shape (...,3,3)")
    trace = np.einsum("...ii", R)                # sum of diagonal entries
    cos_theta = (trace - 1.0) / 2.0
    cos_theta = np.clip(cos_theta, -1.0, 1.0)   # numeric stability
    return np.arccos(cos_theta)


def compute_symmetry_matrices(n_sym: int) -> np.ndarray:
    R_sym = np.zeros((n_sym, 3, 3))
    for idx in range(n_sym):
        R_sym[idx] = R.from_euler('z', idx * 360 / n_sym, degrees=True).as_matrix()
    return R_sym