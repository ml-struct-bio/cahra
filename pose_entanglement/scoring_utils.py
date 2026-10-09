import numpy as np
from scipy.spatial.transform import Rotation as R
import pandas as pd
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, roc_auc_score
from dataclasses import dataclass
from typing import List

from rotation_metrics import rotation_geodesic_angle, geodesic_from_identity


def compute_state_metrics(s_pred: np.ndarray, s_gt: np.ndarray) -> dict:
    """
    s_pred: predicted scores (float, typically in [0,1])
    s_gt: ground truth states (0 or 1)
    Returns dictionary of metrics.
    """
    s_pred = np.asarray(s_pred, dtype=float).ravel()
    s_rounded = np.round(s_pred)
    s_gt = np.asarray(s_gt, dtype=float).ravel()
    if s_pred.shape != s_gt.shape:
        raise ValueError("Prediction and ground-truth vectors must have same shape")
    N = s_pred.size
    S = float(np.mean(s_pred))
    S_star = float(np.mean(s_gt))
    global_soft = abs(S - S_star)
    S_hard = float(np.mean(s_rounded))
    global_hard = abs(S_hard - S_star)
    local_hard = float(np.mean(np.abs(s_rounded - s_gt)))
    local_soft = float(np.mean(np.abs(s_pred - s_gt)))

    tn, fp, fn, tp = confusion_matrix(s_gt, s_rounded).ravel()
    accuracy  = accuracy_score(s_gt, s_rounded)
    precision = precision_score(s_gt, s_rounded)
    recall    = recall_score(s_gt, s_rounded)
    auc = roc_auc_score(s_gt, s_pred)
    # fpr, tpr, thresholds = roc_curve(s_gt, s_pred)
    
    from sklearn.metrics import log_loss

    state_cross_entropy_soft = log_loss(s_gt, s_pred)
    state_cross_entropy_hard = log_loss(s_gt, s_rounded)

    return {
        "S_pred_hard": float(S_hard),
        "S_pred_soft": float(S),
        "S_gt": float(S_star),
        "state_global_l1_soft": global_soft,
        "state_global_l1_hard": global_hard,
        "state_local_l1_hard": local_hard,
        "state_local_l1_soft": local_soft,
        "state_tn": int(tn),
        "state_fp": int(fp),
        "state_fn": int(fn),
        "state_tp": int(tp),
        "state_accuracy": float(accuracy),
        "state_precision": float(precision),
        "state_recall": float(recall),
        "state_auc": float(auc),
        "state_cross_entropy_soft": float(state_cross_entropy_soft),
        "state_cross_entropy_hard": float(state_cross_entropy_hard),

    }


def compute_pose_metrics(rotation_pred, rotation_gt):
    angles = np.rad2deg(rotation_geodesic_angle(rotation_pred, rotation_gt))

    return {'score_rot1': float(np.mean(angles)),
            'score_rot2': float(np.mean(angles)),
            'score_rot3': float(np.mean(angles)),
            'score_rot4': float(np.mean(angles)),
            # 'rot_median_error': float(np.median(angles)),
            # 'rot_75th_percentile_error': float(np.percentile(angles, 75)),
            # 'rot_90th_percentile_error': float(np.percentile(angles, 90)),
            }


def compute_global_rotation(initial_pose_mat, refinement_pose_mat0, refinement_pose_mat1, R_sym, percentage=0.5):
    n_sym = R_sym.shape[0]
    relative_rotation0 = np.zeros((n_sym, len(initial_pose_mat), 3, 3))
    relative_rotation1 = np.zeros((n_sym, len(initial_pose_mat), 3, 3))
    for idx in range(n_sym):
        relative_rotation0[idx] = R_sym[idx].reshape(1,3,3) @ initial_pose_mat @ refinement_pose_mat0.transpose(0,2,1)
        relative_rotation1[idx] = R_sym[idx].reshape(1,3,3) @ initial_pose_mat @ refinement_pose_mat1.transpose(0,2,1)
    theta_batch0 = geodesic_from_identity(relative_rotation0.reshape(-1,3,3)).reshape(relative_rotation0.shape[:-2])
    theta_batch1 = geodesic_from_identity(relative_rotation1.reshape(-1,3,3)).reshape(relative_rotation1.shape[:-2])
    min_symmetry0_idx = np.argmin(theta_batch0, axis=0)
    min_symmetry1_idx = np.argmin(theta_batch1, axis=0)
    min_batch = np.minimum(theta_batch0[min_symmetry0_idx, np.arange(len(min_symmetry0_idx))], 
                        theta_batch1[min_symmetry1_idx, np.arange(len(min_symmetry1_idx))])
    best_symmetry_relative_rotation = relative_rotation0[min_symmetry0_idx, np.arange(len(min_symmetry0_idx))]
    # mean of close rotations (avoid outliers, 1st percentile)
    d = geodesic_from_identity(best_symmetry_relative_rotation)
    argsort_idx = np.argsort(d)
    
    close_relative_rotation = best_symmetry_relative_rotation[argsort_idx[:max(1, int(len(argsort_idx)*percentage))]]
    # project to SO(3)
    mean_relative_rotation = np.mean(close_relative_rotation, axis=0)
    U, _, Vt = np.linalg.svd(mean_relative_rotation)
    G = U @ Vt
    if np.linalg.det(G) < 0:
        U[:, -1] *= -1
        G = U @ Vt
    global_residual_deg = np.rad2deg(geodesic_from_identity(G))
    print('global residual rotation error (degrees):', global_residual_deg)
    return G, global_residual_deg

def compute_geodesic_deg(initial_pose_mat, refinement_pose_mat0, refinement_pose_mat1, R_sym, G):
    n_sym = R_sym.shape[0]
    relative_rotation0 = np.zeros((n_sym, len(initial_pose_mat), 3, 3))
    relative_rotation1 = np.zeros((n_sym, len(initial_pose_mat), 3, 3))
    # apply this global rotation
    for idx in range(n_sym):
        relative_rotation0[idx] = R_sym[idx].reshape(1,3,3) @ initial_pose_mat @ (G @ refinement_pose_mat0).transpose(0,2,1)
        relative_rotation1[idx] = R_sym[idx].reshape(1,3,3) @ initial_pose_mat @ (G @ refinement_pose_mat1).transpose(0,2,1)

    theta_batch0 = geodesic_from_identity(relative_rotation0.reshape(-1,3,3)).reshape(relative_rotation0.shape[:-2])
    theta_batch1 = geodesic_from_identity(relative_rotation1.reshape(-1,3,3)).reshape(relative_rotation1.shape[:-2])
    min_symmetry0_idx = np.argmin(theta_batch0, axis=0)
    min_symmetry1_idx = np.argmin(theta_batch1, axis=0)
    min_batch = np.minimum(theta_batch0[min_symmetry0_idx, np.arange(len(min_symmetry0_idx))], 
                        theta_batch1[min_symmetry1_idx, np.arange(len(min_symmetry1_idx))])
    geodesic_deg = np.rad2deg(min_batch)
    return geodesic_deg

def pose_metric(initial_pose_mat, refinement_pose_mat0, refinement_pose_mat1, R_sym):
    rotation_global = np.eye(3)
    n_global_iterations = 40
    rotation_global = np.eye(3)
    for _ in range(n_global_iterations):
        G, _ = compute_global_rotation(initial_pose_mat, rotation_global @ refinement_pose_mat0, rotation_global @ refinement_pose_mat1, R_sym)
        geodesic_deg = compute_geodesic_deg(initial_pose_mat,  rotation_global @ refinement_pose_mat0,  rotation_global @ refinement_pose_mat1, R_sym, G)
        rotation_global = G @ rotation_global

    global_residual_deg = np.rad2deg(geodesic_from_identity(rotation_global))
    return geodesic_deg, rotation_global, global_residual_deg


def final_scores(df_submissions, metrics_cfg):
    """
    metrics_cfg contains:
        lower_better
        higher_better
        rotation_metrics
        state_metrics
    """

    df = df_submissions.copy()

    lower_better = metrics_cfg.lower_better
    higher_better = metrics_cfg.higher_better

    # -------------------------------------------------
    # 1. metric-level ranking (respect direction)
    # -------------------------------------------------
    rank_parts = []

    if lower_better:
        rank_parts.append(
            df[lower_better].rank(method="average", ascending=True)
        )

    if higher_better:
        rank_parts.append(
            df[higher_better].rank(method="average", ascending=False)
        )

    metric_ranks = pd.concat(rank_parts, axis=1)

    # -------------------------------------------------
    # 2. category scores (mean metric rank)
    # -------------------------------------------------
    rot_cols = metrics_cfg.rotation_metrics
    state_cols = metrics_cfg.state_metrics

    df["rot_score"] = metric_ranks[rot_cols].mean(axis=1)
    df["state_score"] = metric_ranks[state_cols].mean(axis=1)

    # -------------------------------------------------
    # 3. category ranks (1..N)
    # -------------------------------------------------
    df["rot_rank"] = df["rot_score"].rank(
        method="average", ascending=True
    )

    df["state_rank"] = df["state_score"].rank(
        method="average", ascending=True
    )

    # -------------------------------------------------
    # 4. equal weighting (50/50)
    # -------------------------------------------------
    df["combined_rank_score"] = (
        0.5 * df["rot_rank"]
        + 0.5 * df["state_rank"]
    )

    # -------------------------------------------------
    # 5. final leaderboard rank
    # -------------------------------------------------
    df["final_rank"] = df["combined_rank_score"].rank(
        method="average", ascending=True
    )

    return df.sort_values("final_rank")


@dataclass
class MetricsConfig:
    lower_better: List[str]
    higher_better: List[str]
    rotation_metrics: List[str]
    state_metrics: List[str]


