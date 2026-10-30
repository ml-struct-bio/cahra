import time
from omegaconf import OmegaConf, DictConfig
import hydra
import starfile
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial.transform import Rotation as R
from dataclasses import dataclass
import pickle

from scoring_utils import (
    pose_metric,
    MetricsConfig,
    compute_state_metrics,
    final_scores,
)
from rotation_metrics import compute_symmetry_matrices


def ground_truth(fname):
    gt = starfile.read(fname)["particles"]
    s_gt = gt["rlnstate"]
    rotation_gt_mat = R.from_euler(
        "ZYZ",
        angles=gt[["rlnAngleRot", "rlnAngleTilt", "rlnAnglePsi"]].values,
        degrees=True,
    ).as_matrix()
    return gt, s_gt, rotation_gt_mat


@dataclass
class Config:
    precomputed_fname: str
    metrics: MetricsConfig


@hydra.main(version_base=None, config_path=".", config_name="scoring")
def main(config: DictConfig):
    timestr = time.strftime("%Y%m%d")

    OmegaConf.to_object(config)
    print("Config:", config)

    if config["precomputed_fname"] is not None:
        df_submissions = pd.read_csv(config["precomputed_fname"])
    else:
        with open(config["pre_wrangled_submissions_fname"], "rb") as f:
            submissions = pickle.load(f)

        gt, s_gt, rotation_gt_mat = ground_truth(
            config["ground_truth_star_fname"]
        )
        R_sym = compute_symmetry_matrices(n_sym=11)
        print("Group action matrices:", R_sym)

        all_metrics_list = []
        global_rotations = []
        per_particle_df_dict = {}
        for submitter in submissions.keys():
            for col in ["rlnAngleRot", "rlnAngleTilt", "rlnAnglePsi"]:
                print(gt.columns)
                print(submissions[submitter]["starfile_particles"].columns)
                plt.scatter(
                    x=gt[col],
                    y=submissions[submitter]["starfile_particles"][col],
                    alpha=0.03,
                    label=col,
                )
                plt.legend(framealpha=1.0, markerscale=1, scatterpoints=100)
                plt.grid(True, alpha=0.3)
                plt.savefig(
                    f"{config['output_dir']}/plots/pose_marginals_{submitter}_{col}.png"
                )
                plt.clf()

            print(f"Metrics for {submitter}:")
            s_pred = submissions[submitter]["s_pred"]
            rotation_sub_mat = submissions[submitter]["rotation_mat"]
            state_metrics = compute_state_metrics(s_pred, s_gt)
            state_metrics["submitter"] = submitter
            print(state_metrics)
            # state_metrics_list.append(state_metrics)

            geodesic_deg, global_rotation, global_residual_deg = pose_metric(
                rotation_gt_mat, rotation_sub_mat, rotation_sub_mat, R_sym
            )
            plt.hist(geodesic_deg, bins=50)
            plt.xlabel("Geodesic Angle (degrees)")
            plt.ylabel("Count")
            plt.title(f"Rotation Residuals for {submitter}")
            plt.grid(False)
            plt.savefig(
                f"{config['output_dir']}/plots/pose_residuals_{submitter}.png"
            )
            plt.clf()
            pose_metrics = {
                "rot_mean": float(np.mean(geodesic_deg)),
                "rot_median_error": float(np.median(geodesic_deg)),
                "rot_10th_percentile_error": float(
                    np.percentile(geodesic_deg, 10)
                ),
                "rot_90th_percentile_error": float(
                    np.percentile(geodesic_deg, 90)
                ),
                "rot_global_residual": float(global_residual_deg),
            }
            print(pose_metrics)
            all_metrics = {**state_metrics, **pose_metrics}
            print(all_metrics)
            all_metrics_list.append(all_metrics)
            global_rotations.append(global_rotation)

            per_particle_df_dict[submitter] = pd.DataFrame(
                {
                    "submitter": submitter,
                    #   'particle_idx': gt['particle_idx'],
                    "geodesic_deg": geodesic_deg,
                    "state": s_gt,
                    "predicted_state": s_pred,
                    "snr": gt["rlnsnr"],
                }
            )
        # pickle
        with open(
            f"{config['output_dir']}/validation_per_particle_{timestr}.pkl",
            "wb",
        ) as f:
            pickle.dump(per_particle_df_dict, f)

        global_rotations_arr = np.array(global_rotations)
        print(global_rotations_arr)
        np.savez(
            f"{config['output_dir']}/validation_poses_{timestr}.npz",
            global_rotation=global_rotations_arr,
            submission=list(submissions.keys()),
        )

        fname = f"{config['output_dir']}/validation_{timestr}.csv"
        df_submissions = pd.DataFrame(all_metrics_list)
        df_submissions.to_csv(fname, index=False)

    df_submissions = final_scores(df_submissions, config["metrics"])
    df_submissions.to_csv(
        f"{config['output_dir']}/validation_final_{timestr}.csv", index=False
    )


if __name__ == "__main__":
    main()
