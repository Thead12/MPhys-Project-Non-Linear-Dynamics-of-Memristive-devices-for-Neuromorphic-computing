import gymnasium as gym

from stable_baselines3 import PPO
from stable_baselines3.common.env_util import make_vec_env
from stable_baselines3.common.vec_env import SubprocVecEnv

from NeuronSelectivityEnv import NeuronSelectivityEnv

# neuron_env = NeuronSelectivityEnv()
def train(n_envs=4,n_steps=32,time_steps=1):
    total_timesteps=n_envs*n_steps*time_steps

    # Parallel environments
    vec_env = make_vec_env(NeuronSelectivityEnv,
                           n_envs=n_envs,
                           vec_env_cls=SubprocVecEnv)
    
    model = PPO("MlpPolicy", vec_env, n_steps=n_steps, verbose=0, tensorboard_log="./ppo_neuron_logs/")
    model.learn(total_timesteps=total_timesteps, progress_bar=True)
    model.save("ppo_neuron")
    
    vec_env.close()
    del model # remove to demonstrate saving and loading
    print("===== Training Complete =====")

import numpy as np
import csv
from tqdm import tqdm
def evaluate(V_DC_0=6.0):
    print("===== Starting Evaluation =====")

    env = NeuronSelectivityEnv(V_DC_0=V_DC_0)
    model = PPO.load("ppo_neuron")
    obs, info = env.reset()

    terminated = False
    truncated = False

    max_steps = 100
    step = 0
    pbar = tqdm(total=max_steps, desc="Evaluating Agent")
    trajectory_data = []
    while not (terminated or truncated):
        action, _states = model.predict(obs, deterministic=True)
        obs, reward, terminated, truncated, info = env.step(action)
        trajectory = np.insert(obs, 0, reward)
        trajectory = np.insert(trajectory, 0, step)
        trajectory_data.append(trajectory)

        pbar.update(1)
        step+=1
        pbar.set_postfix({"V_DC": f"{obs[4]:.2f}","R_ext": f"{obs[5]:.2f}", "Reward": f"{reward:.2f}"})
    pbar.close()
    env.close()

    log_file = open("final_trajectory.csv", "w", newline='')
    writer = csv.writer(log_file)
    writer.writerow(["step","reward","width","CV2_at_target_frequency","centre_of_phase_lock","depth","V_DC","R_ext"])
    writer.writerows(trajectory_data)

    print("Evaluation finished.")

if __name__ == '__main__':
    train(n_envs=8,n_steps=64,time_steps=100)
    evaluate(V_DC_0=4)
