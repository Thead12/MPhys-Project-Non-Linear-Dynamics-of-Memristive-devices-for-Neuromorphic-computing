import os
import gymnasium as gym

from stable_baselines3 import PPO
from stable_baselines3.common.env_util import make_vec_env
from stable_baselines3.common.vec_env import SubprocVecEnv
from stable_baselines3.common.utils import LinearSchedule

from NeuronSelectivityEnv import NeuronSelectivityEnv

# --------This is the file where the PPO agent is trained 
#          and can be evaluated for a given initial (V_ext,R_ext) combo ---------------------

def train(n_envs=4, n_steps=32, time_steps=1, model_name="ppo_neuron"):
    total_timesteps = n_envs * n_steps * time_steps
    
    # parallel environments
    vec_env = make_vec_env(NeuronSelectivityEnv,
                           n_envs=n_envs,
                           vec_env_cls=SubprocVecEnv)
    
    model_path = f"{"ppo_neuron"}.zip"

    if os.path.exists(model_path):
        print(f"*** Loading existing model: {model_path} ***")
        model = PPO.load(model_path, env=vec_env)
        
    else:
        print("*** No existing model found. Creating a new one. ***")
        lr_schedule = LinearSchedule(5e-5, 1e-6, 1.0)
        model = PPO("MlpPolicy", vec_env, n_steps=n_steps, verbose=1,
                    learning_rate=lr_schedule,
                    ent_coef=0.001)

    print(f"Beginning training for {total_timesteps} steps...")
    model.learn(total_timesteps=total_timesteps, progress_bar=True, reset_num_timesteps=False)
    
    model.save(model_name)
    vec_env.close()
    del model
    print("===== Training Complete =====")

import numpy as np
import csv
from tqdm import tqdm
def evaluate(V_DC_0=None, R_ext_0=None,
             max_steps=100, model_path="ppo_neuron",
             output_file="final_trajectory.csv",
             convergence_threshold=0.1,
             convergence_window=10):
    print(f"===== Starting Evaluation at V_DC: {V_DC_0}, R_ext: {R_ext_0} =====")


    env = NeuronSelectivityEnv(V_DC_0=V_DC_0, R_ext_0=R_ext_0, max_steps=max_steps)
    model = PPO.load(model_path)
    
    obs, info = env.reset(options={"mode": "eval"})

    terminated = False
    truncated = False
    step = 0
    max_steps = env._max_episode_steps
    trajectory_data = []
    recent_rewards = []
   
    pbar = tqdm(total=max_steps, desc="Evaluating Agent")
    while not (terminated or truncated):
        action, _states = model.predict(obs, deterministic=True)

        obs, reward, terminated, truncated, info = env.step(action)

        row = np.concatenate(([step, reward, action[0], action[1]], obs))
        trajectory_data.append(row)

        # convergence check
        recent_rewards.append(reward)
        if len(recent_rewards) > convergence_window:
             recent_rewards.pop()
             reward_std = np.std(recent_rewards)
             print(reward_std)
             if reward_std < convergence_threshold:
                  print(f"Converged at step {step}")
                  break

        
        pbar.update(1)
        step+=1
        pbar.set_postfix({
            "V_DC": f"{obs[0]:.2f}",
            "R_ext": f"{obs[1]:.2f}",
            "Reward": f"{reward:.2f}"})
    pbar.close()
    env.close()

    with open(output_file, "w", newline='') as log_file:
              writer = csv.writer(log_file)
              writer.writerow([
                  "step","reward","dV_DC","dR_ext","V_DC","R_ext",
                  "width","cv1","cv2","f_centre","depth"
              ])
              writer.writerows(trajectory_data)

    print(f"===== Evaluation Complete. Results saved to {output_file} =====")



if __name__ == '__main__':
    train(n_envs=20,n_steps=128,time_steps=500)
    evaluate(V_DC_0=0.8, R_ext_0=0.8, max_steps=100)
