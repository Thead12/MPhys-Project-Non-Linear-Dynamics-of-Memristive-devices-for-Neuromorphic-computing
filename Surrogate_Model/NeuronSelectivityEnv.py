from typing import Optional
import time
import numpy as np
import gymnasium as gym
import sys
import os

from SurrogatePrediction import SurrogatePredictor

class NeuronSelectivityEnv(gym.Env):

    def __init__(self, f_target = 0.064, V_DC_0 = 6.0, R_ext_0 = 500.0):
        super(NeuronSelectivityEnv, self).__init__()

        self.surrogate = SurrogatePredictor()

        self.f_target = f_target
        self.V_DC_0 = V_DC_0
        self.R_ext_0 = R_ext_0

        # self.V_AC = 3.0 
        # self.target_freq = 0.064 
        # self.frequency_range = 0.02
        # self.frequency_step = 0.005
        # self.warmUp_cycles = 5.0
        # self.num_cycles = 10.0
        # self.depth_checking_range = 0.003
        # self.dt = 1e-3

        self.V_AC = 3.0 
        self.target_freq = 0.064 
        self.frequency_range = 0.04
        self.frequency_step = 0.002
        self.warmUp_cycles = 10.0
        self.num_cycles = 10.0
        self.depth_checking_range = 0.003
        self.dt = 1e-3


        # state bounds of [V_DC, R_ext]
        self.state_bounds_min = [3.5, 150.0]
        self.state_bounds_max = [10.0, 1000.0]
        # Observables: [width, CV2_at_target_frequency, centre_of_locking_frequency, depth, V_DC, R_ext]
        obs_low = np.array([0.03, 0, 0.01, 0, *self.state_bounds_min], dtype=np.float32)
        obs_high = np.array([0.107, 2, 0.1, 1.0, *self.state_bounds_max], dtype=np.float32)
        self.observation_space = gym.spaces.Box(low=obs_low, high=obs_high, dtype=np.float32)

        # Action: deltas for [V_DC, R_ext]
        self.action_space = gym.spaces.Box(low=-0.1, high=0.1, shape=(2,), dtype=np.float32)

        # Internal agent state = [V_DC, R_ext]
        self._agent_state = np.array([V_DC_0, R_ext_0], dtype=np.float32)
        self._last_results = None  # Cache last simulation results
        self._step_count = 0
        self._max_episode_steps = 100

    def _get_obs(self):
        import time
        start_time = time.time()
        V_DC, R_ext = np.float32(self._agent_state[0]).item(), np.float32(self._agent_state[1]).item()
        
        metrics = self.surrogate.predict([V_DC, R_ext])
        width = metrics[0]
        CV2_at_target_freq = metrics[1]
        centre_of_phase_locking = metrics[2]  # centre of phase locking region
        depth = metrics[3]

        obs = np.array([width, CV2_at_target_freq, centre_of_phase_locking, depth, *self._agent_state], dtype=np.float32)

        return obs 
    
    def _get_info(self):
        """Compute auxiliary information for debugging.

        Returns:
            dict: Info with current state and simulation results
        """
        res = self._last_results if self._last_results else {}
        return {
            "V_DC": self._agent_state[0],
            "R_ext": self._agent_state[1],
            "f_center": res.get("centre_of_phase_locking", 0),
            "depth": res.get("depth", 0),
            "CV2_at_target": res.get("CV2_at_target_freq", 0),
        }
    
    def reset(self, seed=None, options: Optional[dict] = None):
        super().reset(seed=seed)

        # 1. Set agents initial state
        if options:
            self.V_DC_0 = options.get("V_DC_0", self.V_DC_0)
            self.R_ext_0 = options.get("R_ext_0", self.R_ext_0)

        
        # Adds slight randomness to starting positions
        self.V_DC_0 += self.np_random.uniform(-0.3, 0.3)
        self.R_ext_0 += self.np_random.uniform(-25.0, 25.0)

        self._agent_state = np.array([self.V_DC_0, self.R_ext_0], dtype=np.float32)
        self._last_results = None  # Clear cached results
        self._step_count = 0

        return self._get_obs(), {}

    def step(self, action):
        self._step_count += 1
        self._agent_state += action
        self._agent_state = np.clip(self._agent_state, self.state_bounds_min, self.state_bounds_max)

        # make an observation
        obs_start = time.time()
        obs = self._get_obs()
        obs_elapsed = time.time() - obs_start

        width, CV2_at_target_frequency, centre_of_phase_locking, depth = obs[0:4]
        freq_error = abs(self.f_target - centre_of_phase_locking)

        # calculate reward:
        # 1. Frequency centre close to target
        # 2. Deep selectivity well (high depth)
        # 3. Low CV2 at target frequency (good regularity)
        # 4. Minimise width (encourage sharp selectivity)
        w1 = 1.0  # weight for width
        w2 = 1.0  # weight for CV2_at_target_frequency
        w3 = 1.0  # weight for freq_error
        w4 = 1.0  # weight for depth

        reward = - (w1 * width) + (w2 * CV2_at_target_frequency) - (w3 * freq_error) + (w4 * depth)

        truncated = self._step_count >= self._max_episode_steps
        return obs, reward, False, truncated, self._get_info()
