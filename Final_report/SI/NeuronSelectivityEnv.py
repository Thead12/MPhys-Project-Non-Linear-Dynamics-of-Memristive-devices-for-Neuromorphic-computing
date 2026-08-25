import numpy as np
import gymnasium as gym
from gymnasium import spaces
from typing import Optional, Dict, Any, Tuple
from SurrogatePrediction import SurrogatePredictor

N_OUTPUTS = 5
SIZE = 64

# ----------- Selectivity environment for the PPO agent to interact with ------------------------

def calculate_reward(f_target, metrics):

    # calculate reward
    # 1 frequency centre close to target
    # 2 deep selectivity well (high depth)
    # 3 low CV2 at target frequency (regularity)
    # 4 minimise width (encourage sharp selectivity)
    width, cv1, cv2, f_centre, depth = metrics

    freq_error = abs(f_target - f_centre)
    
    k_f = 15.0
    k_w = 6 
    sigma_f = 0.2
    delta = 0.05
    f_penalty = k_f * (np.sqrt(freq_error**2 + delta**2) - delta)
    
    gate = np.exp(-(freq_error**2)/(2*sigma_f**2))
    w_penalty = k_w*width*gate
    return -(f_penalty + w_penalty)


class NeuronSelectivityEnv(gym.Env):
    def __init__(self, V_DC_0=0.5, R_ext_0=0.5, max_steps = 100, f_target = 0.571):
        super(NeuronSelectivityEnv, self).__init__()
        
        # setup
        self.surrogate = SurrogatePredictor(n_outputs=N_OUTPUTS, size=SIZE)
        self.f_target = f_target
        self.V_DC_0=V_DC_0
        self.R_ext_0=R_ext_0
        self._max_episode_steps=max_steps

        # state bounds of [V_DC, R_ext]
        self.state_bounds_min = np.array([0.0, 0.0], dtype=np.float32)
        self.state_bounds_max = np.array([1.0, 1.0], dtype=np.float32)

        # observation space: [V_DC, R_ext, width, CV1, CV2, f_centre, depth]
        obs_low = np.array([0.0]*7, dtype=np.float32)
        obs_high = np.array([1.0]*7, dtype=np.float32)
        self.observation_space = spaces.Box(low=obs_low, high=obs_high, dtype=np.float32)

        # action space: deltas for [V_DC, R_ext]
        self.action_space = gym.spaces.Box(low=-0.01, high=0.01, shape=(2,), dtype=np.float32)

        # internal agent state = [V_DC, R_ext]
        self._agent_state = np.array([self.V_DC_0, self.R_ext_0], dtype=np.float32)
        self._step_count = 0
        self._last_metrics = np.zeros(5, dtype=np.float32)

    def _get_obs(self):
        V_DC, R_ext = self._agent_state
        self._last_metrics = self.surrogate.predict_without_scaling([V_DC, R_ext])

        # combine action and metrics into single observation 
        obs = np.concatenate(([V_DC, R_ext], self._last_metrics)).astype(np.float32)
        return obs 
    
    def _get_info(self) -> Dict[str, Any]:
        return {
            "V_DC": self._agent_state[0],
            "R_ext": self._agent_state[1],
            "metrics": self._last_metrics,
            "step": self._step_count
        }
    
    def reset(self, seed: Optional[int] = None, options: Optional[dict] = None) -> Tuple[np.array, Dict]:
        super().reset(seed=seed)
        self._step_count = 0
        # set agents initial state
        if options and options.get("mode") == "eval":
            self._agent_state = np.array([self.V_DC_0, self.R_ext_0], dtype=np.float32)
        else:
            self._agent_state = self.np_random.uniform(self.state_bounds_min, self.state_bounds_max).astype(np.float32)

        return self._get_obs(), self._get_info()

    def step(self, action):
        self._step_count += 1 

        # update state and clip
        self._agent_state = np.clip(self._agent_state + action, self.state_bounds_min, self.state_bounds_max)

        # observe new state
        obs = self._get_obs()

        # unpack metrics
        reward = calculate_reward(self.f_target, obs[2:])

        terminated = False
        truncated = self._step_count >= self._max_episode_steps

        return obs, float(reward), terminated, truncated, self._get_info()
