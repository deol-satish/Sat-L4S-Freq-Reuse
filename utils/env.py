import gymnasium
import numpy as np
import matlab.engine
from gymnasium.spaces import MultiDiscrete, Dict, Box
import logging
import json
import math
from datetime import datetime, timedelta, timezone
import time
import os



env_name = "CogSatEnv-v1"

 
class CogSatEnv(gymnasium.Env):
    """Gymnasium environment for MATLAB-based Cognitive Satellite Simulation"""
 
    def __init__(self, env_config=None, render_mode=None):
        super(CogSatEnv, self).__init__()
        if not hasattr(self, 'spec') or self.spec is None:
            self.spec = gymnasium.envs.registration.EnvSpec("CogSatEnv-v1")


        self.tIndex = 0
        self.mat_filename = "A2C"     


        self.no_channels = env_config.get('no_channels')
        self.no_leo_user = env_config.get('no_leo_user')
        self.no_geo_user = env_config.get('no_geo_user')
        self.saved_folder = env_config.get('saved_folder')

        # Logging statements
        logging.info("self.no_channels: %s", self.no_channels)
        logging.info("self.no_leo_user: %s", self.no_leo_user)
        logging.info("self.no_geo_user: %s", self.no_geo_user)
        logging.info("self.saved_folder: %s", self.saved_folder)

        self.episode_number = -1            
 
        # Start MATLAB engine and set path
        self.eng = matlab.engine.start_matlab()
        self.eng.cd(self.eng.pwd(), nargout=0)  # ensure working directory is set

        self.eng.addpath(r'./downlink_matlab_code', nargout=0)

        # Initialize the MATLAB scenario and Save Baseline

        # Define MATLAB code as a string
        matlab_code = f"""
        clear; clc; close all hidden;
        rng(42); % For reproducibility
        P01_Parameters
        NumLeoUser = {self.no_leo_user};
        NumGeoUser = {self.no_geo_user};
        numChannels =  {self.no_channels};
        P02_GStations
        fprintf('Creating satellite scenario...\\n');
        sc = satelliteScenario(startTime, stopTime, sampleTime);
        P03_GeometricSimulatoion
        P04_RxSimulation
        P05_ChannelAllocation

        T = length(ts);
        SINR = NaN(NumGS, T);
        Intf = NaN(NumGS, T);
        Thrpt = NaN(NumGS, T);
        SINR_mW_dict = NaN(NumGS, T);
        Intf_mW_dict = NaN(NumGS, T);
        Noise_mW = 10^(ThermalNoisedBm / 10);
        SE = zeros(1, T);
        FreqAlloc = NaN(NumGS, length(ts));
        t = 1;
        """

        # Evaluate the MATLAB code in the engine
        self.eng.eval(matlab_code, nargout=0)

        # Log the values using logging.info
        logging.info("self.eng.workspace['numChannels']: %s", self.eng.workspace['numChannels'])
        logging.info("self.eng.workspace['NumLeoUser']: %s", self.eng.workspace['NumLeoUser'])
        logging.info("self.eng.workspace['NumGeoUser']: %s", self.eng.workspace['NumGeoUser'])

        self.eng.eval("P06_Intf_Eval", nargout=0)
        
        logging.info("=== Baseline SINR dB === %s", np.array(self.eng.workspace['SINR'])[:, self.tIndex])

        self.eng.eval("resetScenario", nargout=0)
        self.save_npy_data('Baseline')       

        filename = f"{self.saved_folder}/Baseline_workspace_Saved_data.mat"
        save_cmd = f"save('{filename}')"
        self.eng.eval(save_cmd, nargout=0)

        Lat = np.array(self.eng.workspace['Lat'])
        np.save(f'{self.saved_folder}/Baseline_Lat.npy', Lat)


          


        
        self.timelength = self.eng.eval("length(ts)", nargout=1)
        self.NumLeoUser = int(self.eng.workspace['NumLeoUser'])
        self.NumGeoUser = int(self.eng.workspace['NumGeoUser'])
        
        self.LeoChannels = int(self.eng.workspace['numChannels'])
        self.GeoChannels = int(self.eng.workspace['numChannels'])

        self.reward = 0


        self.cur_obs = {
            "utc_time": np.array([0], dtype=np.int64),
            "freq_lgs_leo": np.random.uniform(1.0, self.LeoChannels, size=(self.NumLeoUser,)).astype(np.int64),
            "freq_ggs_geo": np.random.uniform(1.0, self.GeoChannels, size=(self.NumGeoUser,)).astype(np.int64),
            "leo_pos": np.random.uniform(0, 20, size=(self.NumLeoUser*2,)).astype(np.float32),
            
        }
      
        
 
        # Define action and observation space
        self.action_space = gymnasium.spaces.MultiDiscrete([self.LeoChannels] * self.NumLeoUser)  # Select a channel index for each LEO user


        # Observation space structure
        self.observation_space = Dict({
            "utc_time": Box(low=-np.inf, high=np.inf, shape=(1,), dtype=np.int64),
            "freq_lgs_leo": Box(low=1, high=self.LeoChannels+1, shape=(self.NumLeoUser,), dtype=np.int64),
            "freq_ggs_geo": Box(low=1, high=self.GeoChannels+1, shape=(self.NumGeoUser,), dtype=np.int64),
            "leo_pos": Box(low=-np.inf, high=np.inf, shape=(self.NumLeoUser*2,), dtype=np.float32),
        })

        self.ep_start_time = time.time()

    def save_npy_data(self,extra_tag="Original"):
        SINR = np.array(self.eng.workspace['SINR'])
        Intf = np.array(self.eng.workspace['Intf'])
        SINR_mW_dict = np.array(self.eng.workspace['SINR_mW_dict'])
        Intf_mW_dict = np.array(self.eng.workspace['Intf_mW_dict'])
        Thrpt = np.array(self.eng.workspace['Thrpt'])
        SE = np.array(self.eng.workspace['SE'])
        FreqAlloc = np.array(self.eng.workspace['FreqAlloc'])
        berQPSK = np.array(self.eng.workspace['berQPSK'])
        berMQAM = np.array(self.eng.workspace['berMQAM'])
        

        np.save(f'{self.saved_folder}/{extra_tag}_SINR.npy', SINR)
        np.save(f'{self.saved_folder}/{extra_tag}_Intf.npy', Intf)
        np.save(f'{self.saved_folder}/{extra_tag}_SINR_mW_dict.npy', SINR_mW_dict)
        np.save(f'{self.saved_folder}/{extra_tag}_Intf_mW_dict.npy', Intf_mW_dict)
        np.save(f'{self.saved_folder}/{extra_tag}_Thrpt.npy', Thrpt)
        np.save(f'{self.saved_folder}/{extra_tag}_SE.npy', SE)
        np.save(f'{self.saved_folder}/{extra_tag}_FreqAlloc.npy', FreqAlloc)
        np.save(f'{self.saved_folder}/{extra_tag}_berQPSK.npy', berQPSK)
        np.save(f'{self.saved_folder}/{extra_tag}_berMQAM.npy', berMQAM)
    
    
    def get_matlab_ts(self):
        """
        Get the MATLAB timestamp as a list of strings.
        """
        ts_str = self.eng.eval("cellstr(datestr(ts, 'yyyy-mm-ddTHH:MM:SS'))", nargout=1)
        python_datetimes = [datetime.fromisoformat(s) for s in ts_str]
        timestamps = [dt.timestamp() for dt in python_datetimes]
        return timestamps
    

    


    def get_state_from_matlab(self):
        # Log cur_state_from_matlab
        logging.info("=== Current State ===")
        # logging.info(json.dumps(cur_state_from_matlab, indent=2))
        """Reset the environment and initialize the buffer."""

        self.ts = self.get_matlab_ts()

        self.FreqAlloc = np.array(self.eng.workspace['FreqAlloc'])
        self.LEOFreqAlloc = self.FreqAlloc[:self.NumLeoUser,:]
        self.GEOFreqAlloc = self.FreqAlloc[self.NumLeoUser:,:]

        self.cur_obs["utc_time"] = np.array([self.ts[self.tIndex]], dtype=np.int64)
        self.cur_obs["freq_lgs_leo"] = np.array(self.LEOFreqAlloc[:,self.tIndex], dtype=np.int64)
        self.cur_obs["freq_ggs_geo"] = np.array(self.GEOFreqAlloc[:,self.tIndex], dtype=np.int64)
        leo_loc = np.array(self.eng.workspace['LEO_LOC'])
        self.cur_obs["leo_pos"] = np.array(leo_loc[0:self.NumLeoUser,self.tIndex].flatten(), dtype=np.float32)

        # Log freq_lgs_leo
        logging.info("freq_lgs_leo: %s", self.cur_obs["freq_lgs_leo"].tolist())
        logging.info("freq_ggs_geo: %s", self.cur_obs["freq_ggs_geo"].tolist())

        # (Optional) Validate against observation_space
        assert self.observation_space.contains(self.cur_obs), "cur_obs doesn't match the observation space!"

        return self.cur_obs
    

 
    def step(self, action):
        """
        Apply action and return (observation, reward, terminated, truncated, info)
        """
        terminated = False
        truncated = False
        self.eng.workspace['t'] = int(self.tIndex) + 1
        start_time = time.time()

        # Action start from 0 and ends before self.LeoChannels, that means it is not included
        # For example, if self.LeoChannels is 5, action can be 0, 1, 2, 3, or 4.
        # This is because MATLAB uses 1-based indexing, so we need to convert it to 0-based indexing for Python.

        logging.info("")
        logging.info("******************************************************************************************")
        logging.info("=== Step Started ===")
        t = np.array(self.eng.workspace['t'])
        logging.info("===step: t from matlab: %s ===",t)

        action = action + 1 
        #print("After +1 Operation: Action taken: ", action)
        logging.info("=== After +1 Operation: Action Taken === %s", action)

        # Access the variable from MATLAB workspace
        # Convert MATLAB array to NumPy array
        channel_list_leo = np.array(self.eng.workspace['ChannelListLeo'])
        Serv_idxLEO = np.array(self.eng.workspace['Serv_idxLEO'])

        # FreqAlloc = np.array(self.eng.workspace['FreqAlloc'])[:,self.tIndex]
        # logging.info("step Before Assignment: FreqAlloc: %s",FreqAlloc)

        # Store original values for comparison
        before_vals = []
        for i in range(self.NumLeoUser):
            sat = Serv_idxLEO[i, self.tIndex].astype(int) - 1
            before_vals.append(channel_list_leo[i, sat, self.tIndex])

        # Now do the assignment 
        for i in range(self.NumLeoUser):
            sat = Serv_idxLEO[i, self.tIndex].astype(int) - 1
            channel_list_leo[i, sat, self.tIndex] = action[i]

        # Compare
        for i in range(self.NumLeoUser):
            sat = Serv_idxLEO[i, self.tIndex].astype(int) - 1
            logging.info(f"User {i}: Before = {before_vals[i]}, After = {channel_list_leo[i, sat, self.tIndex]}")
        
        self.eng.workspace['ChannelListLeo'] = matlab.double(channel_list_leo)

        self.eng.eval("stepScenario", nargout=0)





        next_observation = self.get_state_from_matlab()  

        FreqAlloc = np.array(self.eng.workspace['FreqAlloc'])[:,self.tIndex]
        # logging.info("step AfterD Assignment: FreqAlloc: %s",FreqAlloc)
    

        # Split into LEO and GEO
        leo_vals = FreqAlloc[:self.NumLeoUser]
        geo_vals = FreqAlloc[self.NumLeoUser:]

        # Find common elements
        common_vals = np.intersect1d(leo_vals, geo_vals)

        # Count
        num_common = len(common_vals)

        # print("Common values:", common_vals)

        
        #print("Next Observation: ", next_observation)
        logging.info("=== Next Observation === %s", next_observation)

        logging.info("===  After Action: FreqAlloc === %s",FreqAlloc)
        logging.info("=== common values === %s", num_common)

        SINR_dB = np.array(self.eng.workspace['SINR'])
        SINR_mW = np.array(self.eng.workspace['SINR_mW_dict'])
        Intf_dB = np.array(self.eng.workspace['Intf'])
        Intf_mw = np.array(self.eng.workspace['Intf_mW_dict'])

        logging.info("=== SINR dB === %s", SINR_dB)
        logging.info("=== Intf dB === %s", Intf_dB)
  
         


        # SINR = np.array(self.eng.workspace['SINR_mW_dict'])
        # Intf = np.array(self.eng.workspace['Intf_mW_dict'])
        Thrpt = np.array(self.eng.workspace['Thrpt'])/(1024*1024)


        

        # Compute reward: sum of SINR across all users at current time step
        mw_Interference_to_leo_users = Intf_mw[:self.NumLeoUser, self.tIndex]
        mw_Interference_to_geo_users = Intf_mw[self.NumLeoUser:self.NumGeoUser+ self.NumLeoUser, self.tIndex]
        db_Interference_to_leo_users = Intf_dB[:self.NumLeoUser, self.tIndex]
        db_Interference_to_geo_users = Intf_dB[self.NumLeoUser:self.NumGeoUser+ self.NumLeoUser, self.tIndex]
        mw_SINR_of_LEO_users = SINR_mW[:self.NumLeoUser, self.tIndex]
        mw_SINR_of_GEO_users = SINR_mW[self.NumLeoUser:self.NumGeoUser+ self.NumLeoUser, self.tIndex]
        db_SINR_of_LEO_users = SINR_dB[:self.NumLeoUser, self.tIndex]
        db_SINR_of_GEO_users = SINR_dB[self.NumLeoUser:self.NumGeoUser+ self.NumLeoUser, self.tIndex]
        Thrpt_of_LEO_users = Thrpt[:self.NumLeoUser, self.tIndex]


        logging.info("=== LEO Thrpt === %s", Thrpt_of_LEO_users)



 

        # self.reward = np.sum(np.log10(mw_SINR_of_LEO_users)) -  num_common

        # self.reward = -1 *( np.median(mw_Interference_to_leo_users) + 0.25 * np.median(mw_Interference_to_geo_users))
        # self.reward = -1 *( np.median(db_Interference_to_leo_users) + 0.25 * np.median(db_Interference_to_geo_users))

        self.reward = np.sum(np.log10(mw_SINR_of_LEO_users)) + 0.25*(np.sum(np.log10(mw_SINR_of_GEO_users)))
        # self.reward = np.sum(np.log10(db_SINR_of_LEO_users)) + 0.25*(np.sum(np.log10(db_SINR_of_GEO_users)))
        # self.reward = np.sum(np.log10(mw_SINR_of_LEO_users)) + (np.sum(np.log10(mw_SINR_of_GEO_users)))
        # self.reward = np.sum(np.log10(db_SINR_of_LEO_users)) + (np.sum(np.log10(db_SINR_of_GEO_users)))



        logging.info("=== Reward === %s", self.reward)

        self.tIndex += 1
        
        if self.tIndex >= self.timelength:
            terminated = True
            logging.info("=== Episode finished after %s timesteps ===", self.tIndex)
            self.save_npy_data(f'Episode_{self.episode_number}')
            # Define MATLAB code as a string
            matlab_code = f"""
            T = length(ts);
            SINR = NaN(NumGS, T);
            Intf = NaN(NumGS, T);
            Thrpt = NaN(NumGS, T);
            SINR_mW_dict = NaN(NumGS, T);
            Intf_mW_dict = NaN(NumGS, T);
            Noise_mW = 10^(ThermalNoisedBm / 10);
            SE = zeros(1, T);
            FreqAlloc = NaN(NumGS, length(ts));
            t = 1;
            """

            # Evaluate the MATLAB code in the engine
            self.eng.eval(matlab_code, nargout=0)

            # Log the values using logging.info
            logging.info("self.eng.workspace['numChannels']: %s", self.eng.workspace['numChannels'])
            logging.info("self.eng.workspace['NumLeoUser']: %s", self.eng.workspace['NumLeoUser'])
            logging.info("self.eng.workspace['NumGeoUser']: %s", self.eng.workspace['NumGeoUser'])
            self.tIndex = 0
            logging.info("=== Baseline SINR dB === %s", np.array(self.eng.workspace['SINR'])[:, self.tIndex])

            self.eng.eval("resetScenario", nargout=0)
            # filename = f"{self.saved_folder}/{self.mat_filename}_{self.episode_number}_workspace_Saved_data_Close.mat"
            # save_cmd = f"save('{filename}')"
            # self.eng.eval(save_cmd, nargout=0)
            # self.eng.eval("P08_SaveData", nargout=0)

        info = {}

        end_time = time.time()
        step_duration = end_time - start_time
        logging.info("=== Time taken for timestep: %.4f seconds ===", step_duration)
        logging.info("")
        logging.info("******************************************************************************************")

 
        return next_observation, self.reward, terminated, truncated, info
 
    def reset(self, *, seed=None, options=None):
        self.episode_number = self.episode_number + 1
        #print("Resetting environment for episode: ", self.episode_number)
        logging.info("=== Resetting Environment for Episode %s ===", self.episode_number)
        super().reset(seed=seed) 
        # Reset the scenario

        self.eng.eval("resetScenario", nargout=0)
        self.eng.eval("stepScenario", nargout=0)
        self.eng.eval("resetScenario", nargout=0)

        self.tIndex = 0
        self.done = 0
        
        self.ep_end_time = time.time()
        self.ep_step_duration = self.ep_end_time - self.ep_start_time
        logging.info("=== Episode Index: {} Time taken for timestep: {:.4f} seconds ===".format(self.episode_number, self.ep_step_duration))

        self.ep_start_time = time.time()

        observation = self.get_state_from_matlab()
        #print("++++===== ENV RESET+++===")
 
        return observation, {}

    def save_env_state(self):
        """
        Saves the current environment state, including MATLAB workspace variables.
        """
        return {
            "tIndex": self.tIndex,
            "ChannelListLeo": np.array(self.eng.workspace["ChannelListLeo"]),
        }

    def restore_env_state(self, state):
        """
        Restores the environment state from the saved state dictionary.
        """
        # Restore tIndex
        self.tIndex = state["tIndex"]

        # Restore MATLAB workspace variables
        self.eng.workspace["ChannelListLeo"] = matlab.double(state["ChannelListLeo"].tolist())
        self.eng.workspace["FreqAlloc"] = matlab.double(state["FreqAlloc"].tolist())

        # Sync observation state
        self.cur_obs["utc_time"] = np.array([self.ts[self.tIndex]], dtype=np.int64)
        self.cur_obs["freq_lgs_leo"] = np.array(state["FreqAlloc"][:self.NumLeoUser, self.tIndex], dtype=np.int64)
        self.cur_obs["freq_ggs_geo"] = np.array(state["FreqAlloc"][self.NumLeoUser:self.NumLeoUser+ self.NumGeoUser, self.tIndex], dtype=np.int64)

        logging.info("=== Environment Restored to tIndex %s ===", self.tIndex)

 
    def render(self):
        pass
 
    def close(self):
        #print("Saving MATLAB Data.")
        logging.info("=== Saving MATLAB Data ===")
        # self.eng.eval("P08_SaveData", nargout=0)
        filename = f"{self.saved_folder}/{self.mat_filename}_workspace_Saved_data_Close.mat"
        save_cmd = f"save('{filename}')"
        self.eng.eval(save_cmd, nargout=0)
        self.eng.quit()
    