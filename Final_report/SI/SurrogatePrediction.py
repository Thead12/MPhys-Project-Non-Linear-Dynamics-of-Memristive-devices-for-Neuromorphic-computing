import torch 
import joblib
import numpy as np
from SurrogateModel import SurrogateModel

# --------- A wrapper class that lets you use the surrogate model -------
# can also be used to apply and reverse minmax scaling
class SurrogatePredictor:
    def __init__(self, model_path='models/surrogate_weights.pth',
                 scaler_x_path='models/scaler_x.pkl',
                 scaler_y_path='models/scaler_y.pkl',
                 n_outputs=4, size=64):
        
        self.scaler_x = joblib.load(scaler_x_path)
        self.scaler_y = joblib.load(scaler_y_path)

        self.model = SurrogateModel(n_outputs=n_outputs, size=size)
        self.model.load_state_dict(torch.load(model_path))
        self.model.eval()

    def predict_without_scaling(self, raw_input):
        if np.array(raw_input).ndim == 1:
            raw_input = np.array(raw_input).reshape(1, -1)

        x_tensor = torch.tensor(raw_input, dtype=torch.float32)
        with torch.no_grad():
            y = self.model(x_tensor).numpy()

        return y.flatten() if y.shape[0] == 1 else y
    
    def predict(self, raw_input, scale_input=True, unscale_output=True):

        if np.array(raw_input).ndim == 1:
            raw_input = np.array(raw_input).reshape(1, -1)

        if scale_input==True:
            x_scaled = self.scaler_x.transform(raw_input)
        x_tensor = torch.tensor(x_scaled, dtype=torch.float32)

        with torch.no_grad():
            y_scaled = self.model(x_tensor).numpy()

        if unscale_output==False:
            return y_scaled.flatten() if y_scaled.shape[0] == 1 else y_scaled
        else:
            y_final = self.scaler_y.inverse_transform(y_scaled)
            return y_final.flatten() if y_final.shape[0] == 1 else y_final