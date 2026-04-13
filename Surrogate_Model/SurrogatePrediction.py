import torch 
import joblib
import numpy as np
from SurrogateModel import SurrogateModel

class SurrogatePredictor:
    def __init__(self, model_path='models/surrogate_weights.pth',
                 scaler_x_path='models/scaler_x.pkl',
                 scaler_y_path='models/scaler_y.pkl'):
        
        self.scaler_x = joblib.load(scaler_x_path)
        self.scaler_y = joblib.load(scaler_y_path)

        self.model = SurrogateModel()
        self.model.load_state_dict(torch.load(model_path))
        self.model.eval()

    def predict(self, raw_input):

        if np.array(raw_input).ndim == 1:
            raw_input = np.array(raw_input).reshape(1, -1)

        x_scaled = self.scaler_x.transform(raw_input)
        x_tensor = torch.tensor(x_scaled, dtype=torch.float32)

        with torch.no_grad():
            y_scaled = self.model(x_tensor).numpy()

        y_final = self.scaler_y.inverse_transform(y_scaled)

        return y_final.flatten() if y_final.shape[0] == 1 else y_final