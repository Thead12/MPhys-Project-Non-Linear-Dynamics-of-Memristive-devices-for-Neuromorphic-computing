import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

import torch
import torch.nn as nn
import torch.optim as optim

import pandas as pd

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from torch.utils.data import Dataset, DataLoader

import SurrogateModel

class SimpleDataset(Dataset):
    def __init__(self, X, y):
        self.X = torch.tensor(X, dtype=torch.float32)
        self.y = torch.tensor(y, dtype=torch.float32)

    def __len__(self):
        return len(self.X)
    
    def __getitem__(self, idx):
        return self.X[idx], self.y[idx]
    

df = pd.read_csv('LH_samples/latin_sample_100.csv')

X = df.iloc[:, :2].values


# V_DC,R_ext,is_phase_lock,width,CV1_at_target_freq,CV2_at_tartget_freq,centre_of_phase_locking,left_depth,righ_depth,depth
# 3 = width, 4 = CV1_at_target_freq, 5 = CV2_at_target_freq, 6 = centre_of_phase_locking, 9 = depth 
y = df.iloc[:, [3, 5, 6, 9]].values


scaler_x = StandardScaler()
scaler_y = StandardScaler()

X_scaled = scaler_x.fit_transform(X)
y_scaled = scaler_y.fit_transform(y)

X_train, X_val, y_train, y_val = train_test_split(X_scaled, y_scaled, test_size=0.2, random_state=42)



train_dataset = SimpleDataset(X_train, y_train)
val_dataset = SimpleDataset(X_val, y_val)

train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
val_x_tensor = torch.tensor(X_val, dtype=torch.float32)
val_y_tensor = torch.tensor(y_val, dtype=torch.float32)


model = SurrogateModel.SurrogateModel()

criterion = nn.MSELoss()
optimiser = optim.Adam(model.parameters(), lr=0.001, weight_decay=1e-4)
scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimiser, patience=150, factor=0.5)

train_losses = []
val_losses = []

epochs = 1600

# training Loop
for epoch in range(epochs):
    model.train()
    epoch_loss = 0
    for batch_x, batch_y in train_loader:

        # forward pass
        predictions = model(batch_x)

        # calculate loss
        loss = criterion(predictions, batch_y)

        # backward pass:
        optimiser.zero_grad()
        loss.backward()

        # update weights
        optimiser.step()

        epoch_loss += loss.item()

    model.eval()
    with torch.no_grad():
        v_loss = criterion(model(val_x_tensor), val_y_tensor)

    scheduler.step(v_loss)
    
    if (epoch + 1) % 20 == 0:
        train_losses.append(epoch_loss/len(train_loader))
        val_losses.append(v_loss.item())
        print(f"Epoch {epoch+1} | Train: {train_losses[-1]:.4f} | Val: {v_loss:.4f}")

import matplotlib.pyplot as plt
plt.plot(train_losses, label='Train')
plt.plot(val_losses, label='Val')
plt.legend()
plt.show()

torch.save(model.state_dict(), 'models/surrogate_weights.pth')

import joblib
joblib.dump(scaler_x, 'models/scaler_x.pkl')
joblib.dump(scaler_y, 'models/scaler_y.pkl')