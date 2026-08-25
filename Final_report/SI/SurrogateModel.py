# ------ MLP architecture build with py torch ------------
import torch.nn as nn

class SurrogateModel(nn.Module):
    def __init__(self, n_outputs=4, size=64):
        super(SurrogateModel, self).__init__()
        self.size = size
        self.model = nn.Sequential(
            nn.Linear(2,self.size),
            nn.Tanh(),
            #nn.Dropout(p=0.1),
            nn.Linear(self.size,self.size),
            nn.Tanh(),
            nn.Linear(self.size, n_outputs),
            nn.ReLU()
        )
    def forward(self, x):
        return self.model(x)
    