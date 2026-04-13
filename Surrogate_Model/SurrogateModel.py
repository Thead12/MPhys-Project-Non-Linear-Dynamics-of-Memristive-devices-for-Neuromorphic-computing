
import torch.nn as nn

class SurrogateModel(nn.Module):
    def __init__(self):
        super(SurrogateModel, self).__init__()
        self.size = 256
        self.model = nn.Sequential(
            nn.Linear(2,self.size),
            nn.Tanh(),
            #nn.Dropout(p=0.1),
            nn.Linear(self.size,self.size),
            nn.Tanh(),
            nn.Linear(self.size, 4)
        )

    def forward(self, x):
        return self.model(x)
    