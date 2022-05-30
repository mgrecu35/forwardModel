import matplotlib
import matplotlib.pyplot as plt
from torch.nn.modules.module import _addindent
import torch
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np

                      
def torch_summarize(model, show_weights=True, show_parameters=True):
    """Summarizes torch model by showing trainable parameters and weights."""
    tmpstr = model.__class__.__name__ + ' (\n'
    for key, module in model._modules.items():
        # if it contains layers let call it recursively to get params and weights
        if type(module) in [
            torch.nn.modules.container.Container,
            torch.nn.modules.container.Sequential
        ]:
            modstr = torch_summarize(module)
        else:
            modstr = module.__repr__()
        modstr = _addindent(modstr, 2)

        params = sum([np.prod(p.size()) for p in module.parameters()])
        weights = tuple([tuple(p.size()) for p in module.parameters()])

        tmpstr += '  (' + key + '): ' + modstr 
        if show_weights:
            tmpstr += ', weights={}'.format(weights)
        if show_parameters:
            tmpstr +=  ', parameters={}'.format(params)
        tmpstr += '\n'   

    tmpstr = tmpstr + ')'
    return tmpstr

class Sequence(nn.Module):
    def __init__(self,ninp,nh,nout,ipia):
        super(Sequence, self).__init__()
        self.lstm_pia = nn.LSTMCell(ninp, nh)
        self.linear_pia = nn.Linear(nh, nout)
        self.lstm1 = nn.LSTMCell(ninp+ipia, nh)
        self.lstm2 = nn.LSTMCell(nh, nh)
        self.linear = nn.Linear(nh, nout)
        self.nh=nh
        self.ninp=ninp
        self.nout=nout
        self.ipia=ipia
        print(ninp)
    #@torch.jit.script    
    def forward(self, input, future = 0):
        outputs = []
        h_t = torch.zeros(input.size(0), self.nh, dtype=torch.double)
        c_t = torch.zeros(input.size(0), self.nh, dtype=torch.double)
        h_t2 = torch.zeros(input.size(0), self.nh, dtype=torch.double)
        c_t2 = torch.zeros(input.size(0), self.nh, dtype=torch.double) #(ns,self.nh)
        if self.ipia==1:
            h_t_pia = torch.zeros(input.size(0), self.nh, dtype=torch.double)
            c_t_pia = torch.zeros(input.size(0), self.nh, dtype=torch.double)
            output_pias=[]
            for input_t in input[:,:-1].split(self.ninp, dim=1):
                h_t_pia, c_t_pia = self.lstm_pia(input_t, (h_t_pia, c_t_pia)) 
                output_pia = self.linear(h_t_pia)
                output_pias+=[output_pia]
            output_pias=torch.cat(output_pias,dim=1)
            output_pias=torch.abs(output_pias)
            lstm_pia=output_pias.sum(axis=1)
            output_pias=torch.multiply(output_pias,input[:,-1,None])

        if self.ipia==1:
            for i,input_t in enumerate(input[:,:-1].split(self.ninp, dim=1)):
                input_t=torch.cat([input_t,output_pias[:,i:i+1]],dim=1)
                h_t, c_t = self.lstm1(input_t, (h_t, c_t)) #(ns,self.nh)
                h_t2, c_t2 = self.lstm2(h_t, (h_t2, c_t2))
                output = self.linear(h_t2)
                outputs += [output]
            outputs+=[output_pia]
        else:
            for i,input_t in enumerate(input[:,:].split(self.ninp, dim=1)):
                h_t, c_t = self.lstm1(input_t, (h_t, c_t)) #(ns,self.nh)
                h_t2, c_t2 = self.lstm2(h_t, (h_t2, c_t2))
                output = self.linear(h_t2)
                outputs += [output]
        outputs = torch.cat(outputs, dim=1)
        return outputs

from torch.utils.data import DataLoader
from torch.utils.data.dataset import Dataset

class MyCustomDataset(Dataset):
    def __init__(self, input_, target):
        self.input=input_
        self.target=target
        
    def __getitem__(self, index):
        # stuff
        return (self.input[index,:], self.target[index,:])

    def __len__(self):
        return self.input.shape[0]

    
from netCDF4 import Dataset
import glob

files=sorted(glob.glob("wrf*nc"))
input_=[]
target=[]

np.random.seed(0)
torch.manual_seed(0)
                  
#idf=0
#ipia=0

from sklearn.preprocessing import StandardScaler
import sys, argparse

parser = argparse.ArgumentParser()
parser.add_argument('--idf', type=int, default=0, help='dual frequency?')
parser.add_argument('--ipia', type=int, default=0, help='srt pia?')
opt = parser.parse_args()
print(opt)
idf=opt.idf
ipia=opt.ipia
print(idf,ipia)
idf,ipia=0,1
#stop
for f in files[:2]:
    fh=Dataset(f)
    zKu=fh["zKu"][:]
    zKa=fh["zKa"][:]
    pRate=fh["pRate"][:]
    dpia=fh["piaKa"][:]-fh["piaKu"][:]
    piaKu=fh["piaKu"][:]
    piaKa=fh["piaKa"][:]
    nz=zKu.shape[1]
    for i in range(zKu.shape[0]):
        nt=nz
        if idf==1:
            nt+=nz
        if ipia==1:
            nt+=1
        x2=np.zeros((nt),float)
        t=np.zeros((nz+ipia),float)
        if idf==1:
            if ipia==1:
                x2[:-1:2]=zKu[i,::-1]
                x2[1:-1:2]=zKa[i,:]
                x2[-1]=dpia[i]
            else:
                x2[::2]=zKu[i,::-1]
                x2[1::2]=zKa[i,:]
        else:
            if ipia==1:
                x2[:-1]=zKu[i,::-1]
                x2[-1]=piaKu[i]
            else:
                x2[:]=zKu[i,::-1]
        input_.append(x2)
        if ipia==1:
            t[:-1]=pRate[i,::-1]
            t[-1]=1.0
            target.append(t)
        else:
            t[:]=pRate[i,::-1]
            target.append(t)

stdScaler_input=StandardScaler()
stdScaler_target=StandardScaler()

nt=len(input_)

r=np.random.rand(nt)

stdScaler_input.fit(input_)
stdScaler_target.fit(target)
a=np.nonzero(r<0.5)
b=np.nonzero(r>0.5)
input_n=stdScaler_input.transform(input_)
target_n=stdScaler_target.transform(target)

t_input_n=torch.tensor(input_n[a[0],:])
t_target_n=torch.tensor(target_n[a[0],:])

v_input_n=torch.tensor(input_n[b[0],:])
v_target_n=torch.tensor(target_n[b[0],:])
                      
training_data=MyCustomDataset(t_input_n,t_target_n)

train_dataloader = DataLoader(training_data, batch_size=64, shuffle=False)

np.random.seed(0)
torch.manual_seed(0)

lstm_model = Sequence(idf+1,25,1,ipia)
lstm_model.double()

criterion = nn.MSELoss()

optimizer = optim.Adam(lstm_model.parameters())#,lr=0.5, momentum=0.1)#, lr=0.8)

for i in range(30):
    print('STEP: ', i)
    loss_av=[]
    for x,y in train_dataloader:
        def closure():
            optimizer.zero_grad()
            out = lstm_model(x)
            loss = criterion(out, y)
            loss_av.append(loss.item())
            #print('loss:', )
            loss.backward()
            return loss
        optimizer.step(closure)
    print(np.mean(loss_av))
    


torch.save(lstm_model.state_dict(),"lstm_model_%2.2i_%2.2i.pt"%(idf,ipia))

y_=lstm_model(v_input_n)
y_=y_[:,-2].detach().numpy()
target2_=v_target_n[:,-2].detach().numpy()
print(np.corrcoef(y_,target2_))

