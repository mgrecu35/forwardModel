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
    def __init__(self,ninp,nh,nout):
        super(Sequence, self).__init__()
        self.lstm_pia = nn.LSTMCell(ninp, nh)
        self.linear_pia = nn.Linear(nh, nout)
        self.lstm1 = nn.LSTMCell(ninp+1, nh)
        self.lstm2 = nn.LSTMCell(nh, nh)
        self.linear = nn.Linear(nh, nout)
        self.nh=nh
        self.ninp=ninp
        self.nout=nout
    #@torch.jit.script    
    def forward(self, input, future = 0):
        outputs = []
        h_t = torch.zeros(input.size(0), self.nh, dtype=torch.double)
        c_t = torch.zeros(input.size(0), self.nh, dtype=torch.double)
        h_t2 = torch.zeros(input.size(0), self.nh, dtype=torch.double)
        c_t2 = torch.zeros(input.size(0), self.nh, dtype=torch.double) #(ns,self.nh)
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
        
        for i,input_t in enumerate(input[:,:-1].split(self.ninp, dim=1)):
            input_t=torch.cat([input_t,output_pias[:,i:i+1]],dim=1)
            h_t, c_t = self.lstm1(input_t, (h_t, c_t)) #(ns,self.nh)
            h_t2, c_t2 = self.lstm2(h_t, (h_t2, c_t2))
            output = self.linear(h_t2)
            outputs += [output]
        outputs+=[output_pia]
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
        return len(input_)

    
from netCDF4 import Dataset
import glob

files=sorted(glob.glob("wrf*nc"))
input_=[]
target=[]
for f in files[:2]:
    fh=Dataset(f)
    zKu=fh["zKu"][:]
    zKa=fh["zKa"][:]
    pRate=fh["pRate"][:]
    dpia=fh["piaKa"][:]-fh["piaKu"][:]
    nz=zKu.shape[1]
    for i in range(zKu.shape[0]):
        x2=np.zeros((2*nz+1),float)
        t=np.zeros((nz+1),float)
        x2[:-1:2]=zKu[i,::-1]
        x2[1:-1:2]=zKa[i,:]
        x2[-1]=dpia[i]
        input_.append(x2)
        t[:-1]=pRate[i,::-1]
        t[-1]=1.0
        target.append(t)

from sklearn.preprocessing import StandardScaler

stdScaler_input=StandardScaler()
stdScaler_target=StandardScaler()

stdScaler_input.fit(input_)
stdScaler_target.fit(target)

input_n=torch.tensor(stdScaler_input.transform(input_))
target_n=torch.tensor(stdScaler_target.transform(target))

training_data=MyCustomDataset(input_n,target_n)

train_dataloader = DataLoader(training_data, batch_size=64, shuffle=False)

np.random.seed(0)
torch.manual_seed(0)

lstm_model = Sequence(2,25,1)
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
    
#test_dataloader = DataLoader(test_data, batch_size=128, shuffle=True)
#data = torch.load('traindata.pt')
#input = torch.from_numpy(data[3:, :-1])
#target = torch.from_numpy(data[3:, 1:])
#


#stop
if __name__ == '__main2__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--steps', type=int, default=15, help='steps to run')
    opt = parser.parse_args()
    # set random seed to 0
    np.random.seed(0)
    torch.manual_seed(0)
    # load data and make training set
    data = torch.load('traindata.pt')
    input = torch.from_numpy(data[3:, :-1])
    target = torch.from_numpy(data[3:, 1:])
    test_input = torch.from_numpy(data[:3, :-1])
    test_target = torch.from_numpy(data[:3, 1:])
    # build the model



    #begin to train
    for i in range(10):
        print('STEP: ', i)
        for x,y in train_dataloader:
            print(x.shape)
            def closure():
                optimizer.zero_grad()
                out = seq(x)
                loss = criterion(out, y)
                print('loss:', loss.item())
                loss.backward()
                return loss
            optimizer.step(closure)
        # begin to predict, no need to track gradient here
        with torch.no_grad():
            future = 1000
            pred = seq(test_input, future=future)
            loss = criterion(pred[:, :-future], test_target)
            print('test loss:', loss.item())
            y = pred.detach().numpy()
        # draw the result
        plt.figure(figsize=(30,10))
        plt.title('Predict future values for time sequences\n(Dashlines are predicted values)', fontsize=30)
        plt.xlabel('x', fontsize=20)
        plt.ylabel('y', fontsize=20)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        def draw(yi, color):
            plt.plot(np.arange(input.size(1)), yi[:input.size(1)], color, linewidth = 2.0)
            plt.plot(np.arange(input.size(1), input.size(1) + future), yi[input.size(1):], color + ':', linewidth = 2.0)
        draw(y[0], 'r')
        draw(y[1], 'g')
        draw(y[2], 'b')
        plt.savefig('predict%d.pdf'%i)
        plt.close()


torch.save(lstm_model.state_dict(),"lstm_model_DPIA.pt")
