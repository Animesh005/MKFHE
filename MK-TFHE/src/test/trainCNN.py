import torch
import torch.nn as nn
import torch.utils.data as data_utils
import torchvision
import torchvision.transforms as transforms
from torchvision import datasets
from torchvision.transforms import ToTensor
from torch.utils.data.sampler import SubsetRandomSampler
from torch.autograd import Variable
import time
import numpy as np
import os.path

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

cuda = torch.cuda.is_available()

torch.manual_seed(1)

if cuda:
    torch.cuda.manual_seed(1)

train_data = datasets.MNIST(
    root = 'data',
    train = True,                         
    transform = ToTensor(), 
    download = True,            
)
test_data = datasets.MNIST(
    root = 'data', 
    train = False, 
    transform = ToTensor(),
    download = True
)

from torch.utils.data import DataLoader
loaders = {
    'train' : torch.utils.data.DataLoader(train_data, 
                                          batch_size=100, 
                                          shuffle=True, 
                                          num_workers=1),
    
    'test'  : torch.utils.data.DataLoader(test_data, 
                                          batch_size=100, 
                                          shuffle=True, 
                                          num_workers=1),
}

# CNN model
class CNN(nn.Module):
    def __init__(self):
        super(CNN, self).__init__()
        self.conv = nn.Sequential(         
            nn.Conv2d(
                in_channels=1,              
                out_channels=5,            
                kernel_size=4,              
                stride=2,                   
                padding=0,                  
            ),                              
            nn.ReLU(),                      
            nn.MaxPool2d(kernel_size=4),    
        )

        # fully connected layer, output 10 classes
        self.fc = nn.Linear(45, 10)

    def forward(self, x):
        x = self.conv(x)
        x = x.view(x.size(0), -1)       
        output = self.fc(x)
        return output



num_epochs = 20

# training function
def train(model, optimizer, loaders, loss_fun):
    average_time = 0
    total = 0
    correct = 0

    model.train()
    # Train the model
    total_step = len(loaders['train'])

    for epoch in range(num_epochs):
      for i, (images, labels) in enumerate(loaders['train']):
          batch_time = time.time()
          images = Variable(images)
          labels = Variable(labels)

          if cuda:
              images, labels = images.cuda(), labels.cuda()

          optimizer.zero_grad()
          outputs = model(images)
          loss = loss_fun(outputs, labels)

          if cuda:
              loss.cpu()

          loss.backward()
          optimizer.step()

          batch_time = time.time() - batch_time
          average_time += batch_time

          _, predicted = torch.max(outputs.data, 1)
          total += labels.size(0)
          correct += (predicted == labels).sum().item()

          if (i + 1) % 10 == 0:
              print('Epoch: [%d/%d], Step: [%d/%d], Loss: %.4f, Accuracy: %.4f, Batch time: %f'
                    % (epoch + 1,
                      num_epochs,
                      i + 1,
                      total_step,
                      loss.item(),
                      correct / total,
                      average_time / 10))


model = CNN()
if cuda:
    model.cuda()

criterion = nn.CrossEntropyLoss()
if cuda:
    criterion.cuda()

optimizer = torch.optim.RMSprop(model.parameters(), lr=0.01)

# train the CNN model
train(model, optimizer, loaders, criterion)

# testing function
def test(model, loaders):
    model.eval()

    total = 0
    correct = 0
    for i, (data, labels) in enumerate(loaders['test']):
        data, labels = Variable(data), Variable(labels)
        if cuda:
            data, labels = data.cuda(), labels.cuda()

        data = data.squeeze(0)
        labels = labels.squeeze(0)

        outputs = model(data)
        if cuda:
            outputs.cpu()

        _, predicted = torch.max(outputs.data, 1)
        total += labels.size(0)
        correct += (predicted == labels).sum().item()

    return correct / total


# test the CNN model
test_acc = test(model, loaders)
print('=> Test set: Accuracy: {:.2f}%'.format(test_acc * 100))

# rounding the values of trained CNN parameters
for name, param in model.named_parameters():
    if 'weight' in name:
        param.data = torch.round(param.data)
    if 'bias' in name:
        param.data = torch.round(param.data)


# CNN inference function in clear
def test_inference(model, loaders):
    model.eval()

    total = 0
    correct = 0
    for i, (data, labels) in enumerate(loaders['test']):
        data = torch.round(data)
        data, labels = Variable(data), Variable(labels)
        if cuda:
            data, labels = data.cuda(), labels.cuda()

        data = data.squeeze(0)
        labels = labels.squeeze(0)

        outputs = model(data)
        if cuda:
            outputs.cpu()

        _, predicted = torch.max(outputs.data, 1)
        total += labels.size(0)
        correct += (predicted == labels).sum().item()

    return correct / total

# Inferencing on rounded parameters of CNN in clear
test_acc = test_inference(model, loaders)
print('=> Test set: Accuracy: {:.2f}%'.format(test_acc * 100))


# saving the parameters for oblivious CNN inference

weight_conv = model.conv[0].weight.cpu().data.numpy().squeeze(1).astype(int)
bias_conv = model.conv[0].bias.cpu().data.numpy().astype(int)

weight_fc = model.fc.weight.cpu().data.numpy().astype(int)
bias_fc = model.fc.bias.cpu().data.numpy().astype(int)

# Save the NumPy array to a text file
np.savetxt("model_params/cnn_weight.txt", weight_conv.reshape(-1, weight_conv.shape[-1]), delimiter=" ", fmt="%d")
np.savetxt("model_params/cnn_bias.txt", bias_conv.reshape(-1, bias_conv.shape[-1]), delimiter=" ", fmt="%d")

np.savetxt("model_params/fc_weight.txt", weight_fc, delimiter=" ", fmt="%d")
np.savetxt("model_params/fc_bias.txt", bias_fc, delimiter=" ", fmt="%d")

