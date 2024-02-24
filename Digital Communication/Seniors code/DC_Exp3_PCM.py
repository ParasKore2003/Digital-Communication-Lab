#!/usr/bin/env python
# coding: utf-8

# ## 20ECE1001 Aaron Lobo
# ## Experiment 3
# ## Pulse Code Modulation

# In[1]:


import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import signal


# ### Uniform PCM

# In[2]:


pi=np.pi
fm=eval(input("Signal Frequency: "))
fs=eval(input("Sampling Frequency: "))
t=np.arange(0,0.001+2/fm,0.001)
ts=np.arange(0,1/fs+2/fm,1/fs)
A=1
y=A*np.sin(2*pi*fm*t)
ys=A*np.sin(2*pi*fm*ts)


# In[3]:


plt.plot(t,y)
plt.stem(ts,ys,"r")
plt.title("Original and Sampled Signal")
plt.ylabel("Amplitude")
plt.xlabel("Time")
plt.show


# In[6]:


L=eval(input("No. of Quantization Levels: "))
N=int(math.log2(L))
delta=2*A/L

levels=np.arange(-A,A+delta,delta)
code=np.arange(-A+delta/2,A,delta)
print(levels)
print(code)
m=len(levels)
b=m-1
q=len(ts)
print(q)
x_quant=np.zeros(q)
x_quant_stair=np.zeros(len(t))
print(len(x_quant_stair))

for i in range(q):
    for j in range(1,m):
        if ys[i]==A:
            x_quant[i]=code[b-1]
        elif ys[i]>=levels[j-1] and ys[i]<levels[j]:
            x_quant[i]=code[j-1]

x_quant_stair=np.repeat(x_quant[0:len(x_quant)-1],50)
x_quant_stair=np.append(x_quant_stair,x_quant[len(x_quant)-1])

print(len(x_quant_stair))
print(x_quant[len(x_quant)-1])


# In[7]:


plt.plot(t,y)
plt.step(t,x_quant_stair)
plt.title("Quantized Waveform")
plt.ylabel("Amplitude - Discrete")
plt.xlabel("Time")
plt.show


# In[8]:


code_len=len(code)
dec=np.zeros(code_len)
for i in range(code_len):
    dec[i]=i

def DectoBin(n):
    return "{0:b}".format(int(n))

binr=[]
for i in range(code_len):
    binr.append(DectoBin(dec[i]))
print(binr)

for i in range(code_len):
    var=binr[i]
    if len(var)<N:
        for j in range(N-len(var)):
            var="0"+var
        binr[i]=var
print(binr)


# In[9]:


x_pcm=[]
vec_x_pcm=""
for i in range(len(x_quant)):
    for j in range(code_len):
        if x_quant[i]==code[j]:
            index=j
            x_pcm.append(binr[index])
            vec_x_pcm+=x_pcm[i]

vec_x_pcm2=np.zeros(len(vec_x_pcm))
for i in range(len(vec_x_pcm)):
    var2=int(vec_x_pcm[i])
    vec_x_pcm2[i]=var2

print(len(vec_x_pcm))
print(len(vec_x_pcm2))

plt.step(np.arange(0,len(vec_x_pcm2)),vec_x_pcm2)
plt.title("Encoded PCM Wave")
plt.ylabel("Amplitude")
plt.xlabel("Time")
plt.show


# In[10]:


def BintoDec(n):
    return int(n,2)

Xrc=[]
pcm_len=len(vec_x_pcm2)
for i in range(0,pcm_len,N):
    var3=vec_x_pcm2[i:i+3]
    string_int=[str(int(int2)) for int2 in var3]
    str_of_ints="".join(string_int)
    var4=BintoDec(str_of_ints)
    Xrc.append(code[var4])
print(Xrc)

plt.plot(np.array(Xrc))
plt.stem(np.array(Xrc))
plt.show


# In[11]:


[num,den]=signal.butter(3,2*fm/(fs/2))
recons=signal.lfilter(num,den,Xrc)
plt.plot(ts,recons)
plt.title("Recovered Signal")
plt.xlabel("Time")
plt.ylabel("Ampltiude")


# In[12]:


print(len(ys))
print(len(x_quant))


# In[13]:


NUM=0
DEN=0
for i in range(len(ys)):
    NUM+=ys[i]**2
for i in range(len(x_quant)):
    DEN+=(ys[i]-x_quant[i])**2


# In[14]:


SQNR=NUM/DEN
print(SQNR)

SQNRdb=10*math.log10(SQNR)
print(SQNRdb)


# ### Non-Uniform PCM

# In[15]:


A=1
fsig=eval(input("Signal Frequency:"))
Tsig=1/fsig
t=np.arange(0,2*Tsig,0.01*Tsig)


# In[16]:


x=A*np.sin(2*pi*fsig*t)
plt.plot(t,x)
plt.title("Original Continuous Signal x(t)")
plt.xlabel("Time")
plt.ylabel("Amplitude")
temp_sig=[i for i in x]
x_org=[i for i in x]


# In[17]:


#μ-Law
u=255
x=np.sign(temp_sig)*np.log(1+(u*np.abs(temp_sig)))/np.log(1+u)

plt.plot(t,x)
plt.title("Signal after application of μ-Law")
plt.xlabel("Time")
plt.ylabel("Amplitude")


# In[18]:


NyqMult=eval(input("Nyquist Mulitiplier:"))
fs=NyqMult*fsig
ts=np.arange(0,2*Tsig,1/fs)
xs=np.zeros(len(ts))
for i in range(len(ts)):
    xs[i]=x[i*5]

plt.plot(t,x)
plt.stem(ts,xs,"r")
plt.title("Sampled version of μ-Law Signal")
plt.ylabel("Amplitude")
plt.xlabel("Time")


# In[22]:


#Quantization
N=eval(input("Number of Bits:"))
L=2**N
delta=2*A/L

levels=np.arange(-A+delta/2,A,delta)
print(levels)

x_qnt=np.zeros(len(xs))
for i in range(len(xs)):
    index=math.floor((xs[i]-(-A))/delta)
    if index==L:
        index-=1
    x_qnt[i]=levels[index]

print(x_quant)

plt.plot(t,x)
plt.step(ts,x_qnt)
plt.title("Quantized Waveform")
plt.xlabel("Time")
plt.ylabel("Amplitude")


# In[23]:


code=[]
for i in range(L):
    num=i
    code_word=""
    for j in range(N):
        if num&(1<<j)==0:
            code_word+="0"
        else:
            code_word+="1"
    code.append(code_word)


# In[24]:


#Encoding
pcm=""
for i in range(len(x_qnt)):
    index=math.floor((x_qnt[i]-(-A))/delta)
    pcm+=code[index]

print(pcm)
pcm_sig=[int(i) for i in list(pcm)]
print(pcm_sig)

plt.step(np.arange(0,len(pcm_sig)),pcm_sig)
plt.title("Encoded PCM Wave")
plt.ylabel("Amplitude")
plt.xlabel("Time")


# In[25]:


#Receiver Reconstruction
code_dict={}
for i in range(len(code)):
    code_dict[code[i]]=levels[i]
x_rec=[]
for i in range(0,len(pcm),N):
    x_rec.append(code_dict[pcm[i:i+N]])
print(x_rec)
plt.plot(x_rec)
plt.stem(x_rec)
plt.title("Receiver Decoded Reconstructed Signal")
plt.xlabel("Time")
plt.ylabel("Amplitude")
temp_sig=[i for i in x_rec]


# In[26]:


#Inverse μ-Law
x_rec=np.sign(temp_sig)*((1+u)**np.abs(temp_sig)-1)/u
plt.plot(x_rec)
plt.stem(x_rec)


# In[27]:


[num,den]=signal.butter(3,2*fsig/(fs/2))
x_rec=signal.lfilter(num,den,x_rec)
plt.plot(ts,x_rec)
plt.title("Reconstructed Signal")
plt.xlabel("Time")
plt.ylabel("Amplitude")


# In[28]:


NUM=0
DEN=0
for i in range(len(xs)):
    NUM+=xs[i]**2
for i in range(len(x_qnt)):
    DEN+=(xs[i]-x_qnt[i])**2


# In[29]:


SQNR=NUM/DEN
print(SQNR)

SQNRdb=10*math.log10(SQNR)
print(SQNRdb)


# In[31]:


power=(np.linalg.norm(x_org,1)**2)/len(x_org)
sqnr=power+6.02*N+4.77
print(sqnr)

sqnr_db=10*math.log10(sqnr)
print(sqnr_db)


# In[ ]:




