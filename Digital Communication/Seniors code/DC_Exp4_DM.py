#!/usr/bin/env python
# coding: utf-8

# ## 20ECE1012 Deep Shetye
# ## Experiment 4
# ## Delta Modulation

# In[41]:


import matplotlib.pyplot as plt
import numpy as np
from scipy import signal


# In[42]:


# Generating the Input Signal and its Sampled Version
fm=100
fs=50*fm
A_m=1
t=np.arange(0,0.001+2/fm,0.0001)
ts=np.arange(0,2/fs+2/fm,1/fs)


# In[43]:


y=A_m*np.sin(2*np.pi*fm*t)
ys=A_m*np.sin(2*np.pi*fm*ts)


# In[44]:


#Calculating the ideal step size (delta)
delta=A_m*2*np.pi*fm*(1/fs)
print(delta)


# In[45]:


#Generating values of error, quantized error and quantized input
e=np.zeros(len(ts))
eq=np.zeros(len(ts))
mq=np.zeros(len(ts))

for i in range(0,len(ts)):
    if i==0:
        e[i]=ys[i]
        eq[i]=delta*np.sign(e[i])
        mq[i]=eq[i]
    else:
        e[i]=ys[i]-mq[i-1]
        eq[i]=delta*np.sign(e[i])
        mq[i]=eq[i]+mq[i-1]


# In[46]:


#Obtaining the reconstructed signal
m_rec=[]
for i in range(0,len(ts)):
    m_rec.append(e[i-1]+mq[i-1])


# In[47]:


#Smoothing of the reconstructed signal using Butterworth Low Pass Filter
[b,a]=signal.butter(3,(1/fm)*0.005,fs=0.0005)
rec=signal.lfilter(b,a,mq)


# In[48]:


plt.plot(t,y)
plt.title("Original Signal m(t)")
plt.xlabel("Time")
plt.ylabel("Amplitude")


# In[49]:


plt.stem(ts,ys)
plt.title("Sampled Signal m[n]")
plt.xlabel("Time")
plt.ylabel("Amplitude")


# ### Ideal Step Size (No Slope Overloading)

# In[50]:


plt.plot(ts,ys)
plt.step(ts,mq)
plt.title("Quantized Waveform for Ideal Step Size")
plt.xlabel("Time")
plt.ylabel("Amplitude")


# In[51]:


plt.step(ts,eq)
plt.title("Encoded Output")
plt.xlabel("Time")
plt.ylabel("Amplitude")


# In[52]:


plt.plot(ts,m_rec)
plt.title("Reconstructed Signal")
plt.xlabel("Time")
plt.ylabel("Amplitude")


# In[53]:


plt.plot(ts,rec)
plt.title("Reconstructed Signal after LPF")
plt.xlabel("Time")
plt.ylabel("Amplitude")


# ### Greater than Ideal Step Size (Positive Slope Overloading)

# In[54]:


delta=0.7

for i in range(0,len(ts)):
    if i==0:
        e[i]=ys[i]
        eq[i]=delta*np.sign(e[i])
        mq[i]=eq[i]
    else:
        e[i]=ys[i]-mq[i-1]
        eq[i]=delta*np.sign(e[i])
        mq[i]=eq[i]+mq[i-1]
        
m_rec=[]
for i in range(0,len(ts)):
    m_rec.append(e[i-1]+mq[i-1])

[b,a]=signal.butter(3,(1/fm)*0.005,fs=0.0005)
rec=signal.lfilter(b,a,mq)


# In[55]:


plt.plot(ts,ys)
plt.step(ts,mq)
plt.title("Quantized Waveform for High Delta")
plt.xlabel("Time")
plt.ylabel("Amplitude")


# In[56]:


plt.step(ts,eq)
plt.title("Encoded Output")
plt.xlabel("Time")
plt.ylabel("Amplitude")


# In[58]:


plt.plot(ts,m_rec)
plt.title("Reconstructed Signal")
plt.xlabel("Time")
plt.ylabel("Amplitude")


# In[59]:


plt.plot(ts,rec)
plt.title("Reconstructed Signal after LPF")
plt.xlabel("Time")
plt.ylabel("Amplitude")


# ### Less than Ideal Step Size (Negative Slope Overloading)

# In[60]:


delta=0.05

for i in range(0,len(ts)):
    if i==0:
        e[i]=ys[i]
        eq[i]=delta*np.sign(e[i])
        mq[i]=eq[i]
    else:
        e[i]=ys[i]-mq[i-1]
        eq[i]=delta*np.sign(e[i])
        mq[i]=eq[i]+mq[i-1]
        
m_rec=[]
for i in range(0,len(ts)):
    m_rec.append(e[i-1]+mq[i-1])

[b,a]=signal.butter(3,(1/fm)*0.005,fs=0.0005)
rec=signal.lfilter(b,a,mq)


# In[61]:


plt.plot(ts,ys)
plt.step(ts,mq)
plt.title("Quantized Waveform for Low Delta")
plt.xlabel("Time")
plt.ylabel("Amplitude")


# In[62]:


plt.step(ts,eq)
plt.title("Encoded Output")
plt.xlabel("Time")
plt.ylabel("Amplitude")


# In[63]:


plt.plot(ts,m_rec)
plt.title("Reconstructed Signal")
plt.xlabel("Time")
plt.ylabel("Amplitude")


# In[64]:


plt.plot(ts,rec)
plt.title("Reconstructed Signal after LPF")
plt.xlabel("Time")
plt.ylabel("Amplitude")

