# FECG_Signal_Processing



This MATLAB program is designed to analyze fetal heart rate and uterine contraction signals from MIT-BIH Arrhythmia Database (http://physionet.org/physiobank/database/) to gain more insights to design the hardware Cardiotocography. 

How to run:

"Fetal_Heart_Rate_Extraction.m" file contains code to simulate the heart rate information from MIT-BIH Database. Ope this file in MATLAB and click run, you should be able to see all the figures. Note: the goal is to be able to identify the Fetal Heart Rate in beats per minute.

"EMG_Signal_Processing.m" contains code to study the Term-Preterm EMG signal from MIT-BIH Database. It's being used to study how to construct our own sensor to adapt the real signal.


Progress:

So far, my senior design team has finished the hardware prototyping on Fetal Heart Rate Monitoring sensor and EMG sensor. We are currently building the proper amplification circuit and BPF for the incoming signals.


Recorded Data file:

"undisturbed.txt" and "wireSwing30s.txt" contain recorded


Note:

the toolboxes are from Peyre's website (http://nbviewer.jupyter.org/github/gpeyre/numerical-tours/blob/master/matlab/denoisingwav_1_wavelet_1d.ipynb)

I included them here to for the convenience of transporting the entire project.
